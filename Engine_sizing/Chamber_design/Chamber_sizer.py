import numpy as np
import cea

class chamber_sizing:
    def __init__(self, oxidiser, fuel, additives, chamber_pressure, nominal_thrust):
        self.oxidiser = oxidiser
        self.fuel = fuel
        self.additives = additives
        self.pc = chamber_pressure  # bar
        self.thrust = nominal_thrust  # N
        self.debug = True

    def run_cea(self, of, cr, pi):
        reac_names = [self.fuel, self.oxidiser]
        if self.additives and self.additives.strip():
            reac_names.append(self.additives)
        
        reac = cea.Mixture(reac_names)
        prod = cea.Mixture(reac_names, products_from_reactants=True)
        
        # Temperatures (K)
        T_fuel, T_ox, T_add = 298.15, 90.0, 298.15
        T_reactant = np.array([T_fuel, T_ox, T_add])
        
        # Mass fractions for the mixture
        fuel_weights = np.array([0.8, 0.0, 0.2]) 
        oxidant_weights = np.array([0.0, 1.0, 0.0])
        
        weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of)
        hc = reac.calc_property(cea.ENTHALPY, weights, T_reactant) / cea.R
        
        solver = cea.RocketSolver(prod, reactants=reac)
        solution = cea.RocketSolution(solver)
        
        try:
            solver.solve(
                solution,
                weights,
                self.pc,
                [pi],           # Pressure ratio Pc/Pe
                ac_at=cr,       # Contraction ratio
                iac=True,       # Finite Area Combustor logic
                hc=hc
            )
            
            # Indexing: 0 is Chamber, -1 is Exit
            Isp = solution.Isp[-1] / 9.80665 
            T = solution.T[0]               # Combustion Temperature
            Mw = solution.MW[0]             # Chamber Molecular Weight
            gamma = solution.cp[0] / solution.cv[0] 
            
            return Isp, gamma, Mw, T
            
        except Exception as e:
            if self.debug: print(f"CEA Error: {e}")
            raise

    def throat_area(self, m_dot, T, pc, gamma, Mw):
        # Calculate the nozzle flow parameter (Gamma)
        # Using the standard equation for sonic flow at the throat
        exp = (gamma + 1) / (gamma - 1)
        Gamma = np.sqrt(gamma * (2 / (gamma + 1))**exp)
        
        # At = (m_dot * sqrt(R_spec * T)) / (Pc * Gamma)
        R_univ = 8314.46
        throat_area = (m_dot * np.sqrt((R_univ / Mw) * T)) / (pc * 1e5 * Gamma)
        return throat_area
    
    def cr(self, throat_area):
        # Empirical relationship for CR based on throat diameter in inches
        Dt_inches = 2 * np.sqrt(throat_area / np.pi) * 39.3701
        return 1.302 * Dt_inches**(-0.481)

    def iterate(self, of_ratio=1.5, cr_init=3.0, pi=10.0, tol=1e-4, max_iter=100, relaxation=0.5):
        print(f"--- Starting Sizing: O/F={of_ratio}, P_ratio={pi} ---")
        current_cr = cr_init
        self.of = of_ratio
        self.pi = pi
        
        for i in range(max_iter):
            isp, gamma, mw, temp = self.run_cea(self.of, current_cr, self.pi)
            m_dot = self.thrust / (isp * 9.80665)
            at = self.throat_area(m_dot, temp, self.pc, gamma, mw)
            new_cr_target = self.cr(at)
            
            if self.debug:
                print(f"Iter {i+1}: CR_in={current_cr:.3f}, CR_out={new_cr_target:.3f}, Isp={isp:.1f}s")
            
            if abs(new_cr_target - current_cr) < tol:
                self.cr_val = current_cr
                self.Isp, self.gamma, self.Mw, self.T = isp, gamma, mw, temp
                self.m_dot, self.At = m_dot, at
                self.Rt = np.sqrt(at / np.pi)
                self.Rc = self.Rt * np.sqrt(self.cr_val)
                print(f"Converged in {i+1} steps.")
                return True

            current_cr =new_cr_target

        return False

    def geometry(self, L_star=1.2, conv_angle=30, div_angle=15):
            """
            L_star: Characteristic length (m)
            conv_angle: Convergence half-angle (degrees) - typical 20-45
            div_angle: Divergence half-angle (degrees) - typical 12-18
            """
            self.L_star = L_star
            self.conv_angle = np.radians(conv_angle)
            self.div_angle = np.radians(div_angle)
            
            # 1. Chamber Volume from L*
            self.Vc = self.At * self.L_star
            
            # 2. Chamber Length (Lc)
            # Using the formula: Vc = V_cyl + V_frustum
            # Lc calculation accounting for the convergent cone volume
            term_frustum = (1/3) * self.Rt * (1/np.tan(self.conv_angle)) * (self.cr_val**(1/3) - 1)
            self.Lc = (self.Vc / (self.At * self.cr_val)) - term_frustum
            
            # 3. Nozzle Exit Geometry
            # Expansion Ratio epsilon = Ae / At = pi (derived from pressure ratio)
            # However, for this simplified sizing, we can approximate epsilon if not returned by CEA
            # In your iterate, we pass 'pi' (Pc/Pe). Let's assume conical nozzle for length:
            self.Ae = self.At * self.pi # Note: 'pi' in your code is Pc/Pe, not epsilon. 
                                        # Ideally, get epsilon from CEA solution.
                                        # Assuming for now self.pi is a placeholder for epsilon for geom:
            self.Re = np.sqrt(self.Ae / np.pi)
            self.Ln = (self.Re - self.Rt) / np.tan(self.div_angle)

            print(f"Geometry: Lc={self.Lc:.2f} m, At={self.At*1e4:.4f} cm^2, Ae={self.Ae*1e4:.4f} cm^2, Ln={self.Ln*1000:.2f} mm")


    def report(self):
        print("\n" + "="*50)
        print("             ENGINE DESIGN SUMMARY")
        print("="*50)
        print(f" {'PARAMETER':<25} | {'VALUE':<15}")
        print("-" * 50)
        print(f" {'Thrust':<25} | {self.thrust:.1f} N")
        print(f" {'Chamber Pressure':<25} | {self.pc:.1f} bar")
        print(f" {'O/F Ratio':<25} | {self.of:.3f}")
        print(f" {'Specific Impulse':<25} | {self.Isp:.2f} s")
        print(f" {'Mass Flow Rate':<25} | {self.m_dot:.4f} kg/s")
        print("-" * 50)
        print(f" {'Throat Radius (Rt)':<25} | {self.Rt*1000:.2f} mm")
        print(f" {'Chamber Radius (Rc)':<25} | {self.Rc*1000:.2f} mm")
        print(f" {'Chamber Length (Lc)':<25} | {self.Lc*1000:.2f} mm")
        print(f" {'Contraction Ratio':<25} | {self.cr_val:.3f}")
        print(f" {'L* (Characteristic)':<25} | {self.L_star:.3f} m")
        print("-" * 50)
        print(f" {'Chamber Temp (Tc)':<25} | {self.T:.1f} K")
        print(f" {'Gamma':<25} | {self.gamma:.3f}")
        print(f" {'Molecular Weight':<25} | {self.Mw:.3f} g/mol")
        print("="*50 + "\n")

if __name__ == "__main__":
    # Define design requirements
    oxidizer = "O2(L)"
    fuel = "C2H5OH(L)"
    additives = "H2O"
    chamber_pressure = 30.0  # bar
    nominal_thrust = 6000.0  # N

    # 1. Initialize the sizing object
    engine = chamber_sizing(
        oxidiser=oxidizer,
        fuel=fuel,
        additives=additives,
        chamber_pressure=chamber_pressure,
        nominal_thrust=nominal_thrust
    )

    # Disable internal debug prints for a cleaner sweep output if desired
    engine.debug = False

    # 2. Define the pressure ratio range
    # From 30 to 100 in steps of 10
    pr_range = np.arange(30, 101, 10)
    
    results = []

    print(f"{'PR':>5} | {'Isp (s)':>10} | {'At (cm2)':>10} | {'CR':>8} | {'Rc (mm)':>8}")
    print("-" * 50)

    for pr in pr_range:
        # Run the iteration for each pressure ratio
        # We use the previous CR as the next starting guess for faster convergence
        success = engine.iterate(
            of_ratio=1.5,
            cr_init=0.1, 
            pi=float(pr),
            tol=1e-4,
            relaxation=0.5
        )

        if success:
            # Store or print the results
            engine.geometry(L_star=0.6, conv_angle=30)
            
            print(f"{pr:5.0f} | {engine.Isp:10.2f} | {engine.At*1e4:10.4f} | {engine.cr_val:8.3f} | {engine.Rc*1000:8.2f}")
            results.append({
                'pr': pr,
                'isp': engine.Isp,
                'cr': engine.cr_val,
                'rc': engine.Rc
            })
        else:
            print(f"{pr:5.0f} | FAILED TO CONVERGE")

    print("-" * 50)