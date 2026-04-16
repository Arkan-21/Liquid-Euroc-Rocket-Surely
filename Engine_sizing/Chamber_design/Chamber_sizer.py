import numpy as np
import cea
import pandas as pd
from scipy.interpolate import CubicSpline

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
            Mw = solution.MW[0]   * 1e-3    # Chamber Molecular Weight in kg/mol
            gamma = solution.cp[0] / solution.cv[0] 
            
            return Isp, gamma, Mw, T
            
        except Exception as e:
            if self.debug: print(f"CEA Error: {e}")
            raise

    def throat_area(self, m_dot, T, pc, gamma, Mw):
        # Calculate the nozzle flow parameter (Gamma)
        # Using the standard equation for sonic flow at the throat
        exp = (1+ gamma) / (1 - gamma )
        Gamma = np.sqrt(gamma * (0.5*(gamma + 1))**exp)
        
        # At = (m_dot * sqrt(R_spec * T)) / (Pc * Gamma)
        R_univ = 8.31446261815324 # J/(mol*K)
        R_spec = R_univ / Mw  # Specific gas constant for the mixture in J/(kg*K)
        throat_area = (m_dot * np.sqrt(R_spec * T)) / (pc * 1e5 * Gamma)
        self.At = throat_area
        self.R_spec = R_spec
    
    def calculate_cr(self, propellant_type='Liq_Liq'):
        """
        Returns the contraction ratio based on throat diameter and propellant type.
        Performs cubic spline interpolation on log-transformed data.
        
        Parameters:
        -----------
        propellant_type : str
            Either 'Liq_Liq' or 'Gas_Liq' (both for 34 bar chamber pressure)
        
        Returns:
        --------
        float : Interpolated contraction ratio
        """
        
        # Data from the CSV files (as provided)
        if propellant_type == 'Liq_Liq':
            # Data from Liq_Liq_34barPC.csv
            data = np.array([
                [0.2039908382403323, 12.531716288281148],
                [0.27404024572728114, 10.802640572601012],
                [0.376283004688613, 9.20218263630338],
                [0.5638987547578291, 7.519642517601275],
                [0.9223036875932982, 5.964958350041419],
                [1.3086365745107398, 5.051136589214079],
                [1.9187085772178882, 4.2268064288737595],
                [3.036945319878672, 3.495241239282929],
                [5.189230005845924, 2.8057296834486847],
                [7.282822830493199, 2.4767495887985858],
                [11.034040200615467, 2.1605281145202024],
                [18.245495243701868, 1.8846805323975047],
                [28.87909714229446, 1.713839479342371],
                [44.23506942114378, 1.6054569578289242]
            ])
        elif propellant_type == 'Gas_Liq':
            # Data from Gas_Liq_34barPC.csv
            data = np.array([
                [0.23258948572929233, 14.798880290379048],
                [0.33001614540460755, 12.383747943992933],
                [0.50549665836649, 9.940787150341539],
                [0.679081617919856, 8.620240295978416],
                [0.9426935038960705, 7.34311440238132],
                [1.3522715402070757, 6.144740419456792],
                [1.816634650678577, 5.3284652584387],
                [2.549555192913489, 4.539030077860336],
                [3.578172242439654, 3.959504885040513],
                [5.852400676896963, 3.2355433986649453],
                [7.862085049632657, 2.9075097347640173],
                [14.188781765858016, 2.447505347432082],
                [21.032134799193297, 2.2522391210399166],
                [32.56977383414044, 2.097315645165334],
                [44.23506942114378, 2.0000000000000004]
            ])
        else:
            raise ValueError("propellant_type must be either 'Liq_Liq' or 'Gas_Liq'")
        
        # Extract diameters and contraction ratios
        diameters = data[:, 0]
        contraction_ratios = data[:, 1]
        
        # Apply log transformation
        log_diameters = np.log(diameters)
        log_contraction_ratios = np.log(contraction_ratios)
        
        # Create interpolation function on log-log scale
        interp_func = CubicSpline(log_diameters, log_contraction_ratios, 
                                    extrapolate=True)
        
        # Interpolate for the requested throat diameter
        throa_D = 2 * np.sqrt(self.At / np.pi) *39.37008 # Convert throat area to diameter in inches
        log_result = interp_func(np.log(throa_D))
        
        # Transform back from log scale
        cr = np.exp(log_result)
        self.cr = 8
        

    def chamber_geometry(self, L_star=1.2, conv_angle=30):
            """
            L_star: Characteristic length (m)
            conv_angle: Convergence half-angle (degrees) - typical 20-45
            
            """
            self.L_star = L_star
            self.conv_angle = np.radians(conv_angle)

            
            # 1. Chamber Volume from L*
            self.Vc = self.At * self.L_star
            self.Ac = self.At * self.cr
            
            # 2. Chamber Length (Lc)
            # Using the formula: Vc = V_cyl + V_frustum
            # Lc calculation accounting for the convergent cone volume
            term_frustum = (1/3) * np.sqrt(self.At / np.pi) * (1/np.tan(self.conv_angle)) * (self.cr**(1/3) - 1)
            self.Lc = (self.Vc / (self.At * self.cr)) - term_frustum
            

            print(f"Geometry: Lc={self.Lc:.2f} m, At={self.At*1e4:.4f} cm^2, Ac={self.Ac*1e4:.4f} cm^2, CR={self.cr:.3f}")


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
        print(f" {'Throat Diameter (Dt)':<25} | {2*self.Rt*1000:.2f} mm")
        print(f" {'Chamber Diameter (Dc)':<25} | {2*self.Rc*1000:.2f} mm")
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

    # Design assumptions
    of = 1.49              # O/F ratio (initial guess)
    expansion_ratio = 30 # Pc/Pe (pressure ratio)
    cr_guess = 5         # Initial contraction ratio guess

    # 1. Initialize the sizing object
    engine = chamber_sizing(
        oxidiser=oxidizer,
        fuel=fuel,
        additives=additives,
        chamber_pressure=chamber_pressure,
        nominal_thrust=nominal_thrust
    )

    # 2. Run CEA
    Isp, gamma, Mw, T = engine.run_cea(of, cr_guess, expansion_ratio)

    engine.Isp = Isp
    engine.gamma = gamma
    engine.Mw = Mw
    engine.T = T
    engine.of = of

    # 3. Compute mass flow from thrust
    g0 = 9.80665
    m_dot = nominal_thrust / (Isp * g0)
    engine.m_dot = m_dot

    # 4. Compute throat area
    engine.throat_area(m_dot, T, chamber_pressure, gamma, Mw)

    # 5. Compute contraction ratio from empirical relation
    engine.calculate_cr()
    engine.cr_val = engine.cr

    # 6. Re-run CEA with updated CR (optional but better)
    Isp, gamma, Mw, T = engine.run_cea(of, engine.cr, expansion_ratio)

    engine.Isp = Isp
    engine.gamma = gamma
    engine.Mw = Mw
    engine.T = T

    # 7. Compute geometry
    engine.chamber_geometry()

    # Derived radii
    engine.Rt = np.sqrt(engine.At / np.pi)
    engine.Rc = np.sqrt(engine.Ac / np.pi)

    # 8. Final report
    engine.report()
    print(engine.R_spec)