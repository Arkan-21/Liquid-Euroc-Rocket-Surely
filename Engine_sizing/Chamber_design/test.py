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

    # -------------------------------------------------
    # CEA
    # -------------------------------------------------
    def run_cea(self, of, cr, eps):
        # Define reactants
        reac_names = [self.fuel, self.oxidiser]
        if self.additives and self.additives.strip(): #if not empty and not nan
            reac_names.append(self.additives)
        
        if self.debug:
            print(f"    CEA reactants: {reac_names}")
            print(f"    O/F: {of}, CR: {cr}, ER: {eps}")
        
        # Create mixture objects
        reac = cea.Mixture(reac_names)
        prod = cea.Mixture(reac_names, products_from_reactants=True)
        
        # Set reactant temperatures
        T_fuel = 298.15   # K 
        T_ox = 120.15     # K
        if self.additives and self.additives.strip():
            T_add = T_fuel      # K only additives to fuel
        T_reactant = np.array([T_fuel, T_ox, T_add])
        fuel_weights = np.array([0.8, 0.0, 0.2]) #percentage of fuel additives
        oxidant_weights = np.array([0.0, 1.0, 0.0])
        

        
        # Convert O/F ratio to weights
        try:
            weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of)
            if self.debug:
                print(f"    Weights: {weights}")
        except Exception as e:
            print(f"    Weight calculation error: {e}")
            raise
        
        # Calculate enthalpy of reactants
        try:
            hc = reac.calc_property(cea.ENTHALPY, weights, T_reactant) / cea.R
            if self.debug:
                print(f"    hc: {hc:.3f}")
        except Exception as e:
            print(f"    Enthalpy calculation error: {e}")
            raise
        
        # Solve
        solver = cea.RocketSolver(prod, reactants=reac)
        solution = cea.RocketSolution(solver)
        
        # Pressure ratios for nozzle expansion (use a single array)
        pressure_ratios = [40.0]
        
        try:
            # The solve method might need different parameters
            solver.solve(
                solution,
                weights,
                self.pc,           # chamber pressure in bar
                pressure_ratios,   # pressure ratios
                supar=[eps],       # exit pressure ratio
                ac_at=cr,          # contraction ratio
                iac=True,
                hc=hc
            )
            
            # Debug: print all available solution attributes
            if self.debug:
                print(f"    Solution attributes: {dir(solution)}")
                print(f"    solution.Isp: {solution.Isp}")
                print(f"    solution.T: {solution.T}")
                print(f"    solution.MW: {solution.MW}")
                print(f"    solution.cp: {solution.cp}")
                print(f"    solution.cv: {solution.cv}")
            
            # Extract results - need to check the correct indices
            # The solution might have arrays for different pressure ratios
            Isp = solution.Isp[-1] / 9.80665  # Use last pressure ratio (chamber conditions)
            T = solution.T[-1]
            Mw = solution.MW[-1]
            gamma = solution.cp[-1] / solution.cv[-1]
            
            if self.debug:
                print(f"    Raw Isp: {solution.Isp[-1]}, Converted: {Isp:.2f}")
            
            # Validate results
            if Isp <= 0 or np.isnan(Isp) or np.isinf(Isp):
                raise ValueError(f"Invalid Isp: {Isp}")
            
            if self.debug:
                print(f"   CEA success: Isp={Isp:.2f}s, T={T:.1f}K, γ={gamma:.3f}, Mw={Mw:.3f}")
            
            return Isp, gamma, Mw, T
            
        except Exception as e:
            print(f"    CEA solver error: {e}")
            print(f"    Trying with simplified parameters...")
            
            # Try simplified approach
            try:
                solver.solve(
                    solution,
                    weights,
                    self.pc,
                    [eps],           # Only expansion ratio
                    ac_at=cr,
                    iac=True,
                    hc=hc
                )
                
                # Try different indexing
                Isp = solution.Isp[0] / 9.80665 if solution.Isp[0] > 0 else solution.Isp[-1] / 9.80665
                T = solution.T[0] if solution.T[0] > 0 else solution.T[-1]
                Mw = solution.MW[0] if solution.MW[0] > 0 else solution.MW[-1]
                gamma = solution.cp[0] / solution.cv[0] if solution.cp[0] > 0 else solution.cp[-1] / solution.cv[-1]
                
                if Isp <= 0 or np.isnan(Isp):
                    raise ValueError(f"Still invalid Isp: {Isp}")
                
                return Isp, gamma, Mw, T
                
            except Exception as e2:
                print(f"     Simplified approach also failed: {e2}")
                raise

    # -------------------------------------------------
    # GEOMETRY
    # -------------------------------------------------
    def compute_geometry(self, cr, At, alpha=np.deg2rad(30)):
        Rt = np.sqrt(At / np.pi)
        Rc = Rt * np.sqrt(cr)
        
        # More robust length calculation
        if np.tan(alpha) > 0:
            Lc = (Rt * (np.sqrt(cr) - 1) + 1.5 * Rt * (1 / np.cos(alpha) - 1)) / np.tan(alpha)
        else:
            Lc = 0.1  # Fallback
            
        Vc = (1 / 3) * np.pi * (Rt**2 + Rt * Rc + Rc**2) * Lc
        
        return Vc, Lc, Rt, Rc
    
    # -------------------------------------------------
    # MAIN ITERATION LOOP
    # -------------------------------------------------
    def iterate(self,
                of=1.9,
                cr_init=3.0,
                eps=10.0,
                Lstar_target=1.2,
                tol=1e-3,
                max_iter=100,
                relaxation=0.5):
        
        cr = cr_init
        At = 1e-4
        Lstar = 0.0
        
        print("\n" + "="*60)
        print("Starting chamber sizing iteration...")
        print("="*60)
        
        for i in range(max_iter):
            print(f"\n{'='*60}")
            print(f"Iteration {i+1}/{max_iter}")
            print(f"{'='*60}")
            
            # --- CEA Calculation ---
            try:
                Isp, gamma, Mw, T = self.run_cea(of, cr, eps)
            except Exception as e:
                print(f" CEA calculation failed at iteration {i+1}")
                print(f"   Error: {e}")
                
                if i == 0:
                    print("\n📝 Troubleshooting suggestions:")
                    print("   1. Check if CEA is properly installed: pip install cea")
                    print("   2. Verify propellant names in CEA database")
                    print("   3. Try alternative naming: 'O2' for oxidizer, 'C2H5OH' for fuel")
                    print("   4. Check O/F ratio range (LOX/ethanol: 1.0-2.5)")
                raise
            
            # --- Flow Properties ---
            g0 = 9.80665
            v_e = Isp * g0
            
            if v_e <= 0:
                raise ValueError(f"Invalid exhaust velocity: {v_e:.2f} m/s")
            
            m_dot = self.thrust / v_e
            
            # --- Throat Area using characteristic velocity ---
            R_gas = 8314.462618 / Mw
            pc_pa = self.pc * 1e5
            
            # C* = sqrt(γ*R*T) / (γ * sqrt((2/(γ+1))^((γ+1)/(γ-1))))
            gamma_term = (2/(gamma+1))**((gamma+1)/(2*(gamma-1)))
            c_star = np.sqrt(gamma * R_gas * T) / gamma_term
            
            At_new = m_dot * c_star / pc_pa
            
            # Validate At_new
            if np.isnan(At_new) or np.isinf(At_new) or At_new <= 0:
                print(f"  Invalid throat area: {At_new}")
                At_new = At
            
            # --- Geometry ---
            try:
                Vc, Lc, Rt, Rc = self.compute_geometry(cr, At_new)
                Lstar_new = Vc / At_new
            except Exception as e:
                print(f"  Geometry calculation error: {e}")
                Lstar_new = Lstar_target
                Rt = np.sqrt(At_new / np.pi)
                Rc = Rt * np.sqrt(cr)
                Lc = 0.1
            
            # --- Print Results ---
            print(f"\n📊 Results:")
            print(f"  O/F         = {of:.4f}")
            print(f"  CR          = {cr:.4f}")
            print(f"  ER          = {eps:.3f}")
            print(f"  Isp         = {Isp:.2f} s")
            print(f"  C*          = {c_star:.1f} m/s")
            print(f"  At          = {At_new:.6e} m²  ({At_new*1e4:.2f} cm²)")
            print(f"  Rt          = {Rt:.4f} m  ({Rt*1000:.1f} mm)")
            print(f"  Rc          = {Rc:.4f} m  ({Rc*1000:.1f} mm)")
            print(f"  Lc          = {Lc:.3f} m  ({Lc*1000:.1f} mm)")
            print(f"  L*          = {Lstar_new:.3f} m  (target: {Lstar_target:.3f})")
            print(f"  m_dot       = {m_dot:.4f} kg/s")
            print(f"  γ           = {gamma:.3f}")
            print(f"  T           = {T:.1f} K")
            print(f"  Mw          = {Mw:.3f} g/mol")
            
            # --- Convergence Check ---
            At_error = abs(At_new - At) / At if At > 0 else 1.0
            Lstar_error = abs(Lstar_new - Lstar_target) / Lstar_target if Lstar_target > 0 else 1.0
            
            print(f"\n📈 Errors:")
            print(f"  At error    = {At_error:.2e}")
            print(f"  L* error    = {Lstar_error:.2e}")
            
            if At_error < tol and Lstar_error < tol:
                print(f"\n✅ Converged in {i+1} iterations!")
                break
            
            # --- Update Variables ---
            if not np.isnan(At_new) and not np.isinf(At_new) and At_new > 0:
                At = relaxation * At_new + (1 - relaxation) * At
            
            if Lstar_new > 0 and not np.isnan(Lstar_new) and not np.isinf(Lstar_new):
                cr_correction = (Lstar_target / Lstar_new) ** (1/3)
                cr_new = cr * cr_correction
                cr = relaxation * cr_new + (1 - relaxation) * cr
                cr = np.clip(cr, 1.5, 15.0)
        
        # Store results
        self.At = At
        self.cr = cr
        self.eps = eps
        self.Isp = Isp
        self.gamma = gamma
        self.Mw = Mw
        self.T = T
        self.Lstar = Lstar_new
        self.Lc = Lc
        self.Rt = Rt
        self.Rc = Rc
        self.m_dot = m_dot
        self.of = of
        self.c_star = c_star
        
        return At, cr, eps
    
    # -------------------------------------------------
    # 📊 FINAL REPORT
    # -------------------------------------------------
    def report(self):
        print("\n" + "="*60)
        print("FINAL ENGINE DESIGN SUMMARY")
        print("="*60)
        
        print("\n📋 OPERATING CONDITIONS:")
        print(f"  Thrust              : {self.thrust:.1f} N")
        print(f"  Chamber Pressure    : {self.pc:.1f} bar")
        print(f"  O/F Ratio           : {self.of:.3f}")
        print(f"  Specific Impulse    : {self.Isp:.2f} s")
        print(f"  Mass flow           : {self.m_dot:.4f} kg/s")
        
        print("\n📐 GEOMETRY:")
        print(f"  Throat area (At)    : {self.At:.6e} m²  ({self.At*1e4:.2f} cm²)")
        print(f"  Throat radius       : {self.Rt:.4f} m  ({self.Rt*1000:.1f} mm)")
        print(f"  Chamber radius      : {self.Rc:.4f} m  ({self.Rc*1000:.1f} mm)")
        print(f"  Chamber length      : {self.Lc:.4f} m  ({self.Lc*1000:.1f} mm)")
        print(f"  Contraction ratio   : {self.cr:.3f}")
        print(f"  Expansion ratio     : {self.eps:.3f}")
        print(f"  Characteristic L*   : {self.Lstar:.3f} m")
        
        print("\n🔥 GAS PROPERTIES (Chamber):")
        print(f"  Gamma (γ)           : {self.gamma:.3f}")
        print(f"  Molecular weight    : {self.Mw:.3f} g/mol")
        print(f"  Temperature         : {self.T:.1f} K  ({self.T-273.15:.1f} °C)")
        print(f"  Characteristic C*   : {self.c_star:.1f} m/s")
        
        print("\n" + "="*60)


# -------------------------------------------------
# ▶ RUN WITH DEBUGGING
# -------------------------------------------------
if __name__ == "__main__":
    
    # Test with gas-phase propellants first (more stable)
    print("Testing with gas-phase propellants...")
    chamber_gas = chamber_sizing(
        oxidiser="O2",
        fuel="C2H5OH",
        additives="",
        chamber_pressure=30,
        nominal_thrust=6000
    )
    
    try:
        chamber_gas.iterate(
            of=1.2,
            cr_init=3.0,
            eps=8.0,
            Lstar_target=1.2,
            tol=5e-3,
            max_iter=5,
            relaxation=0.7
        )
        chamber_gas.report()
    except Exception as e:
        print(f"Gas-phase test failed: {e}")
        
        # If gas-phase works, try liquid
        print("\n" + "="*60)
        print("Testing with liquid-phase propellants...")
        print("="*60)
        
        chamber_liq = chamber_sizing(
            oxidiser="O2(L)",
            fuel="C2H5OH(L)",
            additives="H2O",
            chamber_pressure=30,
            nominal_thrust=6000
        )
        
        try:
            chamber_liq.iterate(
                of=1.7,
                cr_init=3.0,
                eps=8.0,
                Lstar_target=1.2,
                tol=5e-3,
                max_iter=60,
                relaxation=0.7
            )
            chamber_liq.report()
        except Exception as e2:
            print(f"Liquid-phase test also failed: {e2}")
            print("\n  CEA may not have liquid propellant data available")
            print("   Try using gas-phase propellants for initial design")