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
        T_ox = 90     # K
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
                iac=True,          #infinite area combsustor
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


    def throat_area(self, m_dot, T, pc, gamma, Mw):

        Gamma = np.sqrt(gamma*(0.5*(1+gamma))**((gamma+1)/(gamma-1)))
        throat_area = m_dot * np.sqrt(8314.462618/Mw*T) / (pc*1e5*Gamma)
        return throat_area
    
    def cr(self, throat_area):
        Dt = 2 * np.sqrt(throat_area / np.pi) * 39.37008  # Convert m to inches
        cr = 1.302 * Dt**(-0.481)
        return cr

  
    # -------------------------------------------------
    # MAIN ITERATION LOOP
    # -------------------------------------------------
    def iterate(self,
                cr_init=3.0,
                eps=10.0,
                Lstar_target=1.2,
                tol=1e-3,
                max_iter=100,
                relaxation=0.5):
        


    
    # -------------------------------------------------
    #  FINAL REPORT
    # -------------------------------------------------
    def report(self):
        print("\n" + "="*60)
        print("FINAL ENGINE DESIGN SUMMARY")
        print("="*60)
        
        print("\n OPERATING CONDITIONS:")
        print(f"  Thrust              : {self.thrust:.1f} N")
        print(f"  Chamber Pressure    : {self.pc:.1f} bar")
        print(f"  O/F Ratio           : {self.of:.3f}")
        print(f"  Specific Impulse    : {self.Isp:.2f} s")
        print(f"  Mass flow           : {self.m_dot:.4f} kg/s")
        
        print("\n GEOMETRY:")
        print(f"  Throat area (At)    : {self.At:.6e} m²  ({self.At*1e4:.2f} cm²)")
        print(f"  Throat radius       : {self.Rt:.4f} m  ({self.Rt*1000:.1f} mm)")
        print(f"  Chamber radius      : {self.Rc:.4f} m  ({self.Rc*1000:.1f} mm)")
        print(f"  Chamber length      : {self.Lc:.4f} m  ({self.Lc*1000:.1f} mm)")
        print(f"  Contraction ratio   : {self.cr:.3f}")
        print(f"  Expansion ratio     : {self.eps:.3f}")
        print(f"  Characteristic L*   : {self.Lstar:.3f} m")
        
        print("\n GAS PROPERTIES (Chamber):")
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