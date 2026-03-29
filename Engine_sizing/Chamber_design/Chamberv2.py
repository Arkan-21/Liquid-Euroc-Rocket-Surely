import numpy as np
import cea
from matplotlib import pyplot as plt

class chamber_sizing:
    def __init__(self, oxidiser, fuel, additives,  chamber_pressure, nominal_thrust):
        self.oxidiser = oxidiser
        self.fuel = fuel
        self.additives = additives
        self.pc = chamber_pressure
        self.thrust = nominal_thrust
    # -------------------------------------------------
    # CEA
    # -------------------------------------------------
    def run_cea(self, of, cr, eps, T_fuel, T_ox):
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
        pressure_ratios = [30.0, 40.0, 50.0]
        
        try:
            # The solve method might need different parameters
            solver.solve(
                solution,
                weights,
                self.pc,           # chamber pressure in bar
                pressure_ratios,   # pressure ratios
                supar=[eps],       # supersonic area ratio
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
