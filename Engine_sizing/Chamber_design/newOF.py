import numpy as np
import cea
from matplotlib import pyplot as plt

class OFplot:
    def __init__(self, oxidiser, fuel, additives, chamber_pressure, nominal_thrust):
        self.oxidiser = oxidiser
        self.fuel = fuel
        self.additives = additives
        self.pc = chamber_pressure
        self.thrust = nominal_thrust

    def optimal_OF(self):
        """ 
        Returns an arrays of optimal O/F ratios for specific pressure ratios.
        [30, 40, 50, 60, 70, 80, 90, 100] 

        for a chamber pressure of 30bar, that corresponds to max altitude (pressure ratio of 100) of 9km
        

        """


        # Using (L) suffix ensures CEA uses liquid-phase enthalpies for LOX/Ethanol
        reac_names = [self.fuel, self.oxidiser, f"{self.additives}(L)"]
        T_reactant = np.array([298.15, 90.17, 298.15]) # Standard temps for Eth, LOX, Water
        
        fuel_weights = np.array([0.8, 0.0, 0.2]) # 80/20 Ethanol/Water
        oxidant_weights = np.array([0.0, 1.0, 0.0]) # Pure LOX

        reac = cea.Mixture(reac_names)
        prod = cea.Mixture(reac_names, products_from_reactants=True)
        solver = cea.RocketSolver(prod, reactants=reac)
        solution = cea.RocketSolution(solver)        

        # Define the ranges you want to test
        of_ratios = np.arange(0, 7.0, 0.01)
        pi_list = [30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]   # Pressure ratios
        
        plt.figure(figsize=(10, 6))

        optimal_of = []
        cr = 1.0  # Nozzle contraction ratio (not used in this specific plot, but can be added for more complex analyses)
        


        # Nested Loops: Ensure we control the 'Pi' and 'Epsilon' for every solve
        for pi in pi_list:

            isp_results = []
            temp_results = []
                
            for of in of_ratios:
                weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of)
                hc = reac.calc_property(cea.ENTHALPY, weights, T_reactant)/cea.R
                    
                # Passing single values for pi and supar ensures index 2 is ALWAYS the exit
                solver.solve(solution, weights,  self.pc, [pi], cr=cr, hc=hc)
                    
                # solution.Isp[0]=Chamber, [1]=Throat, [2]=Exit
                isp_seconds = solution.Isp[2] / 9.80665
                isp_results.append(isp_seconds)
                temp_results.append(solution.T[0])

            # Find peaks for this specific geometry
            max_isp = max(isp_results)
            best_of = of_ratios[isp_results.index(max_isp)]
            optimal_of.append(best_of)
                
            label = f"Pi={pi}, (Max: {max_isp:.1f}s @ O/F {best_of:.2f})"
            plt.plot(of_ratios, isp_results, label=label)
            print(f"Geometry [Pi:{pi}] -> Best O/F: {best_of:.2f}, Isp: {max_isp:.1f}s, Temp: {temp_results[isp_results.index(max_isp)]:.1f}K")

        #plt.axvline(x=1.6, color='k', linestyle='--', alpha=0.3, label="Typical LOX/Eth Opt")
        plt.xlabel('Oxidizer/Fuel Ratio')
        plt.ylabel(' Isp (s)')
        plt.title(f'LOX/Ethanol Performance @ Pc={self.pc} bar')
        plt.legend(fontsize='small')
        plt.grid(True, which='both', linestyle=':', alpha=0.5)
        plt.show()

        return optimal_of

if __name__ == "__main__":
    # Use exact CEA nomenclature
    of_plot = OFplot(oxidiser="O2(L)", fuel="C2H5OH(L)", additives="H2O", chamber_pressure=30, nominal_thrust=6000)
    of_plot.optimal_OF()