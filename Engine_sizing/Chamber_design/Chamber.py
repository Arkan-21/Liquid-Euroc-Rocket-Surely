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

    def optimal_OF(self):
        reac_names = [self.fuel, self.oxidiser, self.additives]
        T_reactant = np.array([301.15, 150, 301.15])  # Reactant temperatures (K)
        fuel_weights = np.array([0.8, 0.0, 0.2]) #80-20 ethanol water mix
        oxidant_weights = np.array([0.0, 1.0, 0.0]) #pure lox

        #reactants and products
        reac = cea.Mixture(reac_names)
        prod = cea.Mixture(reac_names, products_from_reactants=True)

        #initialize solve for equilibrium
        solver = cea.RocketSolver(prod, reactants=reac)
        solution = cea.RocketSolution(solver)        

        of_ratio = np.arange(0, 5, 0.001)
        T = []
        Isp = []

        pi_p = [30, 35.0, 40.0]   # Pressure ratio
        supar = [10.0, 20.0, 40.0]  # Supersonic area ratio
        ac_at = 1.58  # Area ratio chamber to throat

        for of in of_ratio:
            weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of)
            hc = reac.calc_property(cea.ENTHALPY, weights, T_reactant)/cea.R
            solver.solve(solution, weights, self.pc, pi_p, supar=supar, ac_at=ac_at, iac=False, hc=hc)

            T.append(solution.T)
            Isp.append(solution.Isp)

        # Plot T vs OF
        plt.figure(figsize=(12, 5))

        plt.subplot(1, 2, 1)
        plt.plot(of_ratio, T)
        plt.xlabel('OF Ratio')
        plt.ylabel('Temperature (K)')
        plt.title('Chamber Temperature vs OF Ratio')
        plt.grid(True)

        # Plot Isp vs OF
        plt.subplot(1, 2, 2)
        plt.plot(of_ratio, Isp)
        plt.xlabel('OF Ratio')
        plt.ylabel('Specific Impulse (s)')
        plt.title('Specific Impulse vs OF Ratio')
        plt.grid(True)

        plt.tight_layout()
        plt.show()

if __name__ == "__main__":
    # Example usage
    chamber = chamber_sizing(oxidiser="O2(L)", fuel="C2H5OH(L)", additives="H2O", chamber_pressure=30, nominal_thrust=6000)
    chamber.optimal_OF()