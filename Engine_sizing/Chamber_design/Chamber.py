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

        of_ratio = np.arange(0, 5, 0.01)
        pi_p = [10.0, 12.0, 15.0]   # Pressure ratio
        supar = [8.0, 12.0, 15.0]  # Supersonic area ratio
        ac_at = 1.58  # Area ratio chamber to throat
        num_combinations = len(pi_p) * len(supar)
        T = [[] for _ in range(num_combinations)]
        Isp = [[] for _ in range(num_combinations)]
        Mw = [[] for _ in range(num_combinations)]
        gamma_list = [[] for _ in range(num_combinations)]

        for of in of_ratio:
            weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of)
            hc = reac.calc_property(cea.ENTHALPY, weights, T_reactant)/cea.R
            solver.solve(solution, weights, self.pc, pi_p, supar=supar, ac_at=ac_at, iac=False, hc=hc)

            for i in range(num_combinations):
                T[i].append(solution.T[i])
                Isp[i].append(solution.Isp[i]/10)
                Mw[i].append(solution.MW[i])
                gamma_list[i].append(solution.cp[i]/solution.cv[i])

        # Calculate optimal O/F for each combination (max Isp)
        optimal_ofs = []
        optimal_temps = []
        optimal_isps = []
        optimal_mws = []
        optimal_gammas = []
        for i in range(num_combinations):
            max_T = max(T[i])
            idx = T[i].index(max_T)
            optimal_of = of_ratio[idx]
            optimal_ISP = Isp[i][idx]
            optimal_Mw = Mw[i][idx]
            optimal_gamma = gamma_list[i][idx]
            optimal_ofs.append(optimal_of)
            optimal_temps.append(max_T)
            optimal_isps.append(optimal_ISP)
            optimal_mws.append(optimal_Mw)
            optimal_gammas.append(optimal_gamma)

        # Plot T vs OF
        plt.figure(figsize=(12, 5))

        plt.subplot(1, 2, 1)
        for i in range(num_combinations):
            pi = pi_p[i // len(supar)]
            su = supar[i % len(supar)]
            label = f'pi_p={pi}, supar={su}'
            plt.plot(of_ratio, T[i], label=label)
        plt.xlabel('OF Ratio')
        plt.ylabel('Temperature (K)')
        plt.title('Chamber Temperature vs OF Ratio')
        plt.grid(True)

        # Plot Isp vs OF
        plt.subplot(1, 2, 2)
        for i in range(num_combinations):
            pi = pi_p[i // len(supar)]
            su = supar[i % len(supar)]
            label = f'pi_p={pi}, supar={su}'
            plt.plot(of_ratio, Isp[i], label=label)
        plt.xlabel('OF Ratio')
        plt.ylabel('Specific Impulse (s)')
        plt.title('Specific Impulse vs OF Ratio')
        plt.grid(True)

        plt.tight_layout()
        plt.legend(loc='lower left')
        plt.show()

        # Print individual optimal values
        print("Individual optimal values:")
        for i in range(num_combinations):
            pi = pi_p[i // len(supar)]
            su = supar[i % len(supar)]
            print(f"pi_p={pi}, supar={su}: Optimal OF={optimal_ofs[i]:.3f}, T={optimal_temps[i]:.1f} K, Isp={optimal_isps[i]:.1f} s")

        # Find the best combination (max Isp)
        best_i = optimal_isps.index(max(optimal_isps))
        self.optimal_of = optimal_ofs[best_i]
        self.optimal_T = optimal_temps[best_i]
        self.optimal_Isp = optimal_isps[best_i]
        self.optimal_Mw = optimal_mws[best_i]
        self.optimal_gamma = optimal_gammas[best_i]

        # Store all optimal values
        self.optimal_ofs = optimal_ofs
        self.optimal_temps = optimal_temps
        self.optimal_isps = optimal_isps
        self.optimal_mws = optimal_mws
        self.optimal_gammas = optimal_gammas

    def throat_area(self, of_ratios, isps, temps, mws, gammas):
        g0 = 9.81  # m/s^2

        for i in range(len(isps)):
            OF = of_ratios[i]
            Isp = isps[i]
            T = temps[i]
            Mw = mws[i]
            gamma = gammas[i]

            v_e = Isp * g0  # Exhaust velocity (m/s)
            m_dot = self.thrust / v_e  # Total mass flow rate (kg/s)

            m_dot_ox = m_dot * OF / (OF + 1)  # Oxidizer mass flow rate (kg/s)
            m_dot_fuel = m_dot / (OF + 1)  # Fuel mass flow rate (kg/s)

            # Calculate throat area
            R = 8314 / Mw  # Gas constant (J/kgK)
            p_c = self.pc * 1e5  # Chamber pressure (Pa)

            # Throat area calculation for choked flow
            term1 = (2 / (gamma + 1)) ** ((gamma + 1) / (2 * (gamma - 1)))
            term2 = np.sqrt(gamma * R * T) / p_c
            A_t = m_dot / (p_c * term2 * term1)  # Throat area (m^2)

            # Print results for this combination
            print(f"\nCombination {i+1}:")
            print(f"  O/F: {OF:.3f}")
            print(f"  Temperature: {T:.1f} K")
            print(f"  Specific Impulse: {Isp:.1f} s")
            print(f"  Molecular Weight: {Mw:.3f} kg/mol")
            print(f"  Specific Heat Ratio (gamma): {gamma:.3f}")
            print(f"  Exhaust velocity (v_e): {v_e:.2f} m/s")
            print(f"  Total mass flow rate (m_dot): {m_dot:.2f} kg/s")
            print(f"  Oxidizer mass flow rate: {m_dot_ox:.2f} kg/s")
            print(f"  Fuel mass flow rate: {m_dot_fuel:.2f} kg/s")
            print(f"  Throat area (A_t): {A_t:.6f} m^2")

    def geometry(self, characteristic_length, throat_area):
        self.characteristic_length = characteristic_length
        self.throat_area = throat_area

        contraction_ratio = np.arange(1.5, 20, 0.01)
        contraction_angle = np.arange(np.deg2rad(20), np.deg2rad(45), np.deg2rad(0.01))

        for cr in contraction_ratio:
            for alpha in contraction_angle:

                throat_radius = np.sqrt(self.throat_area/np.pi)

                L_chamber = (throat_radius*(np.sqrt(cr)-1)+1.5*throat_radius*(1/np.cos(alpha)-1))/(np.tan(alpha))
                A_chamber = self.throat_area*cr
                chamber_radius =np.sqrt(A_chamber/np.pi)

                chmaber_volume = 1/3*np.pi()*(throat_radius**2+throat_radius*chamber_radius+chamber_radius**2)*L_chamber
                chamber_area = 2*L_chamber*np.sqrt(np.pi()*cr*self.throat_area)+1/np.sin(alpha)*(cr-1)*self.throat_area
                





if __name__ == "__main__":
    # Example usage
    chamber = chamber_sizing(oxidiser="O2(L)", fuel="C2H5OH(L)", additives="H2O", chamber_pressure=30, nominal_thrust=6000)
    chamber.optimal_OF()
    chamber.throat_area(chamber.optimal_ofs, chamber.optimal_isps, chamber.optimal_temps, chamber.optimal_mws, chamber.optimal_gammas)
