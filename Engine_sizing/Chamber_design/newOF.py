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

        # --- Initialize CEA once ---
        self.reac_names = [self.fuel, self.oxidiser, f"{self.additives}(L)"]
        self.reac = cea.Mixture(self.reac_names)
        self.prod = cea.Mixture(self.reac_names, products_from_reactants=True)
        self.solver = cea.RocketSolver(self.prod, reactants=self.reac)
        self.solution = cea.RocketSolution(self.solver)

        self.T_reactant = np.array([298.15, 90.17, 298.15])

    # --------------------------------------------------
    # 1. CORE FUNCTION (single Pi, single additives %)
    # --------------------------------------------------
    def compute_of_curve(self, pi, of_ratios, additive_pct=0.2):
        """
        Compute Isp and temperature vs OF for ONE pressure ratio and ONE additive %.
        """

        fuel_weights = np.array([1 - additive_pct, 0.0, additive_pct])
        oxidant_weights = np.array([0.0, 1.0, 0.0])

        isp_results = []
        temp_results = []

        for of in of_ratios:
            weights = self.reac.of_ratio_to_weights(oxidant_weights, fuel_weights, of)
            hc = self.reac.calc_property(cea.ENTHALPY, weights, self.T_reactant) / cea.R

            self.solver.solve(self.solution, weights, self.pc, [pi], hc=hc)

            isp_seconds = self.solution.Isp[2] / 9.80665
            isp_results.append(isp_seconds)
            temp_results.append(self.solution.T[0])  # chamber temperature

        return np.array(isp_results), np.array(temp_results)

    # --------------------------------------------------
    # 2. SWEEP MULTIPLE PRESSURE RATIOS (Pi)
    # --------------------------------------------------
    def sweep_pressure_ratios(self, pi_list, of_ratios, additive_pct=0.2):
        results = {}
        optimal_of = []

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        for pi in pi_list:
            isp, temp = self.compute_of_curve(pi, of_ratios, additive_pct)

            max_isp = np.max(isp)
            best_of = of_ratios[np.argmax(isp)]
            optimal_of.append(best_of)

            label = f"Pi={pi}, Max: {max_isp:.1f}s @ OF {best_of:.2f}"

            # --- Isp plot ---
            ax1.plot(of_ratios, isp, label=label)

            # --- Temperature plot ---
            ax2.plot(of_ratios, temp, label=f"Pi={pi}")

            print(f"[Pi:{pi}] -> Best O/F: {best_of:.2f}, "
                  f"Isp: {max_isp:.1f}s, "
                  f"Temp: {temp[np.argmax(isp)]:.1f}K")

            results[pi] = isp

        # --- Formatting ---
        ax1.set_xlabel('O/F ratio')
        ax1.set_ylabel('Isp (s)')
        ax1.set_title(f'Isp vs OF @ Pc={self.pc} bar')
        ax1.grid(True, linestyle=':', alpha=0.5)
        ax1.legend(fontsize='small')

        ax2.set_xlabel('O/F ratio')
        ax2.set_ylabel('Chamber Temperature (K)')
        ax2.set_title(f'Temperature vs OF @ Pc={self.pc} bar')
        ax2.grid(True, linestyle=':', alpha=0.5)
        ax2.legend(fontsize='small')

        plt.tight_layout()
        plt.show()

        return results, optimal_of

    # --------------------------------------------------
    # 3. SWEEP ADDITIVES (FIXED Pi)
    # --------------------------------------------------
    def sweep_additives(self, pi, of_ratios, additive_list):
        results = {}

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        print("\n=== ADDITIVES SWEEP SUMMARY ===")

        for pct in additive_list:
            isp, temp = self.compute_of_curve(pi, of_ratios, pct)

            # --- Isp peak ---
            max_isp = np.max(isp)
            best_of_isp = of_ratios[np.argmax(isp)]

            # --- Temperature peak ---
            max_temp = np.max(temp)
            best_of_temp = of_ratios[np.argmax(temp)]

            label = f"{pct*100:.0f}%, Max: {max_isp:.1f}s @ OF {best_of_isp:.2f}"

            # --- Plot Isp ---
            ax1.plot(of_ratios, isp, label=label)

            # --- Plot Temperature ---
            ax2.plot(of_ratios, temp, label=f"{pct*100:.0f}%")

            # --- Terminal output ---
            print(
                f"[Add:{pct*100:.0f}%] "
                f"Isp max: {max_isp:.1f}s @ OF {best_of_isp:.2f} | "
                f"Temp max: {max_temp:.1f}K @ OF {best_of_temp:.2f}"
            )

            results[pct] = {
                "isp": isp,
                "temp": temp,
                "max_isp": max_isp,
                "of_max_isp": best_of_isp,
                "max_temp": max_temp,
                "of_max_temp": best_of_temp
            }

        # --- Formatting ---
        ax1.set_xlabel('O/F ratio')
        ax1.set_ylabel('Isp (s)')
        ax1.set_title(f'Isp vs OF @ Pi={pi}')
        ax1.grid(True, linestyle=':', alpha=0.5)
        ax1.legend()

        ax2.set_xlabel('O/F ratio')
        ax2.set_ylabel('Chamber Temperature (K)')
        ax2.set_title(f'Temperature vs OF @ Pi={pi}')
        ax2.grid(True, linestyle=':', alpha=0.5)
        ax2.legend()

        plt.tight_layout()
        plt.show()

        return results


# --------------------------------------------------
# MAIN
# --------------------------------------------------
if __name__ == "__main__":

    of_plot = OFplot(
        oxidiser="O2(L)",
        fuel="C2H5OH(L)",
        additives="H2O",
        chamber_pressure=30,
        nominal_thrust=6000
    )

    of_ratios = np.arange(0.5, 3.0, 0.05)

    # --------------------------------------------------
    # 1. SINGLE CASE (side-by-side plot)
    # --------------------------------------------------
    isp, temp = of_plot.compute_of_curve(pi=50, of_ratios=of_ratios)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    ax1.plot(of_ratios, isp)
    ax1.set_title("Isp vs OF (Pi=50)")
    ax1.set_xlabel("O/F")
    ax1.set_ylabel("Isp (s)")
    ax1.grid()

    ax2.plot(of_ratios, temp)
    ax2.set_title("Temperature vs OF (Pi=50)")
    ax2.set_xlabel("O/F")
    ax2.set_ylabel("Temperature (K)")
    ax2.grid()

    plt.tight_layout()
    plt.show()

    # --------------------------------------------------
    # 2. PRESSURE RATIO SWEEP
    # --------------------------------------------------
    pi_list = [30, 40, 50, 60, 70, 80, 90, 100]
    of_plot.sweep_pressure_ratios(pi_list, of_ratios)

    # --------------------------------------------------
    # 3. ADDITIVES SWEEP
    # --------------------------------------------------
    additive_list = [0.0, 0.1, 0.2, 0.3]
    of_plot.sweep_additives(pi=50, of_ratios=of_ratios, additive_list=additive_list)