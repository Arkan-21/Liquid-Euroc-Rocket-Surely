import sys
import os
import numpy as np
import cea
from astropy import units as u

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Design_overview import input_parameters as p


class ChamberSizer:
    """
    Chamber sizing tool for liquid rocket engines.
    Calculates basic performance characteristics and geometry based on:
    - Fuel and oxidiser types
    - Chamber pressure
    - Nominal thrust
    """
    
    def __init__(self, fuel, oxidiser, chamber_pressure, thrust):
        """
        Initialize the chamber sizer.
        
        Parameters
        ----------
        fuel : str
            Fuel type (e.g., "80Ethanol")
        oxidiser : str
            Oxidiser type (e.g., "LOX")
        chamber_pressure : astropy.units.Quantity
            Chamber pressure
        thrust : astropy.units.Quantity
            Nominal thrust
        """
        self.fuel = fuel
        self.oxidiser = oxidiser
        self.chamber_pressure = chamber_pressure.to(u.pascal)
        self.thrust = thrust.to(u.newton)
        
        # Initialize CEA mixture
        self.mixture = cea.Mixture([fuel], [oxidiser])
        
        # Store results
        self.results = {}
    
    def calculate_performance(self, of_ratio=2.0):
        """
        Calculate engine performance characteristics using CEA.
        
        Parameters
        ----------
        of_ratio : float, optional
            Oxidiser-to-fuel mass ratio (default: 2.0)
        
        Returns
        -------
        dict
            Dictionary containing performance characteristics
        """
        # Get CEA properties at chamber conditions
        cea_output = self.mixture.exhaust_at(
            self.chamber_pressure,
            of=of_ratio
        )
        
        # Extract key properties
        tc = cea_output["T"]  # Chamber temperature (K)
        gamma = cea_output["gamma"]  # Specific heat ratio
        molecular_weight = cea_output["M"]  # Molecular weight (g/mol)
        
        # Convert to standard units
        self.results['chamber_temp'] = tc * u.kelvin
        self.results['gamma'] = gamma
        self.results['molecular_weight'] = molecular_weight * u.gram / u.mol
        self.results['of_ratio'] = of_ratio
        
        # Calculate characteristic velocity (c*)
        # c* = sqrt(gamma * R * T_c / (gamma + 1)) where R is specific gas constant
        R_universal = 8.314 * u.joule / (u.mol * u.kelvin)
        R_specific = R_universal / (molecular_weight / 1000)  # Convert g/mol to kg/mol
        
        c_star = np.sqrt(
            gamma * R_specific * self.results['chamber_temp'] / (gamma + 1)
        ).to(u.meter / u.second)
        
        self.results['c_star'] = c_star
        
        # Assume nozzle expansion ratio (typical range 8-40)
        # For basic geometry, use a conservative value
        expansion_ratio = 8.0
        
        # Calculate Isp using thrust equation: Isp = (F / (m_dot * g0))
        # Also calculate from exit velocity for consistency
        g0 = 9.81 * u.meter / u.second**2
        
        # Exit velocity: Ve = c* * sqrt(2*gamma/(gamma-1) * (1 - (1/expansion_ratio)^((gamma-1)/gamma)))
        exponent = (gamma - 1) / gamma
        exit_velocity = c_star * np.sqrt(
            2 * gamma / (gamma - 1) * (1 - expansion_ratio**(-exponent))
        )
        
        self.results['exit_velocity'] = exit_velocity.to(u.meter / u.second)
        
        # Specific impulse
        isp = exit_velocity / g0
        self.results['isp'] = isp.to(u.second)
        
        # Mass flow rate from thrust equation
        mass_flow = self.thrust / exit_velocity
        self.results['mass_flow_rate'] = mass_flow.to(u.kg / u.second)
        
        return self.results
    
    def calculate_geometry(self):
        """
        Calculate chamber geometry based on performance characteristics.
        
        Returns
        -------
        dict
            Dictionary containing geometry characteristics
        """
        if not self.results:
            raise ValueError("Must run calculate_performance() first")
        
        m_dot = self.results['mass_flow_rate']
        c_star = self.results['c_star']
        
        # Critical mass flow rate: m_dot_crit = (P_c * A_t) / c*
        # Rearranged: A_t = (m_dot * c*) / P_c
        
        throat_area = (m_dot * c_star) / self.chamber_pressure
        throat_area = throat_area.to(u.centimeter**2)
        
        # Throat diameter
        throat_diameter = 2 * np.sqrt(throat_area / np.pi)
        throat_diameter = throat_diameter.to(u.centimeter)
        
        self.results['throat_area'] = throat_area
        self.results['throat_diameter'] = throat_diameter
        
        # Chamber diameter (typically 1.5-2.5 x throat diameter)
        # Using 2.0 as baseline
        chamber_diameter = 2.0 * throat_diameter
        self.results['chamber_diameter'] = chamber_diameter.to(u.centimeter)
        
        # Chamber length (typically 1-3 x chamber diameter)
        # Using 1.5 as baseline for compact design
        chamber_length = 1.5 * chamber_diameter
        self.results['chamber_length'] = chamber_length.to(u.centimeter)
        
        # Chamber volume
        chamber_radius = chamber_diameter / 2
        chamber_volume = np.pi * chamber_radius**2 * chamber_length
        self.results['chamber_volume'] = chamber_volume.to(u.centimeter**3)
        
        # Residence time (rough estimate)
        characteristic_length = chamber_volume / throat_area
        residence_time = characteristic_length / c_star
        self.results['characteristic_length'] = characteristic_length.to(u.centimeter)
        self.results['residence_time'] = residence_time.to(u.millisecond)
        
        return self.results
    
    def get_summary(self):
        """Get a formatted summary of all results."""
        summary = {
            "Performance": {
                "Chamber Temperature": self.results.get('chamber_temp'),
                "Specific Impulse": self.results.get('isp'),
                "Characteristic Velocity": self.results.get('c_star'),
                "Exit Velocity": self.results.get('exit_velocity'),
                "Mass Flow Rate": self.results.get('mass_flow_rate'),
                "Gamma": self.results.get('gamma'),
                "O/F Ratio": self.results.get('of_ratio'),
            },
            "Geometry": {
                "Throat Diameter": self.results.get('throat_diameter'),
                "Throat Area": self.results.get('throat_area'),
                "Chamber Diameter": self.results.get('chamber_diameter'),
                "Chamber Length": self.results.get('chamber_length'),
                "Chamber Volume": self.results.get('chamber_volume'),
                "Characteristic Length": self.results.get('characteristic_length'),
                "Residence Time": self.results.get('residence_time'),
            }
        }
        return summary


# Initialize sizer with design parameters
sizer = ChamberSizer(
    fuel=p["fuel"],
    oxidiser=p["oxidiser"],
    chamber_pressure=p["chamber_pressure"],
    thrust=p["nominal_thrust"]
)
