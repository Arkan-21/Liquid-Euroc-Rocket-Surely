"""
Example: Chamber Sizing Tool Usage

This script demonstrates how to use the ChamberSizer class to calculate
engine performance and chamber geometry.
"""

import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Chamber import ChamberSizer
from Design_overview import input_parameters as p
import pprint

# Create a chamber sizer with design parameters
sizer = ChamberSizer(
    fuel=p["fuel"],
    oxidiser=p["oxidiser"],
    chamber_pressure=p["chamber_pressure"],
    thrust=p["nominal_thrust"]
)

# Calculate performance characteristics (O/F ratio can be optimized for Isp)
print("=" * 70)
print("CHAMBER SIZING TOOL - EXAMPLE")
print("=" * 70)
print(f"\nInput Parameters:")
print(f"  Fuel: {p['fuel']}")
print(f"  Oxidiser: {p['oxidiser']}")
print(f"  Chamber Pressure: {p['chamber_pressure']:.2f}")
print(f"  Nominal Thrust: {p['nominal_thrust']:.2f}")

# Scan different O/F ratios to find optimal performance
print("\n" + "=" * 70)
print("SCANNING O/F RATIOS FOR OPTIMAL PERFORMANCE")
print("=" * 70)

of_ratios = [1.5, 2.0, 2.5, 3.0, 3.5]
best_isp = 0
best_of = 2.0

print(f"\n{'O/F Ratio':<12} | {'Tc (K)':<10} | {'Isp (s)':<10} | {'c* (m/s)':<12}")
print("-" * 50)

for of in of_ratios:
    sizer.calculate_performance(of_ratio=of)
    tc = sizer.results['chamber_temp'].value
    isp = sizer.results['isp'].value
    c_star = sizer.results['c_star'].value
    
    print(f"{of:<12.1f} | {tc:<10.0f} | {isp:<10.1f} | {c_star:<12.1f}")
    
    if isp > best_isp:
        best_isp = isp
        best_of = of

# Calculate final design with best O/F ratio
print(f"\n{'→ Optimal O/F Ratio:':<30} {best_of:.1f}")

sizer.calculate_performance(of_ratio=best_of)
sizer.calculate_geometry()

# Display results
print("\n" + "=" * 70)
print("PERFORMANCE CHARACTERISTICS")
print("=" * 70)

results = sizer.results
print(f"\nCombustion Properties:")
print(f"  Chamber Temperature:      {results['chamber_temp']:.2f}")
print(f"  Gamma (γ):                {results['gamma']:.4f}")
print(f"  Molecular Weight:         {results['molecular_weight']:.2f}")

print(f"\nPerformance Parameters:")
print(f"  Characteristic Velocity:  {results['c_star']:.2f}")
print(f"  Specific Impulse:         {results['isp']:.2f}")
print(f"  Exit Velocity:            {results['exit_velocity']:.2f}")
print(f"  Mass Flow Rate:           {results['mass_flow_rate']:.2f}")

print("\n" + "=" * 70)
print("CHAMBER GEOMETRY")
print("=" * 70)

print(f"\nThroat Characteristics:")
print(f"  Throat Diameter:          {results['throat_diameter']:.2f}")
print(f"  Throat Area:              {results['throat_area']:.2f}")

print(f"\nChamber Characteristics:")
print(f"  Chamber Diameter:         {results['chamber_diameter']:.2f}")
print(f"  Chamber Length:           {results['chamber_length']:.2f}")
print(f"  Chamber Volume:           {results['chamber_volume']:.2f}")

print(f"\nFlow Characteristics:")
print(f"  Characteristic Length:    {results['characteristic_length']:.2f}")
print(f"  Residence Time:           {results['residence_time']:.2f}")

print("\n" + "=" * 70)
print("SUMMARY TABLE")
print("=" * 70)

summary = sizer.get_summary()
for category, props in summary.items():
    print(f"\n{category}:")
    for key, value in props.items():
        if value is not None:
            print(f"  {key:<28} {value}")
