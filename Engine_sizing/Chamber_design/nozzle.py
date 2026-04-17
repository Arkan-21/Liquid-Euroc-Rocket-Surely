import ezdxf
import numpy as np
import Chamber_sizer

doc = ezdxf.readfile(r"C:\Users\arkan\OneDrive\Desktop\Liquid Euroc Rocket Surely\Liquid-Euroc-Rocket-Surely\Engine_sizing\Chamber_design\Nozzle.dxf")
msp = doc.modelspace()

coords = []

for entity in msp:
    if entity.dxftype() == "LWPOLYLINE":
        points = entity.get_points()  # (x, y, start_width, end_width, bulge)
        for p in points:
            coords.append([p[0], p[1]])  # 2D profile (x,R)

coords_array = np.array(coords)

print(coords_array)

# Extract radius
r = coords_array[:, 1]

# Find throat (minimum radius)
r_throat = Chamber_sizer.engine.Rt

# Compute area ratio
Area_ratios = (r / r_throat) ** 2
print(Area_ratios)