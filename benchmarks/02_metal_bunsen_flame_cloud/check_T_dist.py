#!/usr/bin/env python3
import glob
import re

total_cells = 0
hot_cells = 0
min_T = 1e9
max_T = -1e9

for f_path in glob.glob("/home/admin123/OpenFOAM/admin123-13/run/MetalCombustion_Case21_BunsenFlame/processor*/0.3/T"):
    with open(f_path, 'r') as f:
        content = f.read()
    
    # Match internalField List
    match = re.search(r'internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\(([\s\S]*?)\);', content)
    if match:
        n_cells = int(match.group(1))
        vals_str = match.group(2).split()
        vals = []
        for v in vals_str:
            try:
                vals.append(float(v))
            except ValueError:
                pass
        total_cells += len(vals)
        for v in vals:
            min_T = min(min_T, v)
            max_T = max(max_T, v)
            if v > 2000.0:
                hot_cells += 1

print(f"Total Mesh Cells: {total_cells}")
print(f"Hot Cells (T > 2000 K): {hot_cells}")
print(f"Min Temperature: {min_T:.1f} K")
print(f"Max Temperature: {max_T:.1f} K")
