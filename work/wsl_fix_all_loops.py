"""Fix both index out-of-bounds issues in pdb_processor.py"""
target = "/home/a/geoevo_work/geoevobuilder/Utils/pdb_processor.py"

with open(target) as f:
    lines = f.readlines()

changes = 0
for i in range(len(lines)):
    line = lines[i]
    # Fix the two for loops that use interval_e without capping
    if 'for i in range(int(interval_b), int(interval_e)):' in line:
        indent = line[:len(line) - len(line.lstrip())]
        lines[i] = f'{indent}for i in range(int(interval_b), min(int(interval_e), len(matrix))):  # len(matrix)\n'
        changes += 1

with open(target, "w") as f:
    f.writelines(lines)
print(f"Fixed {changes} loop(s)")
