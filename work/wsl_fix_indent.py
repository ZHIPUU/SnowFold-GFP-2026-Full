"""Fix GeoEvoBuilder pdb_processor.py - fix indentation of matrix fix"""
target = "/home/a/geoevo_work/geoevobuilder/Utils/pdb_processor.py"

with open(target) as f:
    lines = f.readlines()

# Find and fix the indentation issue
fixed = []
in_matrix_block = False
for i, line in enumerate(lines):
    if "if len(phi_psi_list) < len(sequence_list):" in line:
        lines[i] = "                if len(phi_psi_list) < len(sequence_list):\n"
        lines[i+1] = "                    sequence_list = sequence_list[:len(phi_psi_list)]\n"
        break

with open(target, "w") as f:
    f.writelines(lines)

# Verify
with open(target) as f:
    for i, line in enumerate(f.readlines()):
        if "phi_psi_list" in line:
            print(f"Line {i+1}: {line.rstrip()}")
print("Indentation fixed!")
