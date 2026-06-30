"""Fix indentation in batch_fea function too"""
target = "/home/a/geoevo_work/geoevobuilder/Utils/pdb_processor.py"

with open(target) as f:
    lines = f.readlines()

for i, line in enumerate(lines):
    if i >= 1536 and "sequence_list = sequence_list" in line:
        if not line.startswith("                    "):
            lines[i] = "                    sequence_list = sequence_list[:len(phi_psi_list)]\n"
            print(f"Fixed line {i+1}")
            break

with open(target, "w") as f:
    f.writelines(lines)
print("Done!")
