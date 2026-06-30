#!/home/a/miniconda3/envs/geoevo/bin/python
"""Fix PDB format - add element column for DSSP compatibility"""
import sys, os
sys.path.insert(0, os.path.expanduser("~/geoevo_work"))

with open("1GFL_fixed_Val.pdb") as f:
    lines = f.readlines()

fixed = []
for line in lines:
    if line.startswith("ATOM") or line.startswith("HETATM"):
        if len(line.rstrip("\n")) < 78:
            atom_name = line[12:16].strip()
            element = atom_name[0]
            line = line.rstrip("\n")
            line = f"{line:<78}{element:>2}\n"
    fixed.append(line)

with open("1GFL_fixed_Val.pdb", "w") as f:
    f.writelines(fixed)

n_atom = sum(1 for l in fixed if l.startswith("ATOM"))
print(f"Fixed {n_atom} ATOM lines - element columns added")
