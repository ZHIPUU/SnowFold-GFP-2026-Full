"""Fix: cap interval_e to matrix length in loadpdbs"""
target = "/home/a/geoevo_work/geoevobuilder/Utils/pdb_processor.py"

with open(target) as f:
    code = f.read()

# Fix in loadpdbs: cap interval_e
old1 = """                for i in range(int(interval_b), int(interval_e)):  # len(matrix)
                    # surrounding residues"""
new1 = """                max_idx = len(matrix)
                for i in range(int(interval_b), min(int(interval_e), max_idx)):  # len(matrix)
                    # surrounding residues"""
code = code.replace(old1, new1)

# Fix in triangularexpression: same pattern
old2 = """                for i in range(int(interval_b), int(interval_e)):  # len(matrix
                    # surrounding residues
                    NA = ns.search(chain[int(matrix.iloc[i, 3])]["CA"].coord, 12)  # Neighbor Atoms"""
new2 = """                max_idx = len(matrix)
                for i in range(int(interval_b), min(int(interval_e), max_idx)):  # len(matrix
                    # surrounding residues
                    NA = ns.search(chain[int(matrix.iloc[i, 3])]["CA"].coord, 12)  # Neighbor Atoms"""
code = code.replace(old2, new2)

with open(target, "w") as f:
    f.write(code)

print("Interval capping applied!")
