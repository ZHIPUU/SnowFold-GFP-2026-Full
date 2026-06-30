"""Fix GeoEvoBuilder pdb_processor.py matrix/phi_psi length mismatch"""
target = "/home/a/geoevo_work/geoevobuilder/Utils/pdb_processor.py"

with open(target) as f:
    code = f.read()

# Fix: ensure sequence_list matches phi_psi_list length
old = "    matrix = ["
new = """    if len(phi_psi_list) < len(sequence_list):
        sequence_list = sequence_list[:len(phi_psi_list)]
    matrix = ["""
code = code.replace(old, new)

with open(target, "w") as f:
    f.write(code)
print("Matrix fix applied")
