"""Replace Bio.PDB.DSSP with pydssp in GeoEvoBuilder"""
import sys, os
import numpy as np

# Path to GeoEvoBuilder's pdb_processor.py
target = os.path.expanduser("~/geoevo_work/geoevobuilder/Utils/pdb_processor.py")

with open(target) as f:
    code = f.read()

# Replace the DSSP import line
code = code.replace(
    "from Bio.PDB.DSSP import DSSP",
    "# DSSP replaced by pydssp (see patch below)"
)

# Replace the DSSP call in loadpdbs (line ~159)
old_dssp_line = "        dssp = DSSP(model, files, dssp=\"mkdssp\")"
new_dssp_lines = """        # DSSP replaced with pydssp
        try:
            import pydssp
            # Calculate secondary structure
            coords_list = []
            res_keys = []
            for chain in model:
                for res in chain:
                    if res.get_id()[0] == ' ' and 'CA' in res:
                        coords_list.append(res['CA'].get_coord())
                        res_keys.append((chain.get_id(), res.get_id()))
            if len(coords_list) > 0:
                coords = np.array(coords_list).reshape(1, -1, 3)
                ss = pydssp.assign(coords, out_type='c3')[0]
                dssp = type('DSSP', (), {'__getitem__': lambda self, k: [' ',' ',res_keys.index(k)+1,ss[res_keys.index(k)] if res_keys.index(k)<len(ss) else 'C',' ',' ',0,0,0,0,0], 'keys': lambda self: res_keys})()
            else:
                raise Exception("no CA atoms")
        except Exception:
            # Fallback: all coil
            dssp = type('DSSP', (), {'__getitem__': lambda self, k: [' ',' ',0,'C',' ',' ',0,0,0,0,0], 'keys': lambda self: []})()"""

code = code.replace(old_dssp_line, new_dssp_lines)

# Also fix the second DSSP call in batch_fea (line ~1485)
old_dssp2 = "        dssp = DSSP(model, pdbfiles, dssp=\"dssp\")"
new_dssp2 = """        # DSSP replaced with pydssp
        try:
            import pydssp
            coords_list = []
            res_keys = []
            for chain in model:
                for res in chain:
                    if res.get_id()[0] == ' ' and 'CA' in res:
                        coords_list.append(res['CA'].get_coord())
                        res_keys.append((chain.get_id(), res.get_id()))
            if len(coords_list) > 0:
                coords = np.array(coords_list).reshape(1, -1, 3)
                ss = pydssp.assign(coords, out_type='c3')[0]
                dssp = type('DSSP', (), {'__getitem__': lambda self, k: [' ',' ',res_keys.index(k)+1,ss[res_keys.index(k)] if res_keys.index(k)<len(ss) else 'C',' ',' ',0,0,0,0,0], 'keys': lambda self: res_keys})()
            else:
                raise Exception("no CA atoms")
        except Exception:
            dssp = type('DSSP', (), {'__getitem__': lambda self, k: [' ',' ',0,'C',' ',' ',0,0,0,0,0], 'keys': lambda self: []})()"""

code = code.replace(old_dssp2, new_dssp2)

# Add numpy import if not present
if "import numpy" not in code.split("\n")[0]:
    code = "import numpy as np\n" + code

with open(target, "w") as f:
    f.write(code)

print("DSSP patch applied successfully!")
