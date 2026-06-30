import os
for root, dirs, files in os.walk(r'D:\生信\2026Protein Design\work'):
    # 排除 .git
    dirs[:] = [d for d in dirs if d != '.git']
    for f in files:
        path = os.path.join(root, f)
        size = os.path.getsize(path)
        print(f"{size:>12}  {path}")