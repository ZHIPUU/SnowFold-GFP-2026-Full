#!/usr/bin/env python3
"""Generate upload commands for r17_pipeline.py"""
import base64, json

with open(r'D:\workspace\r17_pipeline.py', 'rb') as f:
    b64 = base64.b64encode(f.read()).decode()

chunks = [b64[i:i+1200] for i in range(0, len(b64), 1200)]
print(f'{len(chunks)} chunks')

bat_lines = []
for i, chunk in enumerate(chunks):
    mode = 'wb' if i == 0 else 'ab'
    # Build remote command with chr() to avoid quote issues
    # Path: /root/autodl-tmp/r17_pipeline.py
    path_chrs = [f'chr({ord(c)})' for c in '/root/autodl-tmp/r17_pipeline.py']
    mode_chrs = [f'chr({ord(c)})' for c in mode]
    
    # Split the base64 data into manageable parts
    b64_parts = [chunk[j:j+200] for j in range(0, len(chunk), 200)]
    b64_expr = '+'.join([f'"{p}"' for p in b64_parts])
    
    cmd = f'python3 -c "import base64; open({\"+\".join(path_chrs)},{\"+\".join(mode_chrs)}).write(base64.b64decode({b64_expr}))"'
    
    full = f'gssh exec 9ca7acb1b94c "{cmd}"'
    
    # Write to bat file
    with open(f'D:/workspace/chunk_{i}.sh', 'w') as f:
        f.write(full)
    
    print(f'  Chunk {i+1}: {len(chunk)} chars -> D:/workspace/chunk_{i}.sh')

print('Done')
