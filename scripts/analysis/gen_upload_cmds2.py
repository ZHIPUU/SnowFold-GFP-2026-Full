#!/usr/bin/env python3
import base64

with open(r'D:\workspace\r17_pipeline.py', 'rb') as f:
    b64 = base64.b64encode(f.read()).decode()

chunks = [b64[i:i+1200] for i in range(0, len(b64), 1200)]
print(f'{len(chunks)} chunks')

for i, chunk in enumerate(chunks):
    mode = 'wb' if i == 0 else 'ab'
    path = '/root/autodl-tmp/r17_pipeline.py'
    
    # Build remote command using string concatenation
    path_chrs = '+'.join([f'chr({ord(c)})' for c in path])
    mode_chrs = '+'.join([f'chr({ord(c)})' for c in mode])
    
    b64_parts = [chunk[j:j+200] for j in range(0, len(chunk), 200)]
    b64_expr = '+'.join([f'"{p}"' for p in b64_parts])
    
    # Simple string building
    inner = 'import base64; open(' + path_chrs + ',' + mode_chrs + ').write(base64.b64decode(' + b64_expr + '))'
    cmd = 'gssh exec 9ca7acb1b94c "python3 -c \\"' + inner + '\\""'
    
    with open(f'D:/workspace/chunk_{i}.sh', 'w') as f:
        f.write(cmd)
    print(f'  Chunk {i+1}: {len(chunk)} chars')

print('Done')
