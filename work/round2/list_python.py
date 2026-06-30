import psutil
import os
for p in psutil.process_iter(['pid', 'cmdline', 'create_time']):
    cl = p.info['cmdline']
    if cl and 'python' in (cl[0] or '').lower():
        cmd_short = ' '.join(cl[:3]) if len(cl) >= 3 else ' '.join(cl)
        print(f"PID={p.info['pid']}, cmd={cmd_short}")