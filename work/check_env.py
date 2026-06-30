import sys
import torch
print("torch:", torch.__version__)
print("cuda:", torch.cuda.is_available())
if torch.cuda.is_available():
    print("device:", torch.cuda.get_device_name(0))
    a = torch.zeros(2, 2, device='cuda')
    print("GPU tensor OK")
try:
    import esm
    print("esm: OK")
except Exception as e:
    print("esm error:", e)
print("DONE")
sys.stdout.flush()