import glob

cutoffs = {}

files = glob.glob("*_00*.recpot")
files.sort(key=lambda x: len(x))

for f in files:
    el = f.split('_')[0]
    
    if el in cutoffs:
        continue
        
    c, m, fn = None, None, None
    try:
        with open(f, 'r', errors='ignore') as file:
            lines = [file.readline().upper() for _ in range(15)]
            for line in lines:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        val = float(parts[0])
                        if 'COARSE' in line: c = val
                        if 'MEDIUM' in line: m = val
                        if 'FINE'   in line: fn = val
                    except ValueError:
                        pass
    except Exception:
        pass
        
    if c and m and fn:
        cutoffs[el] = {'coarse': c, 'medium': m, 'fine': fn}

print("INTERNAL_CUTOFFS = {")
for el, vals in sorted(cutoffs.items()):
    print(f"    '{el}': {{'coarse': {vals['coarse']}, 'medium': {vals['medium']}, 'fine': {vals['fine']}}},")
print("}")
