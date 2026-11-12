#!/usr/bin/env python3
"""sass_split.py <obj-or-cubin> <outdir>: one normalised SASS file per function, plus index.tsv.

Normalisation drops the per-instruction address comments and the raw encoding lines, so two
functions compare equal iff their instruction streams are identical -- across two builds, or
between two kernels of one build whose names differ. index.tsv: name, kind, instruction count,
sha1 of the normalised body, file; kind = SUB if the demangled name mentions SubCoordinates.

Every __global__ carries its own copy of the device functions it calls, so the SASS of a flow TU
is a census of what ptxas really compiled (see docs/NUMTRACER_DEVICE_COMPILE.md). Comparing the
hashes of two builds proves a change kernel-neutral without running it.
"""
import hashlib, os, re, subprocess, sys

obj, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
sass = subprocess.run(['cuobjdump', '-sass', obj], capture_output=True, text=True, check=True).stdout
funcs, cur, body = [], None, []
for line in sass.splitlines():
    m = re.match(r'\s*Function : (\S+)', line)
    if m:
        if cur: funcs.append((cur, body))
        cur, body = m.group(1), []
        continue
    if cur is None: continue
    m = re.match(r'\s*/\*[0-9a-f]+\*/\s+(.*?)\s*/\* 0x[0-9a-f]+ \*/\s*$', line)
    if m:
        body.append(m.group(1))
if cur: funcs.append((cur, body))
names = [f[0] for f in funcs]
dem = subprocess.run(['c++filt'], input='\n'.join(names), capture_output=True, text=True).stdout.splitlines()
with open(os.path.join(out, 'index.tsv'), 'w') as idx:
    for i, ((n, b), d) in enumerate(zip(funcs, dem)):
        txt = '\n'.join(b) + '\n'
        h = hashlib.sha1(txt.encode()).hexdigest()[:16]
        fn = f'{i:03d}.sass'
        open(os.path.join(out, fn), 'w').write(txt)
        kind = 'SUB' if 'SubCoordinates' in d else '-'
        idx.write(f'{n}\t{kind}\t{len(b)}\t{h}\t{fn}\n')
print(f'{len(funcs)} functions, {sum(len(b) for _, b in funcs)} instructions')
