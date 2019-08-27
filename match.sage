"""
Consistency check: compare Frobenius traces and characters of the
dual pairs that we have computed to coefficients of eigenforms.
"""
from glob import glob

labels = [f.replace('.gp', '') for f in glob('GL2-*.gp')]

# Compare Frobenius traces.

forms = {}
for f in glob('GL2-*/forms.txt'):
    g = f.replace('/forms.txt', '')
    _, l, n, _ = f.split('-', 3)
    P = {ZZ(l)}.union(ZZ(n).prime_factors())
    for s in open(f).readlines():
        t = tuple(sorted((p, a) for p, a in sage_eval(s)[1] if p not in P))
        if t in forms and forms[t] != g:
            raise ValueError("key {} already found".format(t))
        forms[t] = g

for f in labels:
    t = tuple(map(tuple, sage_eval(open(f + '-traces.txt').read())))
    g = forms[t]
    if f != g: print('mismatch {} != {}'.format(f, g))

# Compare characters.

for f in labels:
    c1 = open(f + '-character.txt').read().strip()
    c2 = open(f + '/character.txt').read().strip()
    if c1 != c2:
        print('mismatch: {} != {}'.format(c1, c2))
