"""
Consistency check: compare Frobenius traces of the dual pairs that
we have computed to coefficients of eigenforms.
"""
from glob import glob

trace_files = glob('GL2-*-traces.txt')
form_files = glob('GL2-*/forms.txt')

traces = {f.replace('-traces.txt', ''):
          tuple(map(tuple, sage_eval(open(f).read()))) for f in trace_files}

forms = {}
for f in form_files:
    _, l, n, _ = f.split('-', 3)
    P = {ZZ(l)}.union(ZZ(n).prime_factors())
    for s in open(f).readlines():
        T = tuple(sorted((p, a) for p, a in sage_eval(s)[1] if p not in P))
        # if T in forms:
        #     raise ValueError("key {} already found".format(T))
        forms[T] = f.replace('/forms.txt', '')

for f, t in traces.iteritems():
    g = forms[t]
    if f != g: print ('mismatch {} != {}'.format(f, g))
