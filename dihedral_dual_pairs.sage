from dual_pair_from_dihedral_field import dual_pair_from_dihedral_field

pari.read('dual_pair.gp')
dual_pair_init = pari('dual_pair_init')
dual_pair_export = pari('dual_pair_export')

# TODO: make the self-duality visible
def reduce_dual_pair(D):
    polys1 = D.algebra1()._polys
    polys2 = D.algebra2()._polys
    phi = D.phi()
    Dp = dual_pair_init([polys1, polys2, phi], [], 1)
    return dual_pair_export(Dp)

def print_traces(D, bound=100):
    for p in primes(bound):
        try:
            print((p, D.frobenius_matrix(p).trace()))
        except:
            pass

R.<x> = QQ[]
F = GF(3)

print('GL2-3-7-0-1')
D = dual_pair_from_dihedral_field(x^4 - x^3 + 2*x + 1, F)
for y in reduce_dual_pair(D): print(y)

print('GL2-3-13-1-2')
D = dual_pair_from_dihedral_field(x^4 + 2*x^2 - 3*x + 1, F)
for y in reduce_dual_pair(D): print(y)

print('GL2-3-16-0-1')
D = dual_pair_from_dihedral_field(x^4 - 3*x^2 + 3, F)
for y in reduce_dual_pair(D): print(y)

print('GL2-3-19-0-1')
D = dual_pair_from_dihedral_field(x^4 - x^3 - 3*x^2 + 2*x + 4, F)
for y in reduce_dual_pair(D): print(y)

"""
candidates = [x^4 - 3*x^2 + 3, x^4 - 2*x^3 - 2*x + 1, x^4 - 6*x^2 + 12, x^4 - 2*x^2 + 3, x^4 + 2*x^2 + 3, x^4 - 2*x^2 - 2, x^4 + 2*x^2 - 2, x^4 - 6*x^2 + 18, x^4 - 4*x^2 + 6, x^4 + 4*x^2 + 6, x^4 - 3, x^4 - 6*x^2 - 9, x^4 + 6*x^2 + 6, x^4 - 6*x^2 + 6, x^4 - 4*x^2 - 2, x^4 + 4*x^2 - 2, x^4 - 18, x^4 + 18, x^4 + 6*x^2 + 3, x^4 - 6*x^2 + 3, x^4 - 6, x^4 - 24, x^4 + 6, x^4 + 24]

for f in candidates:
    print(f)
    D = dual_pair_from_dihedral_field(f, F)
    print_traces(D)
"""
