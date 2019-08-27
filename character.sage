from dual_pairs.dual_pair_import import dual_pair_import

def dirichlet_character_label(chi):
    r"""
    Return the LMFDB label of :math:`\chi`.

    INPUT:

    - ``chi`` -- a :class:`DirichletCharacter` with values in a finite
      field

    EXAMPLES::

        sage: G = DirichletGroup(39, GF(3))
        sage: chi = G([2, 2])
        sage: dirichlet_character_label(chi)
        '3-39.38'

        sage: G = DirichletGroup(55, GF(11))
        sage: chi = G([1, 2])
        sage: dirichlet_character_label(chi)
        '11-55.46'
    """
    G = chi.parent()
    m = G.modulus()
    l = G.base_ring().characteristic()
    z = G.zeta()
    o = G.zeta_order()
    H = pari('idealstar(,{},2)'.format(m))
    cyc = H.getattr('cyc')
    gen = H.getattr('gen')
    v = [chi(g).log(z) * c/o for c, g in zip(cyc, gen)]
    c = pari('znstar({},1)'.format(m)).znconreyexp(v)
    return "{}-{}.{}".format(l, m, c)

if len(sys.argv) != 2:
    raise RuntimeError('usage: character.sage FILE')

D = dual_pair_import(sys.argv[1])
chi = D.dirichlet_character()
print(dirichlet_character_label(chi))
