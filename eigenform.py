"""
Eigenforms given by `q`-expansions
"""
from sage.misc.cachefunc import cached_method
from sage.modular.arithgroup.all import Gamma1
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.structure.element import Element

from copy import copy

class Eigenform(Element):
    def __init__(self, collection, apdict):
        Element.__init__(self, collection)
        self._apdict = copy(apdict)
        self._minimal_form = None

    def group(self):
        return self.parent().group()

    def level(self):
        return self.parent().level()

    def weight(self):
        return self.parent().weight()

    def base_ring(self):
        return self.parent().base_ring()

    def character(self):
        return self.parent().character()

    def primitive_character(self):
        return self.parent().primitive_character()
    
    @cached_method
    def coefficient(self, m):
        if m <= 0:
            # TODO: we are assuming the form is a cusp form
            return 0
        elif m == 1:
            return 1
        from sage.arith.misc import factor
        fact = factor(m)
        if len(fact) == 1:
            p, j = fact[0]
            a_p = self._apdict[p]
            # prime power
            if j == 1:
                return a_p
            else:
                return (a_p * self.coefficient(p**(j - 1))
                        - self.character()(p) * p**(self.weight() - 1) * self.coefficient(p**(j - 2)))
        else:
            from sage.misc.all import prod
            return prod(self.coefficient(p**k) for p, k in fact)

    def q_expansion(self, m):
        # This is somewhat inefficient, but we shouldn't care.
        return [self.coefficient(j) for j in range(m)]

    def __hash__(self):
        return hash(tuple(sorted(self._apdict.iteritems())))

    def _repr_(self):
        return repr((self.level(), self.weight(), sorted(self._apdict.iteritems())))

    def _cmp_(self, other):
        B = self.parent().sturm_bound() + 1
        return cmp(self.q_expansion(B), other.q_expansion(B))

    def _raise(self, p, level, weight):
        """
        Helper function for :meth:`raise_level` and :meth:`raise_weight`.
        """
        from eigenform_collection import EigenformCollection
        F = self.base_ring()
        R = PolynomialRing(F, 't')
        t = R.gen()
        chi = self.character()
        poly = t**2 - self.coefficient(p)*t + chi(p) * p**(self.weight() - 1)
        roots = poly.roots(multiplicities=False)
        apdict = copy(self._apdict)
        C = EigenformCollection(Gamma1(level), weight, F, chi.extend(level))
        for a in roots:
            apdict[p] = a
            yield C(apdict)

    def raise_level(self, p):
        """
        Return all eigenforms that are linear combinations of the images
        of ``self`` under the degeneracy maps `i_1` and `i_p`.
        """
        return self._raise(p, self.level() * p, self.weight())

    def raise_weight(self):
        """
        Return all eigenforms that are linear combinations of the images
        of ``self`` under the Frobenius map and multiplication by the
        Hasse invariant.
        """
        if self.weight() != 1:
            raise ValueError("weight must be 1")
        l = self.base_ring().characteristic()
        return self._raise(l, self.level(), l)

    def minimal_form(self):
        return self._minimal_form

    def set_minimal_form(self, f):
        self._minimal_form = f
