"""
Parents for eigenforms given by `q`-expansions.
"""
from sage.structure.factory import UniqueFactory
from sage.structure.parent import Parent

from eigenform import Eigenform

class EigenformCollection_class(Parent):
    """
    Parent class for eigenforms with given group, weight, base ring
    and character.
    """

    Element = Eigenform

    def __init__(self, group, k, base_ring, char):
        from sage.categories.finite_sets import FiniteSets
        self._group = group
        self._weight = k
        self._base_ring = base_ring
        if isinstance(char, dict):
            G = DirichletGroup(group.level(), base_ring)
            L = [chi for chi in G if all(chi(x) == y for x, y in char.iteritems())]
            if len(L) > 1:
                raise ValueError('character not uniquely determined')
            chardict = L[0]
        self._character = char
        Parent.__init__(self, category=FiniteSets())

    def group(self):
        return self._group
        
    def level(self):
        return self.group().level()

    def weight(self):
        return self._weight

    def base_ring(self):
        return self._base_ring

    def character(self):
        return self._character

    def primitive_character(self):
        return self.character().primitive_character()

    def sturm_bound(self):
        from sage.modular.arithgroup.all import Gamma0
        # TODO: justify (see article of Buzzard and Stein)
        return Gamma0(self.level()).sturm_bound(self.weight())


class EigenformCollectionFactory(UniqueFactory):
    def create_key(self, group, k, base_ring, char):
        from sage.rings.integer_ring import ZZ
        k = ZZ(k)
        if isinstance(char, dict):
            from sage.modular.dirichlet import DirichletGroup
            G = DirichletGroup(group.level(), base_ring)
            L = [chi for chi in G if all(chi(x) == y for x, y in char.iteritems())]
            if len(L) > 1:
                raise ValueError('character not uniquely determined')
            char = L[0]
        return (group, k, base_ring, char)

    def create_object(self, version, key):
        return EigenformCollection_class(*key)
    
EigenformCollection = EigenformCollectionFactory("EigenformCollection")
