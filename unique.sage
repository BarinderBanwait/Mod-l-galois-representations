"""
Build a database of sets of eigenforms modulo `l` having
isomorphic Galois representations.
"""
from eigenform_collection import EigenformCollection

def convert_form(f, l):
    (n, k, i, d, aplist, charlist) = f
    F = FiniteField(l)
    # We forget about i and d, and we convert aplist and charlist to
    # dictionaries.
    apdict = {p: F(a_p) for p, a_p in aplist}
    chardict = {d: F(chi_d) for d, chi_d in charlist}
    return EigenformCollection(Gamma1(n), k, F, chardict)(apdict)

def import_forms(filename, l):
    forms = [sage_eval(s) for s in file(filename).readlines()]
    # forms is a list of eigenforms mod l in the following format:
    # [n, k, i, d, aplist, charlist]
    return [convert_form(f, l) for f in forms]

def unique_with_char(l, chi, forms, max_level=24):
    print(chi)
    print('>> number of forms: %s' % len(forms))
    levels_weights = {(f.level(), f.weight()) for f in forms}
    forms_dict = {(n, k): [f for f in forms if f.level() == n and f.weight() == k]
                  for n, k in levels_weights}
    for n, k in sorted(levels_weights):
        for f in forms_dict[n, k]:
            print("form %s" % f)
            if f.minimal_form() is None:
                f.set_minimal_form(f)
            # Raise levels
            for p in primes(max_level // n + 1):
                if p == l:
                    continue
                # Raise the level at p.
                forms_np = forms_dict.get((n * p, k), [])
                for g in f.raise_level(p):
                    forms_g = [h for h in forms_np if h == g]
                    if forms_g == []:
                        print("warning: form %s not found" % g)
                        forms_np.append(g)
                        # TODO: check whether the q-expansion is there
                        # in a different weight
                    else:
                        # Modify the form that was already in the
                        # list, not the one we got by level raising.
                        g = forms_g[0]
                        assert g.minimal_form() in [None, f.minimal_form()], \
                            "inconsistent minimal form"
                    g.set_minimal_form(f.minimal_form())
            # TODO: find weight 1 forms.  Does this work for non-cusp forms?
            # Raise weights
            if k == 2:
                # TODO: interaction with k = l + 1
                pass
            # TODO: interaction between k = 3 and k = 4 for l = 2
    print('>> number of unique representations: %s'
          % len({f.minimal_form() for f in forms}))

def unique(l, max_level=24):
    F = FiniteField(l)
    filename = 'mod-{}-forms.txt'.format(l)
    forms = import_forms(filename, l)
    chars = {f.primitive_character() for f in forms}
    print('number of forms: %s' % len(forms))
    print('number of characters: %s' % len(chars))
    count = 0
    for chi in chars:
        # Filter out the forms with the correct character.
        forms_chi = [f for f in forms if f.primitive_character() == chi]
        unique_with_char(l, chi, forms_chi, max_level)
        count += len({f.minimal_form() for f in forms_chi})
    print('number of unique representations: %s' % count)

# forms at lines 8 and 65 

forms = import_forms("mod-5-forms.txt", 5)
f7 = forms[7]
f21 = forms[64]
