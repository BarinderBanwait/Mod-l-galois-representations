# Mod-l-galois-representations

In this repository we store raw data for the collection [reps] in the database [mod_l_galois] for the [LMFDB](https://github.com/LMFDB/lmfdb).

Each line in a data file will contain a single comma-separated list (enclosed in square brackets) correspond to one mod l Galois representation. The elements of the list are as follows, in this exact order:

* **numberfield_label** (string): LMFDB label of the field whose absolute Galois group is being represented;
* **dimension** (integer): the dimension of the representation;
* **alg_group** (string): the algebraic group for the codomain (e.g. "GL", "GSp");
* **[finitefield_char, finitefield_degree, finitefield_poly]** (int, int, list of integers): coefficient field description. In general here there will be the descriptor for a finite ring;
* **conductor** (integer): prime-to-ell Artin conductor of the representation;
* **primes_dividing_conductor_with_exponents** (list of pair of ints): [[p,e]: p prime dividing conductor, e exponent];
* **power_cyclo** (integer): exponent of the cyclotomic character of the representation, integer modulo (finitefield_char - 1), in the determinant of the representations;
* **determinant** (string): label of the mod ell Galois character which is the determinant of the representation;
* **image_type** (string): human sensible descriptor of the image, e.g "big" = contains SL_n; "Borel", "cyclic" etc;
* **image_label** (string): label of group. For character "Cn" where n is the order, for dimension 2 over prime field use Sutherland's label, "C1" for trivial, in other cases we need a labeling scheme, so "to be decided";
* **image_attribute** (list of Booleans): 

      [ is_solvable, is_transitive, is_irred, is_absolutely_irred, is_surjective ], 

 0 = False , 1 = True, transitive is referred to the action on the primitive vectors in the representations space;
* **image_order** (integer): size of the image;
* **degree_proj_field** (integer): degree of minimal coefficient field of the projective representation;
* **projective_type** (string): human sensible descriptor of the projective image, e.g. "A5" or "trivial" similar;
* **projective_label** (string): label of group. For character "Cn" where n is the order, for dimension 2 over prime field use Sutherland's label, "C1" for trivial, in other cases we need a labeling scheme, so "to be decided";
* **projective_image_attribute** (list of Booleans): 
      
      [ is_transitive, is_surjective ], 
      
 0 = False , 1 = True, transitive is referred to the action on the set of lines (submodule generated by a primitive vector) in the representations space;
* **good_prime_list** (list): list of lists of the form 

      [p, trace, determinant, charpoly, factorization_charpoly, order, projective_order, representative_matrix]
   
   where 
   - p is int(p) or string if p is a prime in a number field different from Q;
   - trace, determinant are lists of integers which represent the corresponding elements in the coefficient field;
   - charpoly (list of lists): coefficients of the characteristic polynomial of the Frobenius matrix at p represented as above;
   - factorization_charpoly (list of pairs): every pair is a lists consisting of a list of coefficients of an irreducible factor of charpoly and its exponent in the factorization;
   - order, projective order (integers): order, projective order of the Frobenius matrix at p;
   - representative_matrix (list of lists): list of length dimension^2 consisting of coefficients of a representative for the conjugacy class (inside the image) of a Frobenius matrix at p;

for primes p< 100 (we will decide later about number fields) **unramified** (including finitefield_char is the representation is unramified at finitefield_char);
* **real_infinite_list** (list): same as for good_prime_list where the firts entry is a string representing a real root of the defining polynomial of numberfield_label;
* **bad_prime_list** (list): list of lists. For each ramified prime give a list of 4 elements:

      [p, inertia_invariants, inertia_coinvariants, type] 
 
 where p is a prime, inertia_invariants and inertia_coinvariants are lists of the form
 
          [trace, determinant, charpoly, factorization_charpoly, order, projective_order, representative_matrix]

 everything represented with the same format as good_prime_list, and type is a list of descriptors. For 2-dimensional representation type = [is_reducible, is_decomposable] (list of Booleans).
         
* **poly_ker** (list of lists of strings): if this list has length one, it contains the list of coefficients of an irreducible polynomial whose splitting field is isomorphic to the splitting field of the representation. If the length is greater than 1, it contains lists of coefficients of irreducible polynomials whose common splitting field is isomorphic to the splitting field of the representation;
* **poly_proj_ker** (list of lists of strings): same, but for the projective representation;
* **related_objects** (list of pairs of strings): each pair is of the form ["object_type", "object_label"] for objects related to representation, e.g. ["Dirichlet character", "label_of"]. This is a list containing a finite number of samples.

We insist that all strings be enclosed in straight double quotes, "like so"; unknown strings should look like "".
