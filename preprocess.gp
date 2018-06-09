\\ Filter out reducible representations

\\ charlist = [[a_1, chi(a_1)], ... [a_m, chi(a_m)]]
\\ identify the Conrey character of chi (for the
\\ given root of unity z in (Z/lZ)^\times)
identify_character(H, z, charlist) =
{
   my(n = H.mod, m = znorder(z), o, candidate = []);
   for(c = 0, n - 1,
      if(gcd(c, n) > 1, next);
      o = charorder(H, c);
      if(m % o != 0, next);
      for(j = 1, #charlist,
	 if(z^(chareval(H, c, charlist[j][1]) * m) != charlist[j][2],
	    next(2)));
      if(candidate == [],
	 candidate = c,
	 error("multiple characters with given values")));
   if(candidate == [],
      error("no character with given values"));
   candidate;
}

\\ Non-proven check for reducibility of the representation attached to
\\ the given form.  TODO: for CM forms, does the fact that many
\\ coefficients vanish make it likely that we will guess incorrectly
\\ that the representation is reducible?
check_reducible(form) =
{
   my(l, H, chi, z, k, aplist, n,
      H1, H2, c2, n1, n2, o1, o2,
      c1prim, c2prim, p, pairs = [], pair);
   [l, H, chi, z, k, aplist] = form;
   n = H.mod;
   for(c1 = 0, n - 1,
      if(gcd(c1, n) > 1, next);
      c2 = chardiv(H, chi, c1);
      n1 = zncharconductor(H, c1);
      n2 = zncharconductor(H, c2);
      if(n % (n1 * n2) != 0, next);
      o1 = charorder(H, c1);
      if((l - 1) % o1 != 0, next);
      o2 = charorder(H, c2);
      if((l - 1) % o2 != 0, next);
      [H1, c1prim] = znchartoprimitive(H, c1);
      [H2, c2prim] = znchartoprimitive(H, c2);
      for(j = 1, #aplist,
	 p = aplist[j][1];
	 if(p == l || n % p == 0, next);
	 if(z^(chareval(H1, c1prim, p) * (l - 1))
	    + p^(k-1) * z^(chareval(H2, c2prim, p) * (l - 1))
	    != aplist[j][2], next(2)));
      if(type(c1prim) != "t_INT", c1prim = znconreyexp(H1, c1prim));
      if(type(c2prim) != "t_INT", c2prim = znconreyexp(H2, c2prim));
      \\ The order of the pairs does not matter if the (k-1)st power
      \\ of the cyclotomic character is trivial.
      pairs = setunion(pairs, [if((k - 1) % (l - 1) == 0,
				  concat(vecsort([[n1, c1prim], [n2, c2prim]])),
				  [n1, c1prim, n2, c2prim])]));
   if(#pairs > 1,
      error("multiple candidates: ", pairs));
   pairs;
}

\\ Proven check for irreducibility.  TODO: we should raise an error if
\\ irreducibility cannot be verified.  This happens for a number of CM
\\ forms because they have many coefficients equal to 0.
prove_irreducible(l, H, chi, z, k, n, aplist) =
{
   my(p, a_p, chi_p);
   chi = zncharinduce(H, chi, n);
   H = znstar(n, 1);
   for(i = 1, #aplist,
      [p, a_p] = aplist[i];
      if(n % p == 0, next);
      chi_p = z^(chareval(H, chi, p) * (l - 1));
      if(#polrootsmod('x^2 - a_p * 'x + chi_p * p^(k - 1), l) == 0,
	 return(1)));
   print("warning: cannot prove irreducibility: [n, k, l] = ", [n, k, l]);
}

gather_reducible(pairs) =
{
   my(map = Map(), forms, form, chars,
      l, H, chi, z, k, aplist, key);
   for(i = 1, #pairs,
      [form, chars] = pairs[i];
      [l, H, chi, z, k, aplist] = form;
      key = concat([[l], chars, [k]]);
      mapput(map, key, if(mapisdefined(map, key, &forms),
			  concat(forms, [[H.mod, aplist]]),
			  [[H.mod, aplist]])));
   Col(Mat(map))~;
}

\\ Helper function for gather_irreducible.
aplist_matches(n, ap, bp) =
{
   my(P = setminus(setintersect(Set([c[1] | c <- ap]),
				Set([c[1] | c <- bp])),
		   factor(n)[,1]~));
   ap = Map(matconcat(ap~));
   bp = Map(matconcat(bp~));
   for(i = 1, #P,
      if(mapget(ap, P[i]) != mapget(bp, P[i]), return(0)));
   return(1);
}

gather_irreducible(forms) =
{
   my(l, H, chi, z, k, aplist, nf, n0, chi0,
      key, found, forms2, match, ng, bplist,
      list = []);
   for(i = 1, #forms,
      [l, H, chi, z, k, aplist] = forms[i];
      nf = H.mod;
      n0 = znconreyconductor(H, chi);
      if(type(n0) == "t_VEC", n0 = n0[1]);
      chi0 = chi % n0;
      key = [l, n0, chi0, k];
      found = 0;
      for(j = 1, #list,
	 if(key != list[j][1], next);
	 forms2 = list[j][2];
	 match = 0;
	 for(m = 1, #forms2,
	    [ng, bplist] = forms2[m];
	    if(aplist_matches(lcm(nf, ng) * l, aplist, bplist),
	       match++));
	 if(match == 0, next,
	    match == #forms2,
	    list[j][2] = concat(forms2, [[nf, aplist]]); found++,
	    error("form matches some but not all other forms")));
      if(found == 0, list = concat(list, [[key, [[nf, aplist]]]]),
	 found > 1, error("found multiple matches")));
   Set(list);
}

/*
  Preprocess the given eigenforms: compute characters and decide
  reducibility and pairwise isomorphism.

  looks_reducible: list of pairs of the shape
    [[l, n_1, chi_1, n_2, chi_2, k], [form_1, form_2, ..., form_j]]
  where each chi_i is primitive of conductor n_i
  and each form_i looks like the representation is
  isomorphic to chi_1 + chi_2 * chi_l^(k-1)

  looks_irreducible: list of pairs of the shape
    [[l, n_0, chi_0, k], [form_1, form_2, ..., form_j]]
  where chi_0 is primitive of conductor n_0 and the forms look like
  their representations are mutually isomorphic

  Forms are given as [n, aplist] and are guaranteed to have weight k,
  level n divisible by n_0 and the specified character (mod n).
*/
preprocess(forms) =
{
   my(form, l, H, chi, z, n, k,
      result, chars, aplist, charlist, params, nf,
      looks_reducible = [], looks_irreducible = []);
   \\ First checks for reducibility.
   for(i = 1, #forms,
      form = forms[i];
      [l, n, k, aplist, charlist] = form;
      H = znstar(n, 1);
      z = znprimroot(l);
      chi = identify_character(H, z, charlist);
      form = [l, H, chi, z, k, aplist];
      result = check_reducible(form);
      if(result != [],
	 chars = result[1];
	 looks_reducible = concat(looks_reducible, [[form, chars]]),
	 looks_irreducible = concat(looks_irreducible, [form])));
   \\ Gather the forms with probably reducible (resp. irreducible)
   \\ representations into groups that look like they should
   \\ have the same representations.
   looks_reducible = gather_reducible(looks_reducible);
   looks_irreducible = gather_irreducible(looks_irreducible);
   \\ For each form that looks like its representation is irreducible,
   \\ verify that the representation is indeed irreducible.
   for(i = 1, #looks_irreducible,
      [params, forms] = looks_irreducible[i];
      [l, n, chi, k] = params;
      H = znstar(n, 1);
      z = znprimroot(l);
      for(j = 1, #forms,
	 [nf, aplist] = forms[j];
	 prove_irreducible(l, H, chi, z, k, nf, aplist)));
   \\ TODO: prove that the representations that look reducible are
   \\ really the sum of the "guessed" characters, and that the
   \\ irreducible ones that look isomorphic are really isomorphic.

   \\ Once we have proved that each form is indeed the sum of two
   \\ characters, we can throw away the q-expansions.
   for(i = 1, #looks_reducible,
      [params, forms] = looks_reducible[i];
      looks_reducible[i] = params);
   \\ Pick the form with smallest level in each list.
   for(i = 1, #looks_irreducible,
      [params, forms] = looks_irreducible[i];
      looks_irreducible[i] = concat(params, vecsort(forms, 1)[1]));
   [looks_reducible, looks_irreducible];
}

\r cusp-forms
\r modular-groups
\r init

labelcount = Map();

/*
  Currently, this skips 4 forms for the following [n, k, l]:
  [23, 4, 2]: weight 1 (multiplicity 2)
  [23, 5, 5]: genus 41
  [24, 3, 5]: genus 73
  [24, 3, 5]: genus 73

  Furthermore, the case [13, 2, 3] yields a curve of genus 0
  (see explanation in level_13_weight_2_mod_3.gp).
*/
output_files(forms) =
{
   my(red, irred, l, n_0, chi_0, k, n, aplist);
   [red, irred] = preprocess(forms);
   \\ TODO: output data for reducible forms
   for(i = 1, #red,
      [l, n_1, chi_1, n_2, chi_2, k] = red[i];
      print(red[i]));
   \\ Output data for irreducible forms
   file = fileopen("irreducible.gp", "w");
   for(i = 1, #irred,
      \\ print(irred[i]);
      [l, n_0, chi_0, k, n, aplist] = irred[i];
      if(l == 2 && k == 4,
	 print([n, k, l], ": skipping because of weight 1");
	 next);
      n_1 = if(k == 2, n, n * l);
      Z_0 = znstar(n_0, 1);
      Z_l = znstar(l, 1);
      Z = znstar(n_1, 1);
      chi = zncharinduce(Z_0, chi_0, Z);
      if(k != 2,
	 \\ [1]~ encodes the cyclotomic character, so
	 \\ [k - 2]~ encodes its (k - 2)-nd power
	 if(l != 2,
	    chi = charmul(Z, chi, zncharinduce(Z_l, [k - 2]~, Z)));
	 aplist = [x | x <- aplist, x[1] != l]);
      z = if(l == 2, 1, znconreyexp(Z_l, [1]~));
      H = charker(Z, chi);
      G = Gamma_H(Z, H);
      g = modular_group_genus(G);
      print([n, k, l], ": genus ", g);
      \\ Provisional label: GL2-l-N-beta-c
      \\ (2 stands for the dimension)
      beta = (k - 2) % (l - 1);
      if(!mapisdefined(labelcount, [l, n, beta], &c),
	 c = 1, c++);
      mapput(labelcount, [l, n, beta], c);
      label = concat(["GL2-", l, "-", n, "-", beta, "-", c]);
      eps = [[x, chareval(Z, chi, x, [z, l - 1])] | x <- Z.gen];
      for(j = 0, #aplist,
	 coefs = aplist[1..j];
	 iferr(f = eigenform(G, 2, l, eps, coefs),
	       e, next, errname(e) == "e_USER");
	 s = concat(["eigenform(Gamma_H(znstar(", n_1, ", 1), ", H, "), ",
		     2, ", ", l, ", ", eps]);
	 if(coefs != [], s = concat([s, ", ", coefs]));
	 s = concat(s, ")");
	 filewrite(file, s);
	 extern(concat("mkdir -p ", label));
	 makefile = fileopen(concat(label, "/Makefile"), "w");
	 filewrite(makefile, concat("# genus ", g));
	 filewrite(makefile, concat("FORM = ", s));
	 filewrite(makefile, "PRIMES = 20");
	 filewrite(makefile, concat("MAX_DEGREE = ", l + 1));
	 filewrite(makefile, "include /home/peter/share/modgalrep/Makefile.form");
	 fileclose(makefile);
	 break));
   fileclose(file);
}
