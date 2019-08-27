# Generate LMFDB data for Galois representations attached to modular
# forms over finite fields

# Input for these scripts: the files mod-L-reps.txt for L = 2, 3, 5
# containing eigenforms computed by Andrew Sutherland (not part of
# this repository).

GP ?= gp -q
SAGE ?= sage

all: labels.txt reducible.gp irreducible.gp

.DELETE_ON_ERROR:

reducible.gp: forms.gp
	echo "output_reducible(preprocess(readvec(\"forms.gp\"))[1]);" | $(GP) -q preprocess.gp

labels.txt irreducible.gp: forms.gp
	echo "output_irreducible(preprocess(readvec(\"forms.gp\"))[2]);" | $(GP) -q preprocess.gp

forms.gp: mod-2-forms.txt mod-3-forms.txt mod-5-forms.txt
	> $@
	for l in 2 3 5; do echo "convert(readvec(\"mod-$$l-forms.txt\"),$$l);" | $(GP) convert.gp >> $@; done

mod-%-forms.txt: mod-%-reps.txt
	grep '^<' $< | sed -e 's/</[/g' -e 's/>/]/g' -e 's/,$$//' > $@

GL2-%-traces.txt: GL2-%.gp
	$(SAGE) -c "from dual_pairs.dual_pair_import import dual_pair_import; print(dual_pair_import(\"$<\").frobenius_traces(200))" | sed -e s,\(,[,g -e s,\),],g > $@

.PHONY: traces

traces:
	$(MAKE) `ls GL2-*.gp | sed -e s/.gp/-traces.txt/g`
