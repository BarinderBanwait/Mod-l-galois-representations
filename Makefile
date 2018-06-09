# Generate LMFDB data for Galois representations attached to modular
# forms over finite fields

# Input for these scripts: the files mod-L-reps.txt for L = 2, 3, 5
# containing eigenforms computed by Andrew Sutherland (not part of
# this repository).

GP ?= gp -q

forms.gp: mod-2-reps.gp mod-3-reps.gp mod-5-reps.gp
	> $@
	for l in 2 3 5; do echo "convert(readvec(\"mod-$$l-reps.gp\"),$$l);" | $(GP) convert.gp >> $@; done

mod-%-reps.gp: mod-%-reps.txt
	grep '^<' $< | sed -e 's/<\[/[/' -e 's/]>/]/' -e 's/</[/g' -e 's/>/]/g' -e 's/,$$//' > $@
