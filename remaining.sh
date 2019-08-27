#!/bin/sh

for f in `cut -d: -f1 labels.txt | sort -V`; do
    if [ ! -f $f.gp ]; then
	echo $f "" `grep genus $f/Makefile`
    fi
done
