#!/bin/sh
FASTA=${1}
OUT=${2:-$FASTA}
~/rowbowt/build/pfbwt-f/pfbwt-f64 --non-acgt-to-a --print-docs -ro ${OUT} ${FASTA};
~/rowbowt/build/rb_build -s -o $OUT $OUT
