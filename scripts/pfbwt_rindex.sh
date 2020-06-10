#!/bin/sh

~/r-index2/pfbwt-f/pfbwt-f64 --non-acgt-to-a --print-docs -ro $1 $1;
~/r-index2/build/ri_build $1;
