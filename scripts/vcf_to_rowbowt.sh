#!/bin/bash

FASTA=$1
VCF=$2
SAMPLES=$3
OUT=${4:-out}
WSIZE=${5:-10}

vcf_to_bwt.py --clean --keep_parse -r -m -S $SAMPLES -o ${OUT} --ma_wsize ${WSIZE} ${FASTA} ${VCF}
rb_build -s -m -o $OUT $OUT
