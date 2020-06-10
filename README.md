# rowbowt: run-length BWT tools for working with genomic sequences

Author: Taher Mun

Code heavily derives from Nicola Prezza's [r-index](https://github.com/nicolaprezza/r-index) repository.

## Compiling

Requires a c++17 compliant compiler and a Unix-based operating system.

To download and build:

```
git clone github.com/alshai/rowbowt;
mkdir rowbowt/build;
cd rowbowt/build;
cmake ..
make
```

This produces the exectuables `rb_align` and `rb_build`.


# Building a rowbowt index

1) generate files needed for building the index (exact format specifications will be added below):
 
    - `.bwt` (required) - contains the BWT of the original text. This will let you find the number of occurences of your
      pattern in the original text

    - `.ssa`, `.esa` (optional) - contains the run-starts and run-ends SA samples, respecitively. Use this if you want
      use the index to find pattern locations in the original text.

    - and `.ma` files - contains the markers (ref position + allele) associated with each SA element. Use this if you
      want to find what alleles wrt the reference sequence your pattern overlaps.

    - we recommend using `pfbwt-f` (included in this repo) to generate these required files from a reference sequence
      and a VCF containing a haplotype reference panel. See the README for `pfbwt-f` for details on how to run.

    ```
    cd pfbwt-f
    make vcf_to_bwt
    python3 vcf_to_bwt.py [-s] [-m] -o <output prefix> <reference fasta> <VCF>
    ```
    
2) Run the build command

    ```
    ./rb_build [-s] [-m] -o <output prefix> <rowbowt index prefix> 
    ```
    
    
# Querying sequences against a rowbowt index

```
./rb_align [-s] [-m] <rowbowt index prefix> <fastq>
```
 
- use `-s` and `-m` to query locations and markers, respectively.

- Proper formatting (ie. sam format for locations, VCF format for markers) will be added in the future.
