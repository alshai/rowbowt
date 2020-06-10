r-index: the run-length BWT index
===============


Original Author: Nicola Prezza (nicola.prezza@gmail.com)
Joint work with Travis Gagie and Gonzalo Navarro

Forked by: Taher Mun

Changes: Use of pfBWT as the BWT-construction algorithm (https://gitlab.com/manzai/Big-BWT) and support for indexing+querying DNA sequences. 

cite as:

Gagie T, Navarro G, Prezza N. Optimal-time text indexing in BWT-runs bounded space. In Proceedings of the Twenty-Ninth Annual ACM-SIAM Symposium on Discrete Algorithms, SODA 2018, New Orleans, NA, USA, January 7-10 2017.

Kuhnle A., Mun T., Boucher C., Gagie T., Langmead B., Manzini G. (2019) Efficient Construction of a Complete Index for Pan-Genomics Read Alignment. In: Cowen L. (eds) Research in Computational Molecular Biology. RECOMB 2019. Lecture Notes in Computer Science, vol 11467. Springer, Cham. doi:10.1007/978-3-030-17083-7_10.

### Brief description

The r-index is the first full-text index of size O(r), r being the number of BWT runs of the input text (of size n), supporting fast (almost optimal) locate of pattern occurrences. The r-index employs a novel suffix array sampling of size 2r; in classical FM-indexes, this sampling would result in a locate time of Omega(n/r) per occurrence. The r-index, on the other hand, reduces this time to O(log(n/r)).

Let s be the alphabet size and fix a constant eps>0. The r-index offers the following tradeoffs:

- Space: r * ( log s + (1+eps)log(n/r) + 2log n ) bits
- Count time: O( (m/eps) * (log (n/r) + log s) )
- Locate time: After count, O( log(n/r) ) time per occurrence 

On very repetitive datasets, the r-index locates orders of magnitude faster than the RLCSA (with a sampling rate resulting in the same size for the two indexes).

NEWS: refactored locate strategy. Let (l,r) be the SA range. Now, the index first finds SA[r] and then applies function Phi to locate SA[r-1], SA[r-2], ..., SA[l]. This is both faster and more space efficient than the strategy originally implemented and described in the paper.

### Downloading and Compiling

1) Clone repository

```
git clone --recursive https://github.com/alshai/r-index
```

2) Install required packages:

```
apt-get update -qq && \
apt-get install -y zlib1g-dev git cmake build-essential python3
```

3) compile and install

```
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=<installation dir> ..
make
make install
```

### Usage

To build index from a fasta file (outputs to `input.fa.ri`):

```
ri-buildfasta -b <bigbwt|sais|from_bwt> input.fa
```

To count queries in a fast[a|q] file:

```
ri-align count index_prefix reads.fq
```

To locate queries in a fast[a|q] file:

```
ri-align --max-hits (default:-1) --max-range (default:-1) align index_prefix reads.fq
```
