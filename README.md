[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/scallop-lr/README.html)

# Overview
Scallop-LR is a reference-based transcript assembler for long-reads RNA-seq data.
Its method has been described in the manuscript at [bioRxiv](https://doi.org/10.1101/632703).
The datasets and scripts used in the manuscript to illustrate its performance
are available at (https://github.com/Kingsford-Group/lrassemblyanalysis).

# Release
Latest release of Scallop-LR is [v0.9.2](https://github.com/Kingsford-Group/scallop/releases/tag/isoseq-v0.9.2).

# Installation
Download the source code of Scallop-LR from
[here](https://github.com/Kingsford-Group/scallop/releases/download/isoseq-v0.9.2/scallop-lr-0.9.2.tar.gz).
Scallop-LR uses additional libraries of Boost and htslib (*NOTE:* from v0.9.2 the dependence on Clp is optional). 
If they have not been installed in your system, you first
need to download and install them. You might also need to
export the runtime library path to certain environmental
variable (for example, `LD_LIBRARY_PATH`, for most linux distributions).
After install these dependencies, you then compile the source code of Scallop-LR.
If some of the above dependencies are not installed to the default system 
directories (for example, `/usr/local`, for most linux distributions),
their corresponding installing paths should be specified to `configure` of Scallop-LR.

## Download Boost
If Boost has not been downloaded/installed, download Boost
[(license)](http://www.boost.org/LICENSE_1_0.txt) from (http://www.boost.org).
Uncompress it somewhere (compiling and installing are not necessary).

## Install htslib
If htslib has not been installed, download htslib 
[(license)](https://github.com/samtools/htslib/blob/develop/LICENSE)
from (http://www.htslib.org/) with version 1.5 or higher.
Note that htslib relies on zlib. So if zlib has not been installed in your system,
you need to install zlib first. To do so, download zlib
[(license)](https://zlib.net/zlib_license.html) at (https://zlib.net/).
Use the following commands to install zlib:
```
./configure
make
make install
```
After installing zlib, use the following commands to build htslib:
```
./configure --disable-bz2 --disable-lzma --disable-gcs --disable-s3 --enable-libcurl=no
make
make install
```
The default installation location of htslib is `/usr/lib`.
If you would install it to a different location, replace the above `configure` line with
the following (by adding `--prefix=/path/to/your/htslib` to the end):
```
./configure --disable-bz2 --disable-lzma --disable-gcs --disable-s3 --enable-libcurl=no --prefix=/path/to/your/htslib
```
In this case, you also need to export the runtime library path (note that there
is an additional `lib` following the installation path):
```
export LD_LIBRARY_PATH=/path/to/your/htslib/lib:$LD_LIBRARY_PATH
```
## Install Clp (optional since v0.9.2)
*NOTE:* Clp will be used to solve the linear programming instances
created when decomposing unsplitable vertices. An alternative algorithm
is provided in Scallop-LR from version v0.9.2 (and hence since then the installation of Clp becomes optional).
Our testing shows that these two algorithms
give very similar results.

If Clp has not been installed in your system, 
download Clp [(license)](https://opensource.org/licenses/eclipse-1.0)
from (https://projects.coin-or.org/Clp). 
Use the following to install Clp
```
./configure --disable-bzlib --disable-zlib
make
make install
```
The default installation of Clp is the current directory, rather than `/usr/lib`.
If you would install it to a different location, replace the above `configure` line with
the following (by adding `--prefix=/path/to/your/Clp` to the end):
```
./configure --disable-bzlib --disable-zlib --prefix=/path/to/your/Clp
```
You need to export the runtime library path (note that there
is an additional `lib` following the installation path):
```
export LD_LIBRARY_PATH=/path/to/your/Clp/lib:$LD_LIBRARY_PATH
```

## Compile Scallop-LR

Use the following to compile Scallop-LR (without Clp; therefore the alternative algorithm for decomposing unsplitable vertices will be used; available
for versions newer than v0.9.2):
```
./configure --with-htslib=/path/to/your/htslib --with-boost=/path/to/your/boost
make
```
Use the following to compile Scallop-LR (with Clp; therefore an linear programming formulation will be used to decompose unsplitable vertices):
```
./configure --with-htslib=/path/to/your/htslib --with-boost=/path/to/your/boost --enable-useclp --with-clp=/path/to/your/Clp
make
```

If some of the dependencies are installed in the default system directory (for example, `/usr/lib`),
then the corresponding `--with-` option might not be necessary.
The executable file `scallop-lr` will appear at `src/scallop-lr`.


# Usage

The usage of `scallop` is:
```
./scallop-lr -i <input.bam> -o <output.gtf> -c <ccs-header-file> [options]
```

The `input.bam` is the read alignment file generated by some (long-reads) RNA-seq aligner such as minimap2 or GMAP.
Make sure that it is sorted; otherwise run `samtools` to sort it:
```
samtools sort input.bam > input.sort.bam
```

The assembled transcripts shall be written as gtf format into `output.gtf`.

`ccs-header-file` collects the header lines (i.e., lines started with `>`) of the `.fasta` files
generated with PacBio SMRT pipeline. Scallop-LR uses the information in the header lines
to infer the start/end boundaries of the transcripts.
A typical `ccs-header-file` looks like:

```
>m151676_s1_p0/4994/28_2528_CCS strand=+;fiveseen=1;polyAseen=1;threeseen=1;fiveend=28;polyAend=2528;threeend=2552;primer=1;chimera=0
>m151676_s1_p0/4996/28_2690_CCS strand=+;fiveseen=1;polyAseen=1;threeseen=1;fiveend=28;polyAend=2690;threeend=2714;primer=1;chimera=0
>m151676_s1_p0/5003/30_2471_CCS strand=+;fiveseen=1;polyAseen=1;threeseen=1;fiveend=30;polyAend=2471;threeend=2502;primer=1;chimera=0
...
```


Scallop-LR support the following parameters. 
Please refer to the additional explanation below the table.

 Parameters | Default Value | Description
 ------------------------- | ------------- | ----------
 --help  | | print usage of Scallop-LR and exit
 --version | | print version of Scallop-LR and exit
 --preview | | show the inferred `library_type` and exit
 --verbose | 1 | chosen from {0, 1, 2}
 --library_type               | empty | chosen from {empty, unstranded, first, second}
 --min_transcript_coverage    | 1 | the minimum coverage required to output a multi-exon transcript
 --min_single_exon_coverage   | 10 | the minimum coverage required to output a single-exon transcript
 --min_transcript_length_base      |100 | the minimum base length of a transcript
 --min_transcript_length_increase  | 50 | the minimum increased length of a transcript with each additional exon
 --min_mapping_quality        | 1 | ignore reads with mapping quality less than this value
 --max_num_cigar              | 10000 | ignore reads with CIGAR size larger than this value
 --min_bundle_gap             | 50 | the minimum distances required to start a new bundle
 --min_num_hits_in_bundle     | 1 | the minimum number of reads required in a bundle
 --min_flank_length           | 3 | the minimum match length required in each side for a spliced read
 --min_splice_bundary_hits    | 1 | the minimum number of spliced reads required to support a junction

1. For `--verbose`, 0: quiet; 1: one line for each splice graph; 2: details of graph decomposition.

2. `--library_type` is highly recommended to provide. The `unstranded`, `first`, and `second`
correspond to `fr-unstranded`, `fr-firststrand`, and `fr-secondstrand` used in standard Illumina
sequencing libraries. If none of them is given, i.e., it is `empty` by default, then Scallop-LR
will try to infer the `library_type` by itself (see `--preview`). Notice that such inference is based
on the `XS` tag stored in the input `bam` file. If the input `bam` file do not contain `XS` tag,
then it is essential to provide the `library_type` to Scallop-LR. You can try `--preview` to see
the inferred `library_type`.

3. `--min_transcript_coverage` is used to filter lowly expressed transcripts: Scallop-LR will filter
out transcripts whose (predicted) raw counts (number of moleculars) is less than this number.

4. `--min_transcript_length_base` and `--min_transcript_length_increase` is combined to filter
short transcripts: the minimum length of a transcript is given by `--min_transcript_length_base`
\+ `--min_transcript_length_increase` * num-of-exons-in-this-transcript. Transcripts that are less
than this number will be filtered out.


