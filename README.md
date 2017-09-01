# Overview
Scallop is an accurate reference-based transcript assembler. Scallop features
its high accuracy in assembling multi-exon transcripts as well as lowly
expressed transcripts. Scallop achieves this improvement through a novel
algorithm that can be proved preserving all phasing paths from reads and paired-end reads,
while also achieves both transcripts parsimony and coverage deviation minimization.

**Pre-print** of Scallop is available at [bioRxiv](http://biorxiv.org/content/early/2017/04/03/123612).
The datasets and scripts used in this paper to compare the performance of Scallop
and other assemblers are available at [**scalloptest**](https://github.com/Kingsford-Group/scalloptest).

Please also checkout the **podcast** about Scallop (thanks [Roman Cheplyaka](https://ro-che.info/) for the interview).
It is available at both [the bioinformatics chat](https://bioinformatics.chat/scallop) and
[iTunes](https://itunes.apple.com/us/podcast/the-bioinformatics-chat/id1227281398). 

# Release
Latest release of Scallop is [v0.10.2](https://github.com/Kingsford-Group/scallop/releases/tag/v0.10.2),
including binary 
(for both [linux](https://github.com/Kingsford-Group/scallop/releases/download/v0.10.2/scallop-0.10.2_linux_x86_64.tar.gz)
and [mac](https://github.com/Kingsford-Group/scallop/releases/download/v0.10.2/scallop-0.10.2_macOS_10.10.tar.gz))
and [source code](https://github.com/Kingsford-Group/scallop/releases/download/v0.10.2/scallop-0.10.2.tar.gz).

# Installation
Download the source code of Scallop from
[here](https://github.com/Kingsford-Group/scallop/releases/download/v0.10.2/scallop-0.10.2.tar.gz).
Scallop uses additional libraries of Boost, htslib and Clp. 
If they have not been installed in your system, you first
need to download and install them. You might also need to
export the runtime library path to certain environmental
variable (for example, `LD_LIBRARY_PATH`, for most linux distributions).
After install these dependencies, you then compile the source code of Scallop.
If some of the above dependencies are not installed to the default system 
directories (for example, `/usr/local`, for most linux distributions),
their corresponding installing paths should be specified to `configure` of Scallop.

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

## Install Clp
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

## Compile Scallop

Use the following to compile Scallop:
```
./configure --with-clp=/path/to/your/Clp --with-htslib=/path/to/your/htslib --with-boost=/path/to/your/boost
make
```
If some of the dependencies are installed in the default system directory (for example, `/usr/lib`),
then the corresponding `--with-` option might not be necessary.
The executable file `scallop` will appear at `src/src/scallop`.


# Usage

The usage of `scallop` is:
```
./scallop -i <input.bam> -o <output.gtf> [options]
```

The `input.bam` is the read alignment file generated by some RNA-seq aligner, (for example, TopHat2, STAR, or HISAT2).
Make sure that it is sorted; otherwise run `samtools` to sort it:
```
samtools sort input.bam > input.sort.bam
```

The reconstructed transcripts shall be written as gtf format into `output.gtf`.

Scallop support the following parameters. Please also refer
to the additional explanation below the table.

 Parameters | Default Value | Description
 ------------------------- | ------------- | ----------
 --help  | | print usage of Scallop and exit
 --version | | print version of Scallop and exit
 --preview | | show the inferred `library_type` and exit
 --verbose | 1 | chosen from {0, 1, 2}
 --library_type               | empty | chosen from {empty, unstranded, first, second}
 --min_transcript_coverage    | 1 | the minimum coverage required to output a multi-exon transcript
 --min_single_exon_coverage   | 20 | the minimum coverage required to output a single-exon transcript
 --min_transcript_length_base      |150 | the minimum base length of a transcript
 --min_transcript_length_increase  | 50 | the minimum increased length of a transcript with each additional exon
 --min_mapping_quality        | 1 | ignore reads with mapping quality less than this value
 --min_bundle_gap             | 50 | the minimum distances required to start a new bundle
 --min_num_hits_in_bundle     | 20 | the minimum number of reads required in a bundle
 --min_flank_length           | 3 | the minimum match length required in each side for a spliced read
 --min_splice_bundary_hits    | 1 | the minimum number of spliced reads required to support a junction

1. For `--verbose`, 0: quiet; 1: one line for each splice graph; 2: details of graph decomposition.

2. `--library_type` is highly recommended to provide. The `unstranded`, `first`, and `second`
correspond to `fr-unstranded`, `fr-firststrand`, and `fr-secondstrand` used in standard Illumina
sequencing libraries. If none of them is given, i.e., it is `empty` by default, then Scallop
will try to infer the `library_type` by itself (see `--preview`). Notice that such inference is based
on the `XS` tag stored in the input `bam` file. If the input `bam` file do not contain `XS` tag,
then it is essential to provide the `library_type` to Scallop. You can try `--preview` to see
the inferred `library_type`.

3. `--min_transcript_coverage` is used to filter lowly expressed transcripts: Scallop will filter
out transcripts whose (predicted) raw counts (number of moleculars) is less than this number.

4. `--min_transcript_length_base` and `--min_transcript_length_increase` is combined to filter
short transcripts: the minimum length of a transcript is given by `--min_transcript_length_base`
\+ `--min_transcript_length_increase` * num-of-exons-in-this-transcript. Transcripts that are less
than this number will be filtered out.
