[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/scallop/README.html)

# Overview
Scallop is an accurate reference-based transcript assembler. Scallop features
its high accuracy in assembling multi-exon transcripts as well as lowly
expressed transcripts. Scallop achieves this improvement through a novel
algorithm that can be proved preserving all phasing paths from reads and paired-end reads,
while also achieves both transcripts parsimony and coverage deviation minimization.

Scallop paper has been published at [**Nature Biotechnology**](https://www.nature.com/articles/nbt.4020).
The datasets and scripts used in this paper to compare the performance of Scallop
and other assemblers are available at [**scalloptest**](https://github.com/Kingsford-Group/scalloptest).

Please also checkout the **podcast** about Scallop (thanks [Roman Cheplyaka](https://ro-che.info/) for the interview).
It is available at both [the bioinformatics chat](https://bioinformatics.chat/scallop) and
[iTunes](https://itunes.apple.com/us/podcast/the-bioinformatics-chat/id1227281398). 

# Release
Latest release of Scallop is [v0.10.4](https://github.com/Kingsford-Group/scallop/releases/tag/v0.10.4),
including binary 
(for both [linux](https://github.com/Kingsford-Group/scallop/releases/download/v0.10.4/scallop-0.10.4_linux_x86_64.tar.gz)
and [mac](https://github.com/Kingsford-Group/scallop/releases/download/v0.10.4/scallop-0.10.4_macOS_10.14.tar.gz))
and [source code](https://github.com/Kingsford-Group/scallop/releases/download/v0.10.4/scallop-0.10.4.tar.gz).

Below  we list the systems that have been tested for whether the Scallop binary can run or not.

 Operation System | Version | Code Name | Scallop
 ---------------- | ------- | --------- | ----------
 Debian | 9		| Stretch | [linux](https://github.com/Kingsford-Group/scallop/releases/download/v0.10.4/scallop-0.10.4_linux_x86_64.tar.gz)
 Ubuntu | 14.04 | Trusty Tahr | [linux](https://github.com/Kingsford-Group/scallop/releases/download/v0.10.4/scallop-0.10.4_linux_x86_64.tar.gz)
 Ubuntu | 16.04 | Xenial Xerus | [linux](https://github.com/Kingsford-Group/scallop/releases/download/v0.10.4/scallop-0.10.4_linux_x86_64.tar.gz)
 CentOS | 6.9   | | N/A
 CentOS | 7     | | [linux](https://github.com/Kingsford-Group/scallop/releases/download/v0.10.4/scallop-0.10.4_linux_x86_64.tar.gz)
 Fedora | 24    | | [linux](https://github.com/Kingsford-Group/scallop/releases/download/v0.10.4/scallop-0.10.4_linux_x86_64.tar.gz)
 Mac OS | 10.10 | Yosemite | [mac](https://github.com/Kingsford-Group/scallop/releases/download/v0.10.4/scallop-0.10.4_macOS_10.10.tar.gz)
 Mac OS | 10.11 | El Capitan | [mac](https://github.com/Kingsford-Group/scallop/releases/download/v0.10.4/scallop-0.10.4_macOS_10.10.tar.gz)
 Mac OS | 10.12 | Sierra | [mac](https://github.com/Kingsford-Group/scallop/releases/download/v0.10.4/scallop-0.10.4_macOS_10.10.tar.gz)

# Support
Scallop is, and will continue to be, [freely and actively supported on a
best-effort basis](https://oceangenomics.com/about/#open).

If you need industrial-grade technical support, please consider the options at
[oceangenomics.com/support](http://oceangenomics.com/support).

# Installation
Download the source code of Scallop from
[here](https://github.com/Kingsford-Group/scallop/releases/download/v0.10.4/scallop-0.10.4.tar.gz).
Scallop uses additional libraries of Boost and htslib (*NOTE:* from v0.10.4 the dependence on Clp is optional). 
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
## Install Clp (optional since v0.10.4)
*NOTE:* Clp will be used to solve the linear programming instances
created when decomposing unsplitable vertices. An alternative algorithm
is provided in Scallop from version v0.10.4~(and hence since then the installation of Clp becomes optional).
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

## Compile Scallop

Use the following to compile Scallop (without Clp; therefore the alternative algorithm for decomposing unsplitable vertices will be used; available
for versions newer than v0.10.4):
```
./configure --with-htslib=/path/to/your/htslib --with-boost=/path/to/your/boost
make
```
Use the following to compile Scallop (with Clp; therefore an linear programming formulation will be used to decompose unsplitable vertices):
```
./configure --with-htslib=/path/to/your/htslib --with-boost=/path/to/your/boost --enable-useclp --with-clp=/path/to/your/Clp
make
```

If some of the dependencies are installed in the default system directory (for example, `/usr/lib`),
then the corresponding `--with-` option might not be necessary.
The executable file `scallop` will appear at `src/scallop`.


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

Scallop support the following parameters. Please refer
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
 --max_num_cigar              | 7 | ignore reads with CIGAR size larger than this value
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


# Quantification by Combining Scallop and Salmon

We recommend users to perform RNA-seq quantification using the combination of Scallop and Salmon.
This pipeline involves the following steps:

**Step 1:** Align the reads to a reference genome (for example, with
[TopHat2](https://ccb.jhu.edu/software/tophat/index.shtml),
[STAR](https://github.com/alexdobin/STAR), or
[HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml))
to obtain the (sorted) reads alignment file `sort.bam`.

**Step 2:** Assemble the expressed transcripts with [Scallop](https://github.com/Kingsford-Group/scallop):
```
scallop -i sort.bam -o scallop.gtf
```
The assembled transcripts will be written to `scallop.gtf`.

**Step 3:** Use [gffcompare](http://ccb.jhu.edu/software/stringtie/gff.shtml) to
evaluate the assembled transcripts using a reference annotation:
```
gffcompare -o gffall -r reference.gtf scallop.gtf
```
where `reference.gtf` is the reference annotation file
(for example, [ensembl annotation](https://goo.gl/cifLXe)).
This command will generate a file `gffall.scallop.gtf.map` defining which transcripts in `scallop.gtf`
can be found in the `reference.gtf`.

**Step 4:** Union the assembled transcripts with the reference transcriptome. Specifically,
First, use our tool
[gtfcuff](https://github.com/Kingsford-Group/rnaseqtools) to fetch the transcripts that
are only in `scallop.gtf`:
```
gtfcuff puniq gffall.scallop.gtf.tmap scallop.gtf reference.gtf unique.gtf
```
The uniquely expressed transcripts (i.e., those are in `scallop.gtf` but not in `reference.gtf`)
will be written to `unique.gtf`.
Second, extract the cDNA sequences of the transcripts in `unique.gtf` 
from a reference genome using tool
[gffread](http://ccb.jhu.edu/software/stringtie/gff.shtml):
```
gffread unique.gtf -g genome -w unique.fa
```
where `genome` is the reference genome, for example 
[ensembl reference genome](https://goo.gl/V1V3B1).
The cDNA sequences of the uniquely assembled transcripts (i.e., those in `unique.gtf`)
will be written to `unique.fa`.
Finally, merge `unique.fa` and the reference transcriptome to obtained the extended transcriptome: 
```
cat unique.fa reference.fa > union.fa
```
where `reference.fa` is the reference transcriptome (i.e., the cDNA sequences of the
transcripts in `reference.gtf`), for example,
[ensembl cDNA sequences](https://goo.gl/rJdrEX).
The extended transcriptome will be written to `union.fa`.


**Step 5:**  Run [Salmon](https://combine-lab.github.io/salmon/)
to quantify with respect to the above extended transcriptome.
First, create Salmon index:
```
salmon index -t union.fa -i salmon.index -p 4
```
After that we can quantify:
```
salmon quant -i salmon.index -1 fastq-file1 -2 fastq-file2 -p 4
```
The main quantification file will appear as `salmon.quant/quant.sf`.
Please refer to [Salmon](https://combine-lab.github.io/salmon/) documentation
for more advanced usage.

