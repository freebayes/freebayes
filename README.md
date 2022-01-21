# *freebayes*, a haplotype-based variant detector
## user manual and guide


[![Github-CI](https://github.com/freebayes/freebayes/workflows/CI/badge.svg)](https://github.com/freebayes/freebayes/actions) [![Travis-CI](https://travis-ci.com/freebayes/freebayes.svg?branch=master)](https://travis-ci.com/freebayes/freebayes) [![AnacondaBadge](https://anaconda.org/bioconda/freebayes/badges/installer/conda.svg)](https://anaconda.org/bioconda/freebayes) [![DL](https://anaconda.org/bioconda/freebayes/badges/downloads.svg)](https://anaconda.org/bioconda/freebayes) [![BrewBadge](https://img.shields.io/badge/%F0%9F%8D%BAbrew-freebayes-brightgreen.svg)](https://github.com/brewsci/homebrew-bio) [![GuixBadge](https://img.shields.io/badge/gnuguix-freebayes-brightgreen.svg)](https://www.gnu.org/software/guix/packages/F/) [![DebianBadge](https://badges.debian.net/badges/debian/testing/freebayes/version.svg)](https://packages.debian.org/testing/freebayes) [![Chat on Matrix](https://matrix.to/img/matrix-badge.svg)](https://matrix.to/#/#vcflib:matrix.org)
--------

## Overview

[*freebayes*](http://arxiv.org/abs/1207.3907) is a
[Bayesian](http://en.wikipedia.org/wiki/Bayesian_inference) genetic variant
detector designed to find small polymorphisms, specifically SNPs
(single-nucleotide polymorphisms), indels (insertions and deletions), MNPs
(multi-nucleotide polymorphisms), and complex events (composite insertion and
substitution events) smaller than the length of a short-read sequencing
alignment.

*freebayes* is haplotype-based, in the sense that it calls variants based on
the literal sequences of reads aligned to a particular target, not their
precise alignment.  This model is a straightforward generalization of previous
ones (e.g. PolyBayes, samtools, GATK) which detect or report variants based on
alignments.  This method avoids one of the core problems with alignment-based
variant detection--- that identical sequences may have multiple possible
alignments:

<img src="https://github.com/freebayes/freebayes/raw/v1.3.0/paper/haplotype_calling.png" width=500/>

*freebayes* uses short-read alignments
([BAM](http://samtools.sourceforge.net/SAMv1.pdf) files with
[Phred+33](http://en.wikipedia.org/wiki/Phred_quality_score) encoded quality
scores, now standard) for any number of individuals from a population and a
[reference genome](http://en.wikipedia.org/wiki/Reference_genome) (in
[FASTA](http://en.wikipedia.org/wiki/FASTA_format) format)
to determine the most-likely combination of genotypes for the population at
each position in the reference.  It reports positions which it finds putatively
polymorphic in variant call file ([VCF](http://www.1000genomes.org/node/101))
format.  It can also use an input set of variants (VCF) as a source of prior
information, and a copy number variant map (BED) to define non-uniform ploidy
variation across the samples under analysis.

freebayes is maintained by Erik Garrison and Pjotr Prins. See also [RELEASE-NOTES](./RELEASE-NOTES.md).

## Citing freebayes

A preprint [Haplotype-based variant detection from short-read sequencing](http://arxiv.org/abs/1207.3907) provides an overview of the
statistical models used in freebayes.
We ask that you cite this paper if you use freebayes in work that leads to publication.
This preprint is used for documentation and citation.
freebayes was never submitted for review, but has been used in over 1000 publications.

Please use this citation format:

Garrison E, Marth G. Haplotype-based variant detection from short-read sequencing. *arXiv preprint arXiv:1207.3907 [q-bio.GN]* 2012

If possible, please also refer to the version number provided by freebayes when it is run without arguments or with the `--help` option.

## Install

freebayes is provided as a pre-built 64-bit static Linux binary as part of [releases](https://github.com/freebayes/freebayes/releases).

Debian and Conda packages should work too, see the badges at the top
of this page.

To build freebayes from source check the
[development](#Development) section below. It is important to get the full recursive
git checkout and dependencies.

## Support

Please report any issues or questions to the [freebayes mailing list](https://groups.google.com/forum/#!forum/freebayes). Report bugs on the [freebayes issue tracker](https://github.com/freebayes/freebayes/issues)

## Usage

In its simplest operation, freebayes requires only two inputs: a FASTA reference sequence, and a BAM-format alignment file sorted by reference position.
For instance:

    freebayes -f ref.fa aln.bam >var.vcf

... will produce a VCF file describing all SNPs, INDELs, and haplotype variants between the reference and aln.bam. The CRAM version is

    freebayes -f ref.fa aln.cram >var.vcf

Multiple BAM files may be given for joint calling.

Typically, we might consider two additional parameters.
GVCF output allows us to have coverage information about non-called sites, and we can enable it with `--gvcf`.
For performance reasons we may want to skip regions of extremely high coverage in the reference using the `--skip-coverage` parameter or `-g`.
These can greatly increase runtime but do not produce meaningful results.
For instance, if we wanted to exclude regions of 1000X coverage, we would run:

    freebayes -f ref.fa aln.bam --gvcf -g 1000 >var.vcf

For a description of available command-line options and their defaults, run:

    freebayes --help

## Examples

Call variants assuming a diploid sample:

    freebayes -f ref.fa aln.bam >var.vcf

Call variants on only chrQ:

    freebayes -f ref.fa -r chrQ aln.bam >var.vcf

Call variants on only chrQ, from position 1000 to 2000:

    freebayes -f ref.fa -r chrQ:1000-2000 aln.bam >var.vcf

Require at least 5 supporting observations to consider a variant:

    freebayes -f ref.fa -C 5 aln.bam >var.vcf

Skip over regions of high depth by discarding alignments overlapping positions where total read depth is greater than 200:

    freebayes -f ref.fa -g 200 aln.bam >var.vcf

Use a different ploidy:

    freebayes -f ref.fa -p 4 aln.bam >var.vcf

Assume a pooled sample with a known number of genome copies.  Note that this
means that each sample identified in the BAM file is assumed to have 32 genome
copies.  When running with high --ploidy settings, it may be required to set
`--use-best-n-alleles` to a low number to limit memory usage.

    freebayes -f ref.fa -p 32 --use-best-n-alleles 4 --pooled-discrete aln.bam >var.vcf

Generate frequency-based calls for all variants passing input thresholds. You'd do
this in the case that you didn't know the number of samples in the pool.

    freebayes -f ref.fa -F 0.01 -C 1 --pooled-continuous aln.bam >var.vcf

Use an input VCF (bgzipped + tabix indexed) to force calls at particular alleles:

    freebayes -f ref.fa -@ in.vcf.gz aln.bam >var.vcf

Generate long haplotype calls over known variants:

    freebayes -f ref.fa --haplotype-basis-alleles in.vcf.gz \
                        --haplotype-length 50 aln.bam

Naive variant calling: simply annotate observation counts of SNPs and indels:

    freebayes -f ref.fa --haplotype-length 0 --min-alternate-count 1 \
        --min-alternate-fraction 0 --pooled-continuous --report-monomorphic >var.vcf

## Parallelisation

In general, freebayes can be parallelised by running multiple instances of freebayes on separate regions of the genome, and then concatenating the resulting output.
The wrapper, [freebayes-parallel](https://github.com/ekg/freebayes/blob/master/scripts/freebayes-parallel) will perform this, using [GNU parallel](https://www.gnu.org/software/parallel/).

Example freebayes-parallel operation (use 36 cores in this case):

    freebayes-parallel <(fasta_generate_regions.py ref.fa.fai 100000) 36 \
        -f ref.fa aln.bam > var.vcf

Note that any of the above examples can be made parallel by using the
scripts/freebayes-parallel script.  If you find freebayes to be slow, you
should probably be running it in parallel using this script to run on a single
host, or generating a series of scripts, one per region, and run them on a
cluster. Be aware that the freebayes-parallel script contains calls to other programs using relative paths from the scripts subdirectory; the easiest way to ensure a successful run is to invoke the freebayes-parallel script from within the scripts subdirectory.

A current limitation of the freebayes-parallel wrapper, is that due to variance in job memory and runtimes, some cores can go unused for long periods, as they will not move onto the next job unless all cores in use have completed their respective genome chunk. This can be partly avoided by calculating coverage of the input bam file, and splitting the genome into regions of equal coverage using the [coverage_to_regions.py script](https://github.com/freebayes/freebayes/blob/master/scripts/coverage_to_regions.py). An alternative script [split_ref_by_bai_datasize.py](https://github.com/freebayes/freebayes/blob/master/scripts/split_ref_by_bai_datasize.py) will determine target regions based on the data within multiple bam files, with the option of choosing a target data size. This is useful when submitting to Slurm and other cluster job managers, where use of resources needs to be controlled.

Alternatively, users may wish to parallelise freebayes within the workflow manager [snakemake](https://snakemake.readthedocs.io/en/stable/). As snakemake automatically dispatches jobs when a core becomes available, this avoids the above issue. An example [.smk file](https://github.com/freebayes/freebayes/blob/master/examples/snakemake-freebayes-parallel.smk), and associated [conda environment recipe](https://github.com/freebayes/freebayes/blob/master/examples/freebayes-env.yaml), can be found in the /examples directory.

## Calling variants: from fastq to VCF

You've sequenced some samples.  You have a reference genome or assembled set of
contigs, and you'd like to determine reference-relative variants in your
samples.  You can use freebayes to detect the variants, following these steps:

* **Align** your reads to a suitable reference (e.g. with
[bwa](http://bio-bwa.sourceforge.net/) or
[MOSAIK](https://github.com/wanpinglee/MOSAIK))
* Ensure your alignments have **read groups** attached so their sample may be
identified by freebayes.  Aligners allow you to do this, but you can also use
[bamaddrg](http://github.com/ekg/bamaddrg) to do so post-alignment.
* **Sort** the alignments (e.g. [sambamba sort](https://github.com/biod/sambamba)).
* **Mark duplicates**, for instance with [sambamba markdup](https://github.com/biod/sambamba) (if PCR was used in the preparation of your sequencing library)
* ***Run freebayes*** on all your alignment data simultaneously, generating a
VCF.  The default settings should work for most use cases, but if your samples
are not diploid, set the `--ploidy` and adjust the `--min-alternate-fraction`
suitably.
* **Filter** the output e.g. using reported QUAL and/or depth (DP) or
observation count (AO).
* **Interpret** your results.
* (possibly, **Iterate** the variant detection process in response to insight
gained from your interpretation)

freebayes emits a standard VCF 4.1 output stream.  This format is designed for the
probabilistic description of allelic variants within a population of samples,
but it is equally suited to describing the probability of variation in a single
sample.

Of primary interest to most users is the QUAL field, which estimates the
probability that there is a polymorphism at the loci described by the record.
In freebayes, this value can be understood as 1 - P(locus is homozygous given
the data).  It is recommended that users use this value to filter their
results, rather than accepting anything output by freebayes as ground truth.

By default, records are output even if they have very low probability of
variation, in expectation that the VCF will be filtered using tools such as
[vcffilter](http://github.com/ekg/vcflib#vcffilter) in
[vcflib](http://github.com/ekg/vcflib), which is also included in the
repository under `vcflib/`.  For instance,

    freebayes -f ref.fa aln.bam | vcffilter -f "QUAL > 20" >results.vcf

removes any sites with estimated probability of not being polymorphic less than
phred 20 (aka 0.01), or probability of polymorphism &gt; 0.99.

In simulation, the [receiver-operator
characteristic](https://en.wikipedia.org/wiki/Receiver_operating_characteristic)
 (ROC) tends to have a very sharp inflection between Q1 and Q30, depending on
input data characteristics, and a filter setting in this range should provide
decent performance.  Users are encouraged to examine their output and both
variants which are retained and those they filter out.  Most problems tend to
occur in low-depth areas, and so users may wish to remove these as well, which
can also be done by filtering on the DP flag.


## Calling variants in a population

freebayes is designed to be run on many individuals from the same population
(e.g. many human individuals) simultaneously.  The algorithm exploits a neutral
model of allele diffusion to impute most-confident genotypings
across the entire population.  In practice, the discriminant power of the
method will improve if you run multiple samples simultaneously.  In other
words, if your
study has multiple individuals, you should run freebayes against them at the
same time.  This also ensures consistent reporting of information about
evidence for all samples at any locus where any are apparently polymorphic.

To call variants in a population of samples, each alignment must have a read
group identifier attached to it (RG tag), and the header of the BAM file in
which it resides must map the RG tags to sample names (SM).  Furthermore, read
group IDs must be unique across all the files used in the analysis.  One read
group cannot map to multiple samples.  The reason this is required is that
freebayes operates on a virtually merged BAM stream provided by the BamTools
API.  If merging the files in your analysis using bamtools merge would generate
a file in which multiple samples map to the same RG, the files are not suitable
for use in population calling, and they must be modified.

Users may add RG tags to BAM files which were generated without this
information by using (as mentioned in "Calling variants" above)
[bamaddrg](http://github.com/ekg/bamaddrg).
If you have many files corresponding to
many individuals, add a unique read group and sample name to each, and then
open them all simultaneously with freebayes.  The VCF output will have one
column per sample in the input.


## Performance tuning

If you find freebayes to be slow, or use large amounts of memory, consider the
following options:

- Set `--use-best-n-alleles 4`: this will reduce the number of alleles that are
  considered, which will decrease runtime at the cost of sensitivity to
lower-frequency alleles at multiallelic loci.  Calculating site qualities
requires O(samples\*genotypes) runtime, and the number of genotypes is
exponential in ploidy and the number of alleles that are considered, so this is
very important when working with high ploidy samples (and also
`--pooled-discrete`). By default, freebayes puts no limit on this.

- Remove `--genotype-qualities`: calculating genotype qualities requires
  O(samples\*genotypes) memory.

- Set higher input thresholds. Require that N reads in one sample support an
  allele in order to consider it: `--min-alternate-count N`, or that the allele
fraction in one sample is M: `--min-alternate-fraction M`. This will filter
noisy alleles.  The defaults, `--min-alternate-count 2 --min-alternate-fraction
0.2`, are most-suitable for diploid, moderate-to-high depth samples, and should
be changed when working with different ploidy samples. Alternatively,
`--min-alternate-qsum` can be used to set a specific quality sum, which may be
more flexible than setting a hard count on the number of observations.


## Observation filters and qualities

### Input filters

By default, freebayes doesn't

freebayes may be configured to filter its input so as to ignore low-confidence alignments and alleles which are only supported by low-quality sequencing observations (see `--min-mapping-quality` and `--min-base-quality`).
It also will only evaluate a position if at least one read has mapping quality of `--min-supporting-mapping-quality` and one allele has quality of at least `--min-supporting-base-quality`.

Reads with more than a fixed number of high-quality mismatches can be excluded by specifying `--read-mismatch-limit`.
This is meant as a workaround when mapping quality estimates are not appropriately calibrated.

Reads marked as duplicates in the BAM file are ignored, but this can be disabled for testing purposes by providing `--use-duplicate-reads`.
freebayes does not mark duplicates on its own, you must use another process to do this, such as that in [sambamba](https://github.com/biod/sambamba).

### Observation thresholds

As a guard against spurious variation caused by sequencing artifacts, positions are skipped when no more than `--min-alternate-count` or `--min-alternate-fraction` non-clonal observations of an alternate are found in one sample.
These default to 2 and 0.05 respectively.
The default setting of `--min-alternate-fraction 0.05` is suitable for diploid samples but may need to be changed for higher ploidy.

### Allele type exclusion
freebayes provides a few methods to ignore certain classes of allele, e.g.
`--throw-away-indels-obs` and `--throw-awary-mnps-obs`.  Users are *strongly cautioned against using
these*, because removing this information is very likely to reduce detection
power.  To generate a report only including SNPs, use vcffilter post-call as
such:

    freebayes ... | vcffilter -f "TYPE = snp"

### Normalizing variant representation

If you wish to obtain a VCF that does not contain haplotype calls or complex alleles, first call with default parameters and then decompose the output with tools in vcflib, vt, vcf-tools, bcftools, GATK, or Picard.
Here we use a tool in vcflib that normalizes the haplotype calls into pointwise SNPs and indels:

    freebayes ... | vcfallelicprimitives -kg >calls.vcf

Note that this is not done by default as it makes it difficult to determine which variant calls freebayes completed.
The raw output faithfully describes exactly the calls that were made.

### Observation qualities

freebayes estimates observation quality using several simple heuristics based
on manipulations of the phred-scaled base qualities:

* For single-base observations, *mismatches* and *reference observations*: the
un-adjusted base quality provided in the BAM alignment record.
* For *insertions*: the mean quality of the bases inside of the putatively
inserted sequence.
* For *deletions*: the mean quality of the bases flanking the putatively
deleted sequence.
* For *haplotypes*: the mean quality of allele observations within the
haplotype.

By default, both base and mapping quality are into the reported site quality (QUAL in the VCF) and genotype quality (GQ, when supplying `--genotype-qualities`).
This integration is driven by the "Effective Base Depth" metric first developed in [snpTools](http://www.hgsc.bcm.edu/software/snptools), which scales observation quality by mapping quality: *P(Obs|Genotype) ~ P(MappedCorrectly(Obs))P(SequencedCorrectly(Obs))*.
Set `--standard-gls` to use the model described in the freebayes preprint.

## Stream processing

freebayes can read BAM from standard input `--stdin` instead of directly from
files.  This allows the application of any number of streaming BAM filters and
calibrators to its input.

    bam_merger.sh | streaming_filter_or_process.sh | freebayes --stdin ...

This pattern allows the adjustment of alignments without rewriting BAM files,
which could be expensive depending on context and available storage.  A prime
example of this would be graph-based realignment of reads to known variants as
implemented in [glia](http://github.com/ekg/glia).

Using this pattern, you can filter out reads with certain criteria using
bamtools filter without having to modify the input BAM file.  You can also use
the bamtools API to write your own custom filters in C++.  An example filter is
bamfiltertech
[src/bamfiltertech.cpp](http://github.com/freebayes/freebayes/blob/master/src/bamfilte
rtech.cpp), which could be used to filter out
technologies which have characteristic errors which may frustrate certain types
of variant detection.

## INDELs

In principle, any gapped aligner which is sensitive to indels will
produce satisfactory input for use by freebayes.  Due to potential ambiguity,
indels are
not parsed when they overlap the beginning or end of alignment boundaries.

When calling indels, it is important to homogenize the positional distribution
of insertions and deletions in the input by using left realignment.  This is
now done automatically by freebayes, but the behavior can be turned off via
`--dont-left-align-indels` flag.  You probably don't want to do this.

Left realignment will place all indels in homopolymer and microsatellite
repeats at the same position, provided that doing so does not introduce
mismatches between the read and reference other than the indel.  This method
computationally inexpensive and handles the most common classes of alignment
inconsistency.

## Haplotype calls

As freebayes is haplotype-based, left-alignment is necessary only for the
determination of candidate polymorphic loci.  Once such loci are determined,
haplotype observations are extracted from reads where:

1. putative variants lie within `--haplotype-length` bases of each other
(default 3bp),
2. the reference sequence has repeats (e.g. microsatellites or STRs are called
as one haplotype),
3. the haplotype which is called has Shannon entropy less than
`--min-repeat-entropy`, which is off by default but can be set to ~1 for
optimal genotyping of indels in lower-complexity sequence.

After a haplotype window is determined by greedily expanding the window across
overlapping haplotype observations, all reads overlapping the window are used
to establish data likelihoods, *P(Observations|Genotype)*, for all haplotypes
which have sufficient support to pass the input filters.

Partial observations are considered to support those haplotypes which they
could match exactly.  For expedience, only haplotypes which are contiguously
observed by the reads are considered as putative alleles in this process.  This
differs from other haplotype-based methods, such as
[Platypus](http://www.well.ox.ac.uk/platypus), which consider all possible
haplotypes composed of observed component alleles (SNPs, indels) in a given
region when generating likelihoods.

The primary adantages of this approach are conceptual simplicity and
performance, and it is primarily limited in the case of short reads, an issue
that is mitigated by increasing read lengths.  Also, a hybrid approach must be
used to call haplotypes from high-error rate long reads.

### Re-genotyping known variants and calling long haplotypes

For longer reads with higher error rates, it is possible to generate long
haplotypes in two passes over the data.  For instance, if we had very long
reads (e.g. >10kb) at moderate depth and high error rate (>5%) such as might be
produced by PacBio, we could do something like:

    freebayes -f ref.fa aln.bam | vcffilter -f "QUAL > 20" >vars.vcf

... thus generating candidate variants of suitable quality using the default
detection window.  We can then use these as "basis alleles" for the observation
of haplotypes, considering all other putative variants supported by the
alignment to be sequencing errors:

    freebayes -f ref.fa --haplotype-length 500 \
        --haplotype-basis-alleles vars.vcf aln.bam >haps.vcf

These steps should allow us to read long haplotypes directly from input data
with high error rates.

The high error rate means that beyond a small window each read will contain a
completely different literal haplotype.  To a point, this property improves our
signal to noise ratio and can effectively filter out sequencing errors at the
point of the input filters, but it also decreases the effective observation
depth will prevent the generation of any calls if a long `--haplotype-length`
is combined with high a sequencing error rate.


## Best practices and design philosophy

freebayes follows the patterns suggested by the [Unix philosophy](https://en.wikipedia.org/wiki/Unix_philosophy), which promotes the development of simple, modular systems that perform a single function, and can be combined into more complex systems using stream processing of common interchange formats.

freebayes incorporates a number of features in order to reduce the complexity of variant detection for researchers and developers:

* **Indel realignment is accomplished internally** using a read-independent method, and issues resulting from discordant alignments are dramatically reducedy through the direct detection of haplotypes.
* The need for **base quality recalibration is avoided** through the direct detection of haplotypes. Sequencing platform errors tend to cluster (e.g. at the ends of reads), and generate unique, non-repeating haplotypes at a given locus.
* **Variant quality recalibration is avoided** by incorporating a number of metrics, such as read placement bias and allele balance, directly into the Bayesian model.  (Our upcoming publication will discuss this in more detail.)

A minimal pre-processing pipeline similar to that described in "Calling variants" should be sufficient for most uses.
For more information, please refer to a post by Brad Chapman [on minimal BAM preprocessing methods](http://bcbio.wordpress.com/2013/10/21/updated-comparison-of-variant-detection-methods-ensemble-freebayes-and-minimal-bam-preparation-pipelines/).

## Development

To download freebayes, please use git to download the most recent development tree:

    git clone --recursive https://github.com/freebayes/freebayes.git

If you have a repo, update the submodules with

    git submodule update --init --recursive --progress

On Debian you'll need a gcc compiler and want packages:

- bc
- samtools
- parallel
- meson
- ninja-build
- libvcflib-tools
- vcftools

Build dependencies are listed in [guix.scm](./guix.scm) and
[travis](.travis.yml). Builds have been tested with gcc 7 and clang 9.

## Compilation

Make sure to have dependencies installed and checkout the tree
with `--recursive`.

Freebayes can target AMD64 and ARM64 (with Neon extensions).

Recently we added the meson build system which can be run with

    meson build/ --buildtype debug

or to setup with clang instead

    env CXX=clang++ CC=clang CC_LD=lld meson build --buildtype debug

Next compile and test in the ~build~ directory

    cd build
    ninja
    ninja test

The freebayes binary should be in

    build/freebayes

Tests on ARM may be slow. If you get a TIMEOUT use a multiplier,
e.g.

    meson test -t 4 -C build/

See [meson.build](./meson.build) for more information.

### Compile in a Guix container

After checking out the repo with git recursive create a Guix
container with all the build tools with

    guix shell -C -D -f guix.scm

See also [guix.scm](./guix.scm).
