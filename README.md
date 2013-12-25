# *freebayes*, a haplotype-based variant detector
## user manual and guide

### Erik Garrison <erik.garrison@bc.edu>

--------

## Overview

[*FreeBayes*](http://arxiv.org/abs/1207.3907) is a 
[Bayesian](http://en.wikipedia.org/wiki/Bayesian_inference) genetic variant 
detector designed to find small polymorphisms, specifically SNPs 
(single-nucleotide polymorphisms), indels (insertions and deletions), MNPs 
(multi-nucleotide polymorphisms), and complex events (composite insertion and 
substitution events) smaller than the length of a short-read sequencing 
alignment.

*FreeBayes* is haplotype-based, in the sense that it calls variants based on 
the literal sequences of reads aligned to a particular target, not their 
precise alignment.  This model is a straightforward generalization of previous 
ones (e.g. PolyBayes, samtools, GATK) which detect or report variants based on 
alignments.  This method avoids one of the core problems with alignment-based 
variant detection--- that identical sequences may have multiple possible 
alignments:

<img src="http://hypervolu.me/~erik/freebayes/figures/haplotype_calling.png" 
width=500/>

*FreeBayes* uses short-read alignments 
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


## Citing freebayes

A preprint [Haplotype-based variant detection from short-read 
sequencing](http://arxiv.org/abs/1207.3907) provides an overview of the 
statistical models
used in FreeBayes.  We ask that you cite this paper if you use FreeBayes in
work that leads to publication.

If possible, please also refer to the version number provided by freebayes when 
it is run without arguments or with the `--help` option.  For example, you 
should see something like this:

    version:  v0.9.10-3-g47a713e-dirty

This provides both a point release number and a git commit id, which will 
ensure precise reproduction of results.


## Obtaining

To download FreeBayes, please use git to download the most recent development
tree.  Currently, the tree is hosted on github, and can be obtained via:

    git clone --recursive git://github.com/ekg/freebayes.git

Note the use of --recursive.  This is required, as the project contains some
nested git submodules for external repositories.


## Compilation

FreeBayes requires g++ and the standard C and C++ development libraries.
Additionally, cmake is required for building the BamTools API.

    make

Will build the executable freebayes, as well as the utilities bamfiltertech and 
bamleftalign.  These executables can be found in the `bin/` directory in the 
repository.

Users may wish to install to e.g. /usr/local/bin (default), which is 
accomplished via

    sudo make install


## Usage

In its simplest operation, freebayes requires only two inputs: a FASTA reference
sequence, and a BAM-format alignment file sorted by reference position.  For
instance:

    freebayes --fasta-reference h.sapiens.fasta NA20504.bam

... produce (on standard output) a VCF file on standard out describing
all SNPs, INDELs, MNPs, and Complex events between the reference and the
alignments in NA20504.bam.  In order to produce correct output, the reference
supplied must be the reference to which NA20504.bam was aligned.

Users may specify any number of BAM files on the command line.  FreeBayes uses 
the [BamTools API](http://github.com/pezmaster31/bamtools) to open and parse 
these files in parallel, virtually merging them at runtime into one logical 
file with a merged header.

For a description of available command-line options and their defaults, run:

    freebayes --help


## Calling variants

You've sequenced some samples.  You have a reference genome or assembled set of 
contigs, and you'd like to determine reference-relative variants in your 
samples.  You can use freebayes to detect the variants, following these steps:

* **Align** your reads to a suitable reference (e.g. with 
[bwa](http://bio-bwa.sourceforge.net/) or 
[MOSAIK](https://github.com/wanpinglee/MOSAIK))
* Ensure your alignments have **read groups** attached so their sample may be 
identified by freebayes.  Aligners allow you to do this, but you can also use 
[bamaddrg](http://github.com/ekg/bamaddrg) to do so post-alignment.
* **Sort** the alignments (e.g. bamtools sort).
* **Mark duplicates**, for instance with [samtools 
rmdup](http://samtools.sourceforge.net/) (if PCR was used in the preparation of 
your sequencing library)
* ***Run freebayes*** on all your alignment data simultaneously, generating a 
VCF.  The default settings should work for most use cases, but if your samples 
are not diploid, set the `--ploidy` and adjust the `--min-alternate-fraction` 
suitably.
* **Filter** the output e.g. using reported QUAL and/or depth (DP) or 
observation count (AO).
* **Interpret** your results.
* (possibly, **Iterate** the variant detection process in response to insight 
gained from your interpretation)

FreeBayes a standard VCF 4.1 outut stream.  This format is designed for the
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

FreeBayes is designed to be run on many individuals from the same population
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


## Observation filters and qualities

### Input filters
FreeBayes filters its input so as to ignore low-confidence alignments and
alleles which are only supported by low-quality sequencing observations (see
`--min-mapping-quality` and `--min-base-quality`).  It also will only evaluate a
position if at least one read has mapping quality of
`--min-supporting-mapping-quality` and one allele has quality of at least
`--min-supporting-base-quality`.

Reads with more than a fixed number of high-quality mismatches can be excluded
by specifying `--read-mismatch-limit`.  This is meant as a workaround when 
mapping quality estimates are not appropriately calibrated.

Reads marked as duplicates in the BAM file are ignored, but this can be 
disabled for testing purposes by providing `--use-duplicate-reads`.  FreeBayes 
does not mark duplicates on its own, you must use another process to do this.

### Observation thresholds
As a guard against spurious variation caused by sequencing artifacts, positions
are skipped when no more than `--min-alternate-count` or 
`--min-alternate-fraction`
non-clonal observations of an alternate are found in one sample.  These default 
to 2 and 0.2 respectively.  The default setting of `--min-alternate-fraction 
0.2` is suitable for diploid samples but should be changed for ploidy > 2.

### Allele type exclusion
FreeBayes provides a few methods to ignore certain classes of allele, e.g. 
`--no-indels` and `--no-mnps`.  Users are *strongly cautioned against using 
these*, because removing this information is very likely to reduce detection 
power.  To generate a report only including SNPs, use vcffilter post-call as 
such:

    freebayes ... | vcffilter -f "TYPE = snp"

### Observation qualities

FreeBayes estimates observation quality using several simple heuristics based 
on manipulations of the phred-scaled base qualities:

* For single-base observations, *mismatches* and *reference observations*: the 
un-adjusted base quality provided in the BAM alignment record.
* For *insertions*: the mean quality of the bases inside of the putatively 
inserted sequence.
* For *deletions*: the mean quality of the bases flanking the putatively 
deleted sequence.
* For *haplotypes*: the mean quality of allele observations within the 
haplotype.

### Effective base depth

By default, filters are left completely open, as both mapping quality and base 
quality are incorporated into the reported site quality (QUAL in the VCF) and 
genotype quality (GQ, when supplying `--genotype-qualities`).  This integration 
is driven by the "Effective Base Depth" metric first developed in 
[snpTools](http://www.hgsc.bcm.edu/software/snptools), which scales observation 
quality by mapping quality.  In short, *P(Obs|Genotype) ~ 
P(MappedCorrectly(Obs))P(SequencedCorrectly(Obs))*.


## Stream processing

FreeBayes can read BAM from standard input `--stdin` instead of directly from
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
[src/bamfiltertech.cpp](http://github.com/ekg/freebayes/blob/master/src/bamfilte
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

1. putative variants lie within `--haplotype-window` bases of each other 
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

    freebayes -f ref.fa --haplotype-window 500 \
        --haplotype-basis-alleles vars.vcf aln.bam >haps.vcf

These steps should allow us to read long haplotypes directly from input data 
with high error rates.

The high error rate means that beyond a small window each read will contain a 
completely different literal haplotype.  To a point, this property improves our 
signal to noise ratio and can effectively filter out sequencing errors at the 
point of the input filters, but it also decreases the effective observation 
depth will prevent the generation of any calls if a long `--haplotype-window` 
is combined with high a sequencing error rate.


## Best practices and design philosophy

FreeBayes follows the patterns suggested by the [Unix 
philosophy](https://en.wikipedia.org/wiki/Unix_philosophy), which promotes the 
development of simple, modular systems that perform a single function, and can 
be combined into more complex systems using stream processing of common 
interchange formats.

FreeBayes incorporates a number of features in order to reduce the complexity 
of variant detection for researchers and developers:

* **Indel realignment is accomplished internally** using a read-independent 
method, and issues resulting from discordant alignments are dramatically 
reducedy through the direct detection of haplotypes.
* The need for **base quality recalibration is avoided** through the direct 
detection of haplotypes. Sequencing platform errors tend to cluster (e.g. at 
the ends of reads), and generate unique, non-repeating haplotypes at a given 
locus.
* **Variant quality recalibration is avoided** by incorporating a number of 
metrics, such as read placement bias and allele balance, directly into the 
Bayesian model.  (Our upcoming publication will discuss this in more detail.)

A minimal pre-processing pipeline similar to that described in "Calling 
variants" should be sufficient for most uses.  For more information, please 
refer to a recent post by Brad Chapman [on minimal BAM preprocessing 
methods](http://bcbio.wordpress.com/2013/10/21/updated-comparison-of-variant-det
ection-methods-ensemble-freebayes-and-minimal-bam-preparation-pipelines/).

For a push-button solution to variant detection, from reads to variant calls, 
look no further than the [gkno genome analysis platform](http://gkno.me/).


## Support

Please report any issues or questions to the [freebayes mailing 
list](https://groups.google.com/forum/#!forum/freebayes), [freebayes issue 
tracker](https://github.com/ekg/freebayes/issues), or by email to 
<erik.garrison@bc.edu>.

Note that if you encounter issues with the development HEAD and you would like 
a quick workaround for an issue that is likely to have been reintroduced 
recently, you can use `git checkout` to step back a few revisions.

    git checkout [git-commit-id]

It will also help with debugging to know if a problem has arisen in recent 
commits!
