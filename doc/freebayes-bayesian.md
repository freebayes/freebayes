# FreeBayes Bayesian Variant Calling

## Overview

FreeBayes uses Bayesian inference to calculate the probability that different genotypes are present at each position in the genome, given the observed sequencing reads.

## The Core Bayesian Framework

FreeBayes applies Bayes' theorem to compute the posterior probability of each possible genotype:

```
P(genotype | data) ∝ P(data | genotype) × P(genotype)
```

Where:
- **P(genotype | data)** is the posterior probability - what we want to know
- **P(data | genotype)** is the likelihood - how probable the observed reads are given a particular genotype
- **P(genotype)** is the prior probability - what we expect before seeing the data

## The Likelihood Component

For the likelihood, FreeBayes models how likely the observed reads are under each possible genotype hypothesis. It considers factors like:

- Base quality scores (probability of sequencing error)
- Mapping quality (confidence the read is correctly aligned)
- Allele balance (ratio of reads supporting different alleles)
- The probability of observing each read given the diploid genotype at that position

## The Prior Component

FreeBayes uses priors based on:

- Expected heterozygosity in the population (theta parameter, typically ~0.001)
- Ploidy model (diploid, haploid, or pooled samples)
- Population genetic expectations (e.g., Hardy-Weinberg equilibrium)

## Why Bayesian?

This approach naturally handles uncertainty by producing probability distributions rather than binary yes/no calls. FreeBayes can integrate information across multiple samples simultaneously and properly account for sequencing errors, mapping uncertainty, and biological variation in a principled probabilistic framework.

## Example: Computing the Prior for Diploid Genomes

### Setup

Consider a position where the reference allele is **A**, and we've observed some reads with a **G**. For a diploid individual, the possible genotypes are:

- **A/A** (homozygous reference)
- **A/G** (heterozygous)
- **G/G** (homozygous alternate)

### Computing the Prior

FreeBayes uses the population mutation rate parameter θ (theta), which represents the expected heterozygosity. A typical value is θ = 0.001.

From population genetics theory (Wright-Fisher model for diploid populations), the prior probabilities are approximately:

**P(A/A)** = 1 - θ + θ²/3 ≈ 1 - θ ≈ 0.999

**P(A/G)** = θ - θ²/3 ≈ θ ≈ 0.001

**P(G/G)** = θ²/3 ≈ 0.000000333

### The Intuition

With θ = 0.001:

- There's a ~99.9% prior probability the site matches the reference (A/A)
- There's a ~0.1% prior probability of heterozygosity (A/G)
- There's a negligible prior probability of a homozygous variant (G/G)

This makes sense for diploid genomes: most positions in any individual genome match the reference, heterozygous sites are uncommon, and homozygous alternate genotypes at previously unseen sites are very rare.

### Why These Specific Formulas?

For a diploid individual, if the population mutation rate is θ:

- Most of the time, both alleles match the reference: P(A/A) ≈ 1 - θ
- Sometimes you inherit one mutant allele: P(A/G) ≈ θ
- Rarely do you inherit two independent mutant alleles at the same position: P(G/G) ≈ θ²

The θ²/3 terms are corrections from the exact Wright-Fisher model, but they're small enough that the approximations work well.

## Multiple Samples

When calling variants on multiple samples jointly, FreeBayes considers more complex priors that account for allele frequencies across the population, making variants more likely if they're seen in multiple individuals.

## Key Insight

These priors get updated by the likelihood from the actual read data, so strong evidence (many high-quality reads supporting a variant) can overcome even a low prior probability.

## Haploid vs Diploid

### Diploid (e.g., autosomes)
- Two alleles per position
- Three possible genotypes at a biallelic site: A/A, A/G, G/G

### Haploid (e.g., mitochondrial DNA, Y chromosome, bacteria)
- One allele per position
- Two possible genotypes: A or G
- Simpler priors: P(A) ≈ 1 - θ and P(G) ≈ θ

## Local realignment

FreeBayes performs realignment "on the fly" for each candidate variant site:

- Identify candidate haplotypes: At each position, FreeBayes examines the reads and generates candidate haplotypes (possible sequences) that could explain the observed data, including different combinations of SNPs, indels, and complex variants.
- Realign reads to haplotypes: For each read overlapping the region, FreeBayes realigns it against each candidate haplotype using a pair-HMM (Hidden Markov Model) or similar alignment algorithm. This produces an alignment likelihood for each read-haplotype pair.
- Compute likelihoods: The likelihood P(data | genotype) in the Bayesian framework is computed based on these haplotype-aware alignments, not the original BAM alignments.
- Call variants: The Bayesian genotype calling uses these realignment-based likelihoods.

This integrated approach means:

- You don't need to run a separate realignment preprocessing step (unlike older GATK workflows that required IndelRealigner)
- Reads are realigned specifically in the context of the variants being evaluated
- It naturally handles complex variants where multiple nearby SNPs and indels interact
- The uncertainty from realignment is properly propagated into the variant calls

How it works:

- Candidate Generation Phase: FreeBayes uses heuristics and initial observations from the reads to generate a set of candidate haplotypes (possible variant combinations). This is influenced by:
- Read evidence (what variants appear in reads)
- Prior expectations (minimum allele frequency thresholds)
- Complexity limits (how many variants to consider jointly)
- Realignment Phase: Reads are realigned against each candidate haplotype independently. This produces likelihoods P(read | haplotype).
- Bayesian Calling Phase: The posterior probabilities P(genotype | data) are computed using these likelihoods plus priors.

The Bayesian calling therefore doesn't iterate back to refine the realignments themselves. However, the set of candidate haplotypes that reads are realigned against is influenced by:

- Evidence strength in the reads
- Prior probability considerations (FreeBayes won't generate wildly implausible haplotypes)
- Computational tractability (limiting the haplotype space)
