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
