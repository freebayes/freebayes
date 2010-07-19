# calculates data likelihoods for sets of alleles

import multiset
import sys
import cjson
import phred
import json
import math
import operator
from logsumexp import logsumexp

"""
This module attempts to find the best method to approximate the integration of
data likelihoods for the bayesian variant caller we're currently working on.

stdin should be a stream of newline-delimited json records each encoding a list
of alleles which have been parsed out of alignment records.  alleles.cpp in
this distribution provides such a stream.

Erik Garrison <erik.garrison@bc.edu> 2010-07-15
"""

#potential_alleles = [
#        {'type':'reference'}, 
#        {'type':'snp', 'alt':'A'},
#        {'type':'snp', 'alt':'T'},
#        {'type':'snp', 'alt':'G'},
#        {'type':'snp', 'alt':'C'}
#        ]

ploidy = 2
potential_alleles = ['A','T','G','C']

def list_genotypes_to_count_genotypes(genotypes):
    count_genotypes = []
    for genotype in genotypes:
        counts = {}
        for allele in genotype:
            if counts.has_key(allele):
                counts[allele] += 1
            else:
                counts[allele] = 1
        count_genotypes.append(counts.items())
    return count_genotypes

# genotypes are expressed as sets of allele frequencies
genotypes = list_genotypes_to_count_genotypes(list(multiset.multichoose(ploidy, potential_alleles)))


# TODO
# update this so that we aren't just using the 'alternate' field from the alleles
# and are also incorporating the type of allele (ins, deletion, ref, snp)


def group_alleles(alleles):
    groups = {}
    for allele in alleles:
        alt = allele['alt']
        if groups.has_key(alt):
            groups[alt].append(allele)
        else:
            groups[alt] = [allele]
    return groups

def alleles_quality_to_lnprob(alleles):
    for allele in alleles:
        allele['quality'] = phred.phred2ln(allele['quality'])
    return alleles

def product(l):
    return reduce(operator.mul, l)

def observed_alleles_in_genotype(genotype, allele_groups):
    in_genotype = {}
    not_in_genotype = {}
    for key in allele_groups.keys():
        found = False
        for allele, count in genotype:
            if allele == key:
                in_genotype[key] = allele_groups[key]
                found = True
                break
        if not found:
            not_in_genotype[key] = allele_groups[key]
    return in_genotype, not_in_genotype

#def scaled_sampling_prob(genotype, alleles):
#    """The probability of drawing the observations in the allele_groups out of
#    the given genotype, scaled by the number of possible multiset permutations
#    of the genotype (we scale because we don't phase our genotypes under
#    evaluation)."""
#    allele_groups = group_alleles(alleles)
#    if len(allele_groups.items()) == 0:
#        return 0
#    genotype_allele_frequencies = [x[1] for x in genotype]
#    multiplicity = sum(genotype_allele_frequencies)
#    genotype_allele_probabilities = [float(x)/multiplicity for x in genotype_allele_frequencies]
#    observed_allele_frequencies = [len(x) for x in allele_groups.items()]
#    observation_product = 1
#    for allele, count in genotype:
#        if allele_groups.has_key(allele):
#            observation_product *= math.pow(float(count) / multiplicity, len(allele_groups[allele]))
#    return float(math.pow(math.factorial(multiplicity), 2)) \
#        / (product([math.factorial(x) for x in genotype_allele_frequencies]) * 
#                sum([math.factorial(x) for x in observed_allele_frequencies])) \
#        * observation_product
#

def sampling_prob(genotype, alleles):
    """The probability of drawing the observations in alleles out of the given
    genotype, follows the multinomial probability distribution."""
    allele_groups = group_alleles(alleles)
    genotype_allele_frequencies = [x[1] for x in genotype]
    multiplicity = sum(genotype_allele_frequencies)
    observed_allele_frequencies = [len(x) for x in allele_groups.items()]
    observation_product = 1
    for allele, count in genotype:
        if allele_groups.has_key(allele):
            observation_product *= math.pow(float(count) / multiplicity, len(allele_groups[allele]))
    return float(math.factorial(multiplicity)) \
        / sum([math.factorial(x) for x in observed_allele_frequencies]) \
        * observation_product

def likelihood_given_true_alleles(observed_alleles, true_alleles):
    prob = 0
    for o, t in zip(observed_alleles, true_alleles):
        if o['alt'] == t['alt']:
            prob += math.log(1 - math.exp(o['quality']))
        else:
            prob += o['quality']
    return prob

def data_likelihood_exact(genotype, observed_alleles):
    """'Exact' data likelihood, sum of sampling probability * join Q score for
    the observed alleles over all possible underlying 'true allele'
    combinations."""
    observation_count = len(observed_alleles)
    probs = []
    # for all true allele combinations X permutations
    for true_allele_combination in multiset.multichoose(observation_count, [x[0] for x in genotype]):
        for true_allele_permutation in multiset.permutations(true_allele_combination):
            # this mapping allows us to use sampling_prob the same way as we do when we use JSON allele observation records
            true_alleles = [{'alt':allele} for allele in true_allele_permutation]
            sprob = sampling_prob(genotype, true_alleles)
            lnsampling_prob = math.log(sprob)
            probs.append(lnsampling_prob + \
                    likelihood_given_true_alleles(observed_alleles, true_alleles))
    # sum the individual probability of all combinations
    return math.exp(logsumexp(probs))

def data_likelihood_estimate(genotype, alleles_in_genotype, alleles_not_in_genotype, sampling_prob, error_depth=2):
    """Estimates the data likelihood, which is a sum over all possible error
    profiles, or underlying 'true alleles', motivating the observations."""
    # for up to error_depth errors
    pass

def genotype_probabilities(genotypes, alleles):
    return [[str(genotype), data_likelihood_exact(genotype, alleles)] for genotype in genotypes]

if __name__ == '__main__':
    for line in sys.stdin:
        position = cjson.decode(line)
        samples = position['samples']
        for samplename, sample in samples.items():
            # only process snps and reference alleles
            alleles = [allele for allele in sample['alleles'] if allele['type'] in ['reference', 'snp']]
            alleles = alleles_quality_to_lnprob(alleles)
            #groups = group_alleles(alleles)
            sample['alleles'] = alleles
            sample['genotypes'] = genotype_probabilities(genotypes, alleles)
        print json.dumps(position, indent=2)
