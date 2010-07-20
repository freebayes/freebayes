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
            prob = lnsampling_prob + likelihood_given_true_alleles(observed_alleles, true_alleles)
            #print genotype, math.exp(prob), sprob, true_allele_permutation, [o['alt'] for o in observed_alleles]
            probs.append(prob)
    # sum the individual probability of all combinations
    return logsumexp(probs)

def data_likelihood_estimate(genotype, alleles):
    """Estimates the data likelihood, which is a sum over all possible error
    profiles, or underlying 'true alleles', motivating the observations."""
    # for up to error_depth errors
    pass

def genotype_combination_sampling_probability(genotype_combination, observed_alleles):
    multiplicity = math.log(ploidy * len(genotype_combination))
    result = 1 - multiplicity
    allele_groups = group_alleles(observed_alleles)
    for allele, observations in allele_groups.iteritems():
        result += math.log(math.factorial(len(observations)))
    # scale by product of multiset permutations of all genotypes in combo
    for combo in genotype_combination:
        for genotype in combo:
            m_i = sum([a[1] for a in genotype])
            result += math.log(math.factorial(m_i))
            result -= sum([math.log(math.factorial(allele[1])) for allele in genotype])
    return result

def count_frequencies(genotype_combo):
    counts = {}
    alleles = {}
    for genotype in genotype_combo:
        for allele, count in genotype:
            if alleles.has_key(allele):
                alleles[allele] += count
            else:
                alleles[allele] = count
    for allele, count in alleles.iteritems():
        if counts.has_key(count):
            counts[count] += 1
        else:
            counts[count] = 1
    return counts

def allele_frequency_probability(allele_frequency_counts, theta=0.01):
    """Implements Ewens' Sampling Formula.  allele_frequency_counts is a
    dictionary mapping count -> number of alleles with this count in the
    population."""
    M = sum([frequency * count for frequency, count in allele_frequency_counts.iteritems()])
    return math.factorial(M) \
        / (theta * product([theta + h for h in range(1, M)])) \
        * product([math.pow(theta, count) / math.pow(frequency, count) * math.factorial(count) \
            for frequency, count in allele_frequency_counts.iteritems()])

def genotype_probabilities(genotypes, alleles):
    return [[str(genotype), data_likelihood_exact(genotype, alleles)] for genotype in genotypes]

def genotype_probabilities_heuristic(genotypes, alleles):
    groups = group_alleles(alleles)
    # group genotypes relative to the groups of observed alleles
    # take the first member of each group and apply our data likelihood calculation
    # then apply it to the rest
    if len(groups.keys()) is 1:
        # we can cleanly do all-right, part-right, all-wrong
        pass
    if len(groups.keys()) is 2:
        # we can do all-right, two types of 'part-right', and all-wrong
        pass

def multiset_banded_genotype_combinations(sample_genotypes, bandwidth):
    for index_combo in multiset.multichoose(len(samples), range(bandwidth)):
        for index_permutation in multiset.permutations(index_combo):
            yield [genotypes[index] for index, genotypes in zip(index_combo, sample_genotypes)]

# TODO you should implement gabor's banding solution; the above multiset method
# is comically large and produces incorrect results despite the computational load
def banded_genotype_combinations(sample_genotypes, bandwidth):
    for i in range(bandwidth):
        pass

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
            sample['genotypes'] = [[genotype, data_likelihood_exact(genotype, alleles)] for genotype in genotypes]
            #sample['genotypes_estimate'] = [[str(genotype), data_likelihood_estimate(genotype, alleles)] for genotype in genotypes]
        # estimate the posterior over all genotype combinations within some indexed bandwidth of optimal
        sample_genotypes = [sorted(sample['genotypes'], key=lambda genotype: genotype[1], reverse=True) for sample in samples.values()]
        genotype_combo_probs = []
        #for combo in multiset_banded_genotype_combinations(sample_genotypes, 2):
        for combo in banded_genotype_combinations(sample_genotypes, 2):
            #print combo
            probability_observations_given_genotypes = logsumexp([prob for genotype, prob in combo])
            prior_probability_of_genotype = allele_frequency_probability(count_frequencies([genotype for genotype, prob in combo]))
            genotype_combo_probs.append([combo, math.log(prior_probability_of_genotype) + probability_observations_given_genotypes])
        # sum, use to normalize
        # apply bayes rule
        posterior_normalizer = logsumexp([prob for combo, prob in genotype_combo_probs])
        # XXX incorrect... this depends on the priors
        #best_genotype_combo = [s[0] for s in sample_genotypes]
        best_genotype_combo_prob = None
        best_genotype_combo = []
        for combo, prob in genotype_combo_probs:
            if not best_genotype_combo_prob:
                best_genotype_combo_prob = prob
            if prob <= best_genotype_combo_prob:
                best_genotype_combo_prob = prob
                best_genotype_combo = combo
        best_genotype_probability = math.exp(logsumexp([prob for genotype, prob in best_genotype_combo]) - posterior_normalizer)
        position['genotyping'] = {'prob':best_genotype_probability, 'best_combo':[[str(genotype), math.exp(prob)] for genotype, prob in best_genotype_combo]}
        # TODO estimate marginal probabilities of genotypings
        # here we cast everything into float-space
        for samplename, sample in samples.items():
            sample['genotypes'] = [[str(genotype), math.exp(prob)] for genotype, prob in sample['genotypes']]
        print cjson.encode(position)
