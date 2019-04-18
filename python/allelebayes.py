#!/usr/bin/env python

# calculates data likelihoods for sets of alleles

import multiset
import sys
import cjson
import phred
import json
import math
import operator
from logsumexp import logsumexp
from dirichlet import dirichlet_maximum_likelihood_ratio, dirichlet, multinomial, multinomialln
from factorialln import factorialln

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

"""
ploidy = 2
potential_alleles = ['A','T','G','C']

# genotypes are expressed as sets of allele frequencies
genotypes = list_genotypes_to_count_genotypes(list(multiset.multichoose(ploidy, potential_alleles)))
"""


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


# TODO XXX
# Yes, this is the sampling probability.  It is the multinomial sampling
# probability, which is the specific probability of a specific set of
# categorical outcomes.  Unfortunately, this is not what we really want here.
# What we want is the prior probability that a given set of draws come out of a
# given multiset (genotype, in our case).  I believe that this is given by the
# Dirichlet distribution.  Investigate.
def sampling_prob(genotype, alleles):
    """The specific probability of drawing the observations in alleles out of the given
    genotype, follows the multinomial probability distribution."""
    allele_groups = group_alleles(alleles)
    multiplicity = sum([x[1] for x in genotype])
    print genotype, multiplicity, alleles
    for allele, count in genotype:
        if allele_groups.has_key(allele):
            print allele, count, math.pow(float(count) / multiplicity, len(allele_groups[allele]))
    print product([math.factorial(len(obs)) for obs in allele_groups.values()])
    print allele_groups.values()
    return float(math.factorial(len(alleles))) \
        / product([math.factorial(len(obs)) for obs in allele_groups.values()]) \
        * product([math.pow(float(count) / multiplicity, len(allele_groups[allele])) \
                    for allele, count in genotype if allele_groups.has_key(allele)])

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
    #print "probability that observations", [o['alt'] for o in observed_alleles], "arise from genotype", genotype
    observation_count = len(observed_alleles)
    ploidy = sum([count for allele, count in genotype])
    allele_probs = [count / float(ploidy) for allele, count in genotype]
    probs = []
    # for all true allele combinations X permutations
    for true_allele_combination in multiset.multichoose(observation_count, [x[0] for x in genotype]):
        for true_allele_permutation in multiset.permutations(true_allele_combination):
            # this mapping allows us to use sampling_prob the same way as we do when we use JSON allele observation records
            true_alleles = [{'alt':allele} for allele in true_allele_permutation]
            allele_groups = group_alleles(true_alleles)
            observations = []
            for allele, count in genotype:
                if allele_groups.has_key(allele):
                    observations.append(len(allele_groups[allele]))
                else:
                    observations.append(0)
            #sprob = dirichlet_maximum_likelihood_ratio(allele_probs, observations) # distribution parameter here
            lnsampling_prob = multinomialln(allele_probs, observations)
            prob = lnsampling_prob + likelihood_given_true_alleles(observed_alleles, true_alleles)
            #print math.exp(prob), sprob, genotype, true_allele_permutation
            #print genotype, math.exp(prob), sprob, true_allele_permutation, [o['alt'] for o in observed_alleles]
            probs.append(prob)
    # sum the individual probability of all combinations
    p = logsumexp(probs)
    #print math.exp(p)
    return p

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

def allele_frequency_probability(allele_frequency_counts, theta=0.001):
    """Implements Ewens' Sampling Formula.  allele_frequency_counts is a
    dictionary mapping count -> number of alleles with this count in the
    population."""
    M = sum([frequency * count for frequency, count in allele_frequency_counts.iteritems()])
    return math.factorial(M) \
        / (theta * product([theta + h for h in range(1, M)])) \
        * product([math.pow(theta, count) / math.pow(frequency, count) * math.factorial(count) \
            for frequency, count in allele_frequency_counts.iteritems()])

def powln(n, m):
    """Power of number in log space"""
    return sum([n] * m)

def allele_frequency_probabilityln(allele_frequency_counts, theta=0.001):
    """Log space version to avoid inevitable overflows with coverage >100.
    Implements Ewens' Sampling Formula.  allele_frequency_counts is a
    dictionary mapping count -> number of alleles with this count in the
    population."""
    thetaln = math.log(theta)
    M = sum([frequency * count for frequency, count in allele_frequency_counts.iteritems()])
    return factorialln(M) \
        - (thetaln + sum([math.log(theta + h) for h in range(1, M)])) \
        + sum([powln(thetaln, count) - powln(math.log(frequency), count) + factorialln(count) \
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
            yield [genotypes[index] for index, genotypes in zip(index_permutation, sample_genotypes)]

# TODO you should implement gabor's banding solution; the above multiset method
# is comically large and produces incorrect results despite the computational load
def banded_genotype_combinations(sample_genotypes, bandwidth, band_depth):
    # always provide the 'best' case
    yield [(sample, genotypes[0]) for sample, genotypes in sample_genotypes]
    for i in range(1, bandwidth):
        for j in range(1, band_depth):  # band_depth is the depth to which we explore the bandwith... TODO explain better
            indexes = j * [i] + (len(sample_genotypes) - j) * [0]
            for index_permutation in multiset.permutations(indexes):
                yield [(sample, genotypes[index]) for index, (sample, genotypes) in zip(index_permutation, sample_genotypes)]

def genotype_str(genotype):
    return reduce(operator.add, [allele * count for allele, count in genotype])

if __name__ == '__main__':

    ploidy = 2 # assume ploidy 2 for all individuals and all positions

    potential_alleles = ['A','T','G','C']

    # genotypes are expressed as sets of allele frequencies
    genotypes = list_genotypes_to_count_genotypes(list(multiset.multichoose(ploidy, potential_alleles)))

    for line in sys.stdin:
        position = cjson.decode(line)
        #print position['position']
        samples = position['samples']

        position['coverage'] = sum([len(sample['alleles']) for samplename, sample in samples.iteritems()])

        #potential_alleles = ['A','T','G','C']
        potential_alleles = set()
        for samplename, sample in samples.items():
            # only process snps and reference alleles
            alleles = [allele for allele in sample['alleles'] if allele['type'] in ['reference', 'snp']]
            alleles = alleles_quality_to_lnprob(alleles)
            sample['alleles'] = alleles
            potential_alleles = potential_alleles.union(set([allele['alt'] for allele in alleles]))

        position['filtered coverage'] = sum([len(sample['alleles']) for samplename, sample in samples.iteritems()])

        # genotypes are expressed as sets of allele frequencies
        #genotypes = list_genotypes_to_count_genotypes(list(multiset.multichoose(ploidy, list(potential_alleles))))

        for samplename, sample in samples.items():
            alleles = sample['alleles']
            groups = group_alleles(alleles)
            sample['genotypes'] = [[genotype, data_likelihood_exact(genotype, alleles)] for genotype in genotypes]
            #sample['genotypes_estimate'] = [[str(genotype), data_likelihood_estimate(genotype, alleles)] for genotype in genotypes]
        # estimate the posterior over all genotype combinations within some indexed bandwidth of optimal
        # TODO preserve sample names in the genotype comos
        sample_genotypes = [(name, sorted(sample['genotypes'], key=lambda genotype: genotype[1], reverse=True)) for name, sample in samples.iteritems()]
        genotype_combo_probs = []
        #for combo in multiset_banded_genotype_combinations(sample_genotypes, 2):
        #for combo in banded_genotype_combinations(sample_genotypes, min(len(genotypes), 2), len(samples)):
        # now marginals time...
        marginals = {}
        for name, sample in samples.iteritems():
            marginals[name] = {}

        combos_tested = 0
        for combo in banded_genotype_combinations(sample_genotypes, min(len(genotypes), 2), 2):
            combos_tested += 1
            probability_observations_given_genotypes = sum([prob for name, (genotype, prob) in combo])
            frequency_counts = count_frequencies([genotype for name, (genotype, prob) in combo])
            prior_probability_of_genotype = allele_frequency_probabilityln(frequency_counts)
            combo_prob = prior_probability_of_genotype + probability_observations_given_genotypes
            for name, (genotype, prob) in combo:
                gstr = genotype_str(genotype)
                if marginals[name].has_key(gstr):
                    marginals[name][gstr].append(combo_prob)
                else:
                    marginals[name][gstr] = [combo_prob]
            genotype_combo_probs.append([combo, combo_prob])

        genotype_combo_probs = sorted(genotype_combo_probs, key=lambda c: c[1], reverse=True)
        #for line in [json.dumps({'prob':prior_probability_of_genotype, 'combo':combo}) for combo, prior_probability_of_genotype in genotype_combo_probs]:
        #    print line

        # sum, use to normalize
        # apply bayes rule

        #print genotype_combo_probs
        #print [prob for combo, prob in genotype_combo_probs]
        #for combo, prob in genotype_combo_probs:
        #    print prob
        posterior_normalizer = logsumexp([prob for combo, prob in genotype_combo_probs])

        # handle marginals
        for sample, genotype_probs in marginals.iteritems():
            for genotype, probs in genotype_probs.iteritems():
                marginals[sample][genotype] = logsumexp(probs) - posterior_normalizer

        best_genotype_combo = genotype_combo_probs[0][0]
        best_genotype_combo_prob = genotype_combo_probs[0][1]

        #best_genotype_probability = math.exp(sum([prob for name, (genotype, prob) in best_genotype_combo]) \
        #        + allele_frequency_probabilityln(count_frequencies([genotype for name, (genotype, prob) in best_genotype_combo])) \
        #        - posterior_normalizer)
        best_genotype_probability = math.exp(best_genotype_combo_prob - posterior_normalizer)
        position['best_genotype_combo'] = [[name, genotype_str(genotype), math.exp(marginals[name][genotype_str(genotype)])] 
                                                  for name, (genotype, prob) in best_genotype_combo]
        position['best_genotype_combo_prob'] = best_genotype_probability
        position['posterior_normalizer'] = math.exp(posterior_normalizer)
        position['combos_tested'] = combos_tested
        #position['genotype_combo_probs'] = genotype_combo_probs
        # TODO estimate marginal probabilities of genotypings
        # here we cast everything into float-space
        for samplename, sample in samples.items():
            sample['genotypes'] = sorted([[genotype_str(genotype), math.exp(prob)] for genotype, prob in sample['genotypes']], 
                                            key=lambda c: c[1], reverse=True)

        print cjson.encode(position)
        #print position['position']
