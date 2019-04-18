from dirichlet import multinomial, multinomialln, multinomial_coefficient, multinomial_coefficientln
import math
import operator

def hwe_expectation(genotype, allele_counts):
    """@genotype is counts of A,B,C etc. alleles, e.g. (2,0) is AA and (1,1) is AB
    @allele_counts is the counts of the alleles in the population"""
    population_total_alleles = sum(allele_counts)
    ploidy = sum(genotype)
    genotype_coeff = multinomial_coefficient(ploidy, genotype)
    allele_frequencies = [count / float(population_total_alleles) for count in allele_counts]
    genotype_expected_frequency = genotype_coeff * reduce(operator.mul, [math.pow(freq, p) for freq, p in zip(allele_frequencies, genotype)])
    return genotype_expected_frequency


def add_tuple(a, b):
    l = []
    for i,j in zip(a,b):
        l.append(i + j)
    return tuple(l)

# TODO handle NULL case, genotype has frequency 0
def hwe_sampling_probln(genotype, genotypes, ploidy):
    """@genotype:  counts of A,B,C etc. alleles, e.g. (2,0) is AA and (1,1) is AB
    @genotypes: counts of the genotypes in the population,
                     e.g. (((2,0),1), ((1,1),2), ((0,2),1)) would be 1xAA, 2xAB, 1xBB
    @return: the probability that there are exactly as many "genotype" in the
             population as suggested by the genotype counts, given HWE"""
    population_total_alleles = sum([sum(g[0]) * g[1] for g in genotypes])
    #print "population_total_alleles", population_total_alleles
    allele_counts = reduce(add_tuple, [[a * g[1] for a in g[0]] for g in genotypes])
    #print "allele_counts", allele_counts
    genotype_counts = [g[1] for g in genotypes]
    #print "genotype_counts", genotype_counts
    population_total_genotypes = sum([gtc[1] for gtc in genotypes])
    #print "total genotypes:", population_total_genotypes
    #print "ploidy", ploidy
    # number of arrangements of the alleles in the sample
    arrangements_of_alleles_in_sample = multinomial_coefficientln(population_total_alleles, allele_counts)
    #print "arrangements_of_alleles_in_sample", math.exp(arrangements_of_alleles_in_sample)
    # number of arrangements which contain exactly count genotypes
    #arrangements_with_exactly_count_genotype = multinomial_coefficientln(population_total_genotypes, genotype_counts)
    #print math.exp(multinomial_coefficientln(ploidy, genotype))
    arrangements_with_exactly_count_genotype = \
            multinomial_coefficientln(ploidy, genotype) \
            + multinomial_coefficientln(population_total_genotypes, genotype_counts)
    #print "arrangements_with_exactly_count_genotype", math.exp(arrangements_with_exactly_count_genotype)
    return arrangements_with_exactly_count_genotype - arrangements_of_alleles_in_sample;

def inbreeding_coefficient(genotype, genotypes):
    population_total_alleles = sum([sum(g[0]) * g[1] for g in genotypes])
    allele_counts = reduce(add_tuple, [[a * g[1] for a in g[0]] for g in genotypes])
    genotype_counts = [g[1] for g in genotypes]
    population_total_genotypes = sum([gtc[1] for gtc in genotypes])
    expected = hwe_expectation(genotype, allele_counts) * population_total_genotypes
    observed = 0
    for g in genotypes:
        if g[0] == genotype:
            observed = g[1]
            break
    if observed == 0:
        print "error, no observations of genotype, cannot calculate inbreeding coefficient"
        return 0
    return 1 - (observed / expected)
