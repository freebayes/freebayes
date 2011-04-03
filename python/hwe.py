from dirichlet import multinomial, multinomialln, multinomial_coefficient, multinomial_coefficientln
import math
import operator

# @genotype is counts of A,B,C etc. alleles, e.g. (2,0) is AA and (1,1) is AB
# @allele_counts is the counts of the alleles in the population
# @population_total_alleles is the number of alleles at the locus
def hwe_expectation(genotype, allele_counts, population_total_alleles):
    # how many counts of this genotype do we expect under equilibrium?
    # get the multinomial coefficient
    genotype_coeff = multinomial_coefficient(len(genotype), genotype)
    allele_frequencies = [count / float(population_total_alleles) for count in allele_counts]
    genotype_expected_frequency = genotype_coeff * reduce(operator.mul, [math.pow(freq, p) for freq, p in zip(allele_frequencies, genotype)])
    return genotype_expected_frequency

