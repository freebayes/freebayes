#include "GenotypePriors.h"

/*
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
*/

map<Allele, int> countAlleles(vector<Genotype*>& genotypeCombo) {
    map<Allele, int> alleleCounts;
    for (vector<Genotype*>::iterator g = genotypeCombo.begin(); g != genotypeCombo.end(); ++g) {
        for (Genotype::iterator a = (*g)->begin(); a != (*g)->end(); ++a) {
            map<Allele, int>::iterator c = alleleCounts.find(a->allele);
            if (c != alleleCounts.end()) {
                c->second += a->count;
            } else {
                alleleCounts.insert(make_pair(a->allele, a->count));
            }
        }
    }
    return alleleCounts;
}

map<int, int> countFrequencies(vector<Genotype*>& genotypeCombo) {
    map<int, int> frequencyCounts;
    map<Allele, int> alleles = countAlleles(genotypeCombo);
    for (map<Allele, int>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
        map<int, int>::iterator c = frequencyCounts.find(a->second);
        if (c != frequencyCounts.end()) {
            c->second += 1;
        } else {
            frequencyCounts[a->second] = 1;
        }
    }
    return frequencyCounts;
}

/*
def allele_frequency_probability(allele_frequency_counts, theta=0.001):
    """Implements Ewens' Sampling Formula.  allele_frequency_counts is a
    dictionary mapping count -> number of alleles with this count in the
    population."""
    M = sum([frequency * count for frequency, count in allele_frequency_counts.iteritems()])
    return math.factorial(M) \
        / (theta * product([theta + h for h in range(1, M)])) \
        * product([math.pow(theta, count) / math.pow(frequency, count) * math.factorial(count) \
            for frequency, count in allele_frequency_counts.iteritems()])
*/

long double alleleFrequencyProbability(map<int, int> alleleFrequencyCounts, long double theta) {

    int M = 0;
    long double p = 1;

    for (map<int, int>::iterator f = alleleFrequencyCounts.begin(); f != alleleFrequencyCounts.end(); ++f) {
        int frequency = f->first;
        int count = f->second;
        M += frequency * count;
        p *= pow(theta, count) / pow(frequency, count) * factorial(count);
    }

    long double thetaH = 1;
    for (int h = 1; h < M; ++h)
        thetaH *= theta + h;

    return factorial(M) / (theta * thetaH) * p;

}

/*
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
*/

long double alleleFrequencyProbabilityln(map<int, int> alleleFrequencyCounts, long double theta) {

    int M = 0; // multiplicity of site
    long double p = 0;
    long double thetaln = log(theta);

    for (map<int, int>::iterator f = alleleFrequencyCounts.begin(); f != alleleFrequencyCounts.end(); ++f) {
        int frequency = f->first;
        int count = f->second;
        M += frequency * count;
        p += powln(thetaln, count) - powln(log(frequency), count) + factorialln(count);
    }

    long double thetaH = 0;
    for (int h = 1; h < M; ++h)
        thetaH += log(theta + h);

    return factorialln(M) - (thetaln + thetaH) + p;

}

int alleleFrequency(GenotypeCombo& genotypeCombo, Allele& allele) {

    int frequency = 0;
    for (GenotypeCombo::iterator gc = genotypeCombo.begin(); gc != genotypeCombo.end(); ++gc) {
        frequency += gc->genotype->alleleFrequency(allele);
    }
    return frequency;

}

// TODO XXX only works for the diploid case!!!
long double probabilityDiploidGenotypeComboGivenAlleleFrequencyln(GenotypeCombo& genotypeCombo, Allele& allele) {

    int n = genotypeCombo.size();
    int h = 0; // number of heterozygotes
    int f = alleleFrequency(genotypeCombo, allele);

    for (GenotypeCombo::iterator gc = genotypeCombo.begin(); gc != genotypeCombo.end(); ++gc) {
        if (!gc->genotype->homozygous())
            ++h;
    }

    return powln(log(2), h) - (factorialln(2 * n) - (factorialln(f) + factorialln(2 * n - f)));

}


// XXX although this is intended to replace the above diploid-specific
// function, and will do so with identical results on diploid data, it is
// unproven theoretically for non-uniform copy number
long double probabilityGenotypeComboGivenAlleleFrequencyln(GenotypeCombo& genotypeCombo, Allele& allele) {

    int n = genotypeCombo.numberOfAlleles();
    int h = 0; // number of heterozygotes
    int f = alleleFrequency(genotypeCombo, allele);
    long double lnploidyscalar = 0;

    for (GenotypeCombo::iterator gc = genotypeCombo.begin(); gc != genotypeCombo.end(); ++gc) {
        if (!gc->genotype->homozygous()) {
            lnploidyscalar += log(gc->genotype->ploidy);
            ++h;
        }
    }

    return lnploidyscalar - (factorialln(n) - (factorialln(f) + factorialln(n - f)));

}
