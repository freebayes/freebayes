#include "DataLikelihood.h"
#include "multichoose.h"
#include "multipermute.h"



// log probability of a given matching of true and observed alleles
// true alleles are the 'real' alleles
// observed alleles are whet we observe
// a mismatch between these two sets implies an error, which is scored accordingly here
long double likelihoodGivenTrueAlleles(vector<Allele*>& observedAlleles, vector<Allele*>& trueAlleles) {
    long double prob = 0;
    vector<Allele*>::iterator o = observedAlleles.begin();
    vector<Allele*>::iterator t = trueAlleles.begin();
    for ( ; o != observedAlleles.end() && t != trueAlleles.end(); ++o, ++t)
    {
        Allele* observedAllele = *o;
        Allele* trueAllele = *t;
        if (observedAllele == trueAllele) {
            prob += log(1 - exp(observedAllele->lncurrentQuality()));
        } else {
            prob += observedAllele->lncurrentQuality();
        }
    }
    return prob;
}

/*
  'Exact' data likelihood, sum of sampling probability * joint Q score for the
  observed alleles over all possible underlying 'true allele' combinations."""
*/
long double
probObservedAllelesGivenGenotype(
        vector<Allele*>& observedAlleles,
        Genotype& genotype) {

    int observationCount = observedAlleles.size();
    vector<long double> alleleProbs = genotype.alleleProbabilities();
    vector<long double> probs;
    vector<vector<Allele*> > combos = multichoose_ptr(observationCount, genotype.alleles);
    for (vector<vector<Allele*> >::iterator combo = combos.begin(); combo != combos.end(); ++combo) {
        vector<vector<Allele*> > trueAllelePermutations = multipermute(*combo);
        for (vector<vector<Allele*> >::iterator perm = trueAllelePermutations.begin(); perm != trueAllelePermutations.end(); ++perm) {
            vector<Allele*>& trueAlleles = *perm;
            map<Allele, int> trueAlleleCounts = countAlleles(trueAlleles);
            vector<int> observationCounts; // counts of 'observations' of true alleles, ordered according to the genotype's internal ordering
            for (Genotype::const_iterator g = genotype.begin(); g != genotype.end(); ++g) {
                map<Allele, int>::const_iterator count = trueAlleleCounts.find(g->second);
                if (count != trueAlleleCounts.end()) {
                    observationCounts.push_back(count->second);
                } else {
                    observationCounts.push_back(0);
                }
            }
            probs.push_back(multinomialln(alleleProbs, observationCounts) + likelihoodGivenTrueAlleles(observedAlleles, trueAlleles));
        }
    }
    return logsumexp(probs);
}

vector<pair<Genotype, long double> >
probObservedAllelesGivenGenotypes(
        vector<Allele*>& observedAlleles,
        vector<Genotype>& genotypes) {
    vector<pair<Genotype, long double> > results;
    for (vector<Genotype>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
        results.push_back(make_pair(*g, probObservedAllelesGivenGenotype(observedAlleles, *g)));
    }
    return results;
}
