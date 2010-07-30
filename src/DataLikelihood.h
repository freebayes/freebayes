#include <iostream>
#include <vector>
#include <utility> // pair
#include <algorithm>
#include <numeric>
#include <vector>
#include <iterator>
#include "Allele.h"
#include "Genotype.h"
#include "Utility.h"
#include "Multinomial.h"


long double
probObservedAllelesGivenGenotypeExact(
        Genotype& genotype,
        vector<Allele*>& observedAlleles);


long double
likelihoodGivenTrueAlleles(
        const vector<Allele*>& observedAlleles,
        const vector<Allele*>& trueAlleles);
