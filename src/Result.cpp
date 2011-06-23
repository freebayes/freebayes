#include "Result.h"

void Result::sortDataLikelihoods(void) {
    SampleDataLikelihoodCompare datalikelihoodCompare;
    sort(begin(), end(), datalikelihoodCompare);
}
