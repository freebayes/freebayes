#include "ResultData.h"

using namespace std;


void json(ostream& out, Results& results, AlleleParser* parser) {
    out << "{";
    for (Results::iterator r = results.begin(); r != results.end(); ++r) {
        ResultData& sample = r->second;
        if (r != results.begin()) out << ",";
        out << "\"" << sample.name << "\":{"
            << "\"coverage\":" << sample.observations->observationCount() << ","
            << "\"genotypes\":[";
        for (map<Genotype*, long double>::iterator g = sample.marginals.begin(); 
                g != sample.marginals.end(); ++g) {
            if (g != sample.marginals.begin()) cout << ",";
            out << "[\"" << *(g->first) << "\"," << safe_exp(g->second) << "]";
        }
        out << "]";
        if (parser->parameters.outputAlleles)
            out << ",\"alleles\":" << sample.observations->json();
        out << "}";

    }
    out << "}";
}

// current date string in YYYYMMDD format
string dateStr(void) {

    time_t rawtime;
    struct tm* timeinfo;
    char buffer[80];

    time(&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer, 80, "%Y%m%d", timeinfo);

    return string(buffer);

}

void vcfHeader(ostream& out,
        string referenceName,
        vector<string>& samples,
        Parameters& parameters) {

    out << "##fileformat=VCFv4.0" << endl
        << "##fileDate=" << dateStr() << endl
        << "##source=freeBayes version " << FREEBAYES_VERSION << endl
        << "##reference=" << referenceName << endl
        << "##phasing=none" << endl
        << "##commandline=\"" << parameters.commandline << "\"" << endl
        << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">" << endl
        << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth at the locus\">" << endl
        << "##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Total number of alternate alleles in called genotypes\">" << endl
        << "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">" << endl
        //<< "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Estimated allele frequency in the range (0,1]\">" << endl
        << "##INFO=<ID=BCF,Number=2,Type=Integer,Description=\"Forward-strand base count: the number of observations on the forward strand for, the reference, the first alternate, the second alternate, and so on for each alternate\">" << endl
        << "##INFO=<ID=BCR,Number=2,Type=Integer,Description=\"Reverse-strand base count: the number of observations on the reverse strand for, the reference, the first alternate, the second alternate, and so on for each alternate\">" << endl
        << "##INFO=<ID=SB,Number=1,Type=Float,Description=\"Strand bias for the alternate allele: a number between 0 and 1 representing the ratio of forward strand sequence reads showing the alternate allele to all reads, considering only reads from individuals called as heterozygous\">" << endl
        << "##INFO=<ID=AB,Number=1,Type=Integer,Description=\"Allele balance at heterozygous sites: a number beween 0 and 1 representing the ratio of reads showing the reference allele to all reads, considering only reads from individuals called as heterozygous\">" << endl
        << "##INFO=<ID=ABB,Number=2,Type=Integer,Description=\"Allele balance counts: two numbers giving the numbers of sequence reads from apparent heterozygotes which show reference and alternate alleles for the site\">" << endl
        << "##INFO=<ID=RUN,Number=1,Type=Integer,Description=\"Homopolymer run length: the number of consecutive nucleotides in the reference genome matching the alternate allele prior to the current position\">" << endl
        << "##INFO=<ID=SNP,Number=0,Type=Flag,Description=\"SNP allele\"" << endl
        << "##INFO=<ID=TS,Number=0,Type=Flag,Description=\"transition SNP\"" << endl
        << "##INFO=<ID=TV,Number=0,Type=Flag,Description=\"transversion SNP\"" << endl
        << "##INFO=<ID=CpG,Number=0,Type=Flag,Description=\"CpG site (either CpG, TpG or CpA)\"" << endl
        << "##INFO=<ID=MNP,Number=0,Type=Integer,Description=\"Lengeth of MNP allele, if present\"" << endl
        << "##INFO=<ID=INS,Number=1,Type=Integer,Description=\"Length of insertion allele, if present\"" << endl
        << "##INFO=<ID=DEL,Number=1,Type=Integer,Description=\"Length of deletion allele, if present\"" << endl
        << "##FORMAT=<ID=GT,Number=1,Type=String,\"Genotype\">" << endl
        << "##FORMAT=<ID=GQ,Number=1,Type=Integer,\"Genotype Quality\">" << endl
        << "##FORMAT=<ID=DP,Number=1,Type=Integer,\"Read Depth\">" << endl
        << "##FORMAT=<ID=RA,Number=1,Type=Integer,\"Reference allele observations\">" << endl
        << "##FORMAT=<ID=AA,Number=1,Type=Integer,\"Alternate allele observations\">" << endl
        << "##FORMAT=<ID=BCF,Number=2,Type=Integer,\"Number of forward-strand reference and alternate allele observations\">" << endl
        << "##FORMAT=<ID=BCR,Number=2,Type=Integer,\"Number of reverse-strand reference and alternate allele observations\">" << endl
        << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (vector<string>::iterator s = samples.begin(); s != samples.end(); ++s) {
        out << "\t" << *s;
    }
    out << endl;

}

string vcf(
        long double comboProb,
        //long double alleleSamplingProb,
        Samples& samples,
        string refbase,
        string altbase,
        Allele& altAllele,
        vector<string>& sampleNames,
        int coverage,
        GenotypeCombo& genotypeCombo,
        Results& results,
        AlleleParser* parser) {

    stringstream out;

    // TODO make it so you use the genotypeCombo... 

    GenotypeComboMap comboMap;
    genotypeCombo2Map(genotypeCombo, comboMap);

    int samplesWithData = 0;
    // count alternate alleles in the best genotyping
    int alternateCount = 0;
    int alleleCount = 0;
    // reference / alternate base counts by strand
    map<string, pair<int, int> > altAndRefCountsBySample;
    int alternateObsCount = 0;
    int referenceObsCount = 0;
    // het counts
    int hetReferenceObsCount = 0;
    int hetAlternateObsCount = 0;
    int hetAllObsCount = 0;
    pair<int, int> baseCountsForwardTotal = make_pair(0, 0);
    pair<int, int> baseCountsReverseTotal = make_pair(0, 0);
    map<string, pair<int, int> > baseCountsForwardBySample;
    map<string, pair<int, int> > baseCountsReverseBySample;
    for (vector<string>::iterator sampleName = sampleNames.begin(); sampleName != sampleNames.end(); ++sampleName) {
        GenotypeComboMap::iterator gc = comboMap.find(*sampleName);
        //cerr << "alternate count for " << altbase << " and " << *genotype << " is " << genotype->alleleCount(altbase) << endl;
        if (gc != comboMap.end()) {
            Genotype* genotype = gc->second.first;

            Sample& sample = samples[*sampleName];

            // check that we actually have observations for this sample
            int observationCount = sample.observationCount();
            if (observationCount == 0) {
                continue;
            } else {
                ++samplesWithData;
            }

            alternateCount += genotype->alleleFrequency(altbase);
            alleleCount += genotype->ploidy;

            if (!genotype->homozygous) {
                hetAllObsCount += observationCount;
                hetReferenceObsCount += sample[refbase].size();
                hetAlternateObsCount += sample[altbase].size();
            }

            pair<int, int> altAndRefCounts = make_pair(sample.observationCount(altbase), sample.observationCount(refbase));
            altAndRefCountsBySample[*sampleName] = altAndRefCounts;
            alternateObsCount += altAndRefCounts.first;
            referenceObsCount += altAndRefCounts.second;

            pair<pair<int,int>, pair<int,int> > baseCounts = sample.baseCount(refbase, altbase);
            baseCountsForwardBySample[*sampleName] = baseCounts.first;
            baseCountsReverseBySample[*sampleName] = baseCounts.second;
            baseCountsForwardTotal.first += baseCounts.first.first;
            baseCountsForwardTotal.second += baseCounts.first.second;
            baseCountsReverseTotal.first += baseCounts.second.first;
            baseCountsReverseTotal.second += baseCounts.second.second;
        }
    }

    int allObsCount = alternateObsCount + referenceObsCount;

    string referenceSequence;
    string alternateSequence;
    switch (altAllele.type) {
        case ALLELE_SNP:
            referenceSequence = refbase;
            alternateSequence = altAllele.alternateSequence;
            break;
        case ALLELE_MNP:
            referenceSequence = parser->referenceSubstr(parser->currentPosition, altAllele.length);
            alternateSequence = altAllele.alternateSequence;
            break;
        case ALLELE_DELETION:
            referenceSequence = parser->referenceSubstr(parser->currentPosition, altAllele.length + 1);
            alternateSequence = refbase;
            break;
        case ALLELE_INSERTION:
            referenceSequence = refbase;
            alternateSequence = refbase + altAllele.alternateSequence.substr(1); // strip leading "I"
            break;
        default:
            cerr << "Unhandled allele type: " << altAllele.typeStr() << endl;
            break;
    }

    //string refbase = parser->currentReferenceBase();
    // positional information
    // CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT
    out << parser->currentTarget->seq << "\t"
        << (long unsigned int) parser->currentPosition + 1 << "\t"
        << "." << "\t"
        << referenceSequence << "\t"
        << alternateSequence << "\t"
        << float2phred(1 - comboProb) << "\t"
        << "." << "\t" // filter, no filter applied
        << "NS=" << samplesWithData << ";"
        << "DP=" << coverage << ";"
        << "AC=" << alternateCount << ";"
        << "AN=" << alleleCount << ";"
        //<< "AF=" << (double) alternateCount / (double) alleleCount << ";" // estimated alternate allele frequency in the range (0,1]
        // strand specific base counts, forward strand, reference and alternate, colon separated, comma separated for each alternate
        << "BCF=" << baseCountsForwardTotal.first << "," << baseCountsForwardTotal.second << ";"
        // strand specific base counts, reverse strand, reference and alternate, colon separated, comma separated for each alternate
        << "BCR=" << baseCountsReverseTotal.first << "," << baseCountsReverseTotal.second << ";"
        // strand bias for the alternate allele, a number between 0 and 1 representing the (weighted) ratio 
        // of forward strand sequence reads showing the alternate
        // allele to all sequence reads showing the alternate allele
        << "SB=" << (double) baseCountsForwardTotal.second / (double) alternateObsCount << ";"
        // allele balance at heterozygous sites, a number between 0 and 1 representing the weighted ratio of 
        // reads showing the reference allele to all reads,
        // restricted to sequence data from individuals who appear
        // to be heterozygous at this site.
        << "AB=" << ((hetAllObsCount == 0) ? 0 : (double) hetReferenceObsCount / (double) hetAllObsCount ) << ";"
        // allele balance counts, two numbers, comma separated, giving the numbers of sequence reads
        // from apparent heterozygotes which show the reference and
        // alternate alleles for this site.  ignores mapping strand
        // information.  these could be totals of numbers from the
        // bc fields for individual genotypes, counting only
        // individuals who are heterozygous.
        << "ABB=" << hetReferenceObsCount << "," << hetAlternateObsCount <<  ";"
        //<< "PRSQ=;" // TODO predicted R-squared number between 0 and 1 which estimates the correlation at 
                    // each site between imputed genotypes and the actual
                    // alternate allele counts that would be found if
                    // experimental genotyping were done on this panel of
                    // individuals
        //<< "AVPR=;" // TODO average posterior probability for genotype calls,
                    // average across all individuals of the posterior probability for the
                    // most probable discrete genotype call
        // homopolymer run length.  number of consecutive nucleotides (prior to this position?) in the genome
        // reference sequence matching the alternate allele, after substituting the
        // alternate in place of the reference sequence allele
        << "RUN=" << parser->homopolymerRunLeft(altbase) + 1 + parser->homopolymerRunRight(altbase) << ";";

    // allele class
    if (altAllele.type == ALLELE_DELETION) {
        out << "DEL=" << altAllele.length;
    } else if (altAllele.type == ALLELE_INSERTION) {
        out << "INS=" << altAllele.length;
    } else if (altAllele.type == ALLELE_SNP) {
        out << "SNP";
        // ts/tv
        if (isTransition(refbase, altbase)) {
            out << ";TS";
        } else {
            out << ";TV";
        }

        // CpG
        if (parser->isCpG(altbase)) {
            out << ";CpG";
        }
    } else if (altAllele.type == ALLELE_MNP) {
        out << "MNP=" << altAllele.length;
    }


    out << "\t" << "GT:GQ:DP:RA:AA:BCF:BCR";
    // TODO GL, un-normalized data likelihoods for genotypes

    // samples
    for (vector<string>::iterator sampleName = sampleNames.begin(); sampleName != sampleNames.end(); ++sampleName) {
        GenotypeComboMap::iterator gc = comboMap.find(*sampleName);
        Results::iterator s = results.find(*sampleName);
        if (gc != comboMap.end() && s != results.end()) {
            ResultData& sample = s->second;
            Genotype* genotype = gc->second.first;
            pair<int, int> altAndRefCounts = altAndRefCountsBySample[*sampleName]; // alternateAndReferenceCount(sample.observations, refbase, altbase);
            out << "\t"
                << genotype->relativeGenotype(refbase, altbase)
                << ":" << float2phred(1 - safe_exp(sample.marginals[genotype]))
                << ":" << sample.observations->observationCount()
                << ":" << altAndRefCounts.second
                << ":" << altAndRefCounts.first
                << ":" << sample.observations->baseCount(refbase, STRAND_FORWARD) << "," << sample.observations->baseCount(altbase, STRAND_FORWARD) // TODO
                << ":" << sample.observations->baseCount(refbase, STRAND_REVERSE) << "," << sample.observations->baseCount(altbase, STRAND_REVERSE) // TODO
                ;
                //<< ":" << "GL"  // TODO
        } else {
            out << "\t.";
        }
    }

    return out.str();
}

pair<Genotype*, long double> ResultData::bestMarginalGenotype(void) {
    map<Genotype*, long double>::iterator g = marginals.begin();
    pair<Genotype*, long double> best = make_pair(g->first, g->second); ++g;
    for ( ; g != marginals.end() ; ++g) {
        if (g->second > best.second) {
            best = make_pair(g->first, g->second);
        }
    }
    return best;
}

// TODO vcf output


