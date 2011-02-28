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
        << "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Estimated allele frequency in the range (0,1]\">" << endl
        << "##INFO=<ID=SR,Number=1,Type=Integer,Description=\"Number of reference observations by strand, delimited by |: [forward]|[reverse]\">" << endl
        << "##INFO=<ID=SA,Number=1,Type=Integer,Description=\"Number of alternate observations by strand, delimited by |: [forward]|[reverse]\">" << endl
        << "##INFO=<ID=RA,Number=1,Type=Integer,Description=\"Reference allele observations\">" << endl
        << "##INFO=<ID=AA,Number=1,Type=Integer,Description=\"Alternate allele observations\">" << endl
        << "##INFO=<ID=SB,Number=1,Type=Float,Description=\"Strand bias for the alternate allele: a number between 0 and 1 representing the ratio of forward strand sequence reads showing the alternate allele to all reads, considering only reads from individuals called as heterozygous\">" << endl
        << "##INFO=<ID=AB,Number=1,Type=Float,Description=\"Allele balance at heterozygous sites: a number between 0 and 1 representing the ratio of reads showing the reference allele to all reads, considering only reads from individuals called as heterozygous\">" << endl
        << "##INFO=<ID=ABR,Number=1,Type=Integer,Description=\"Reference allele balance count: the number of sequence reads from apparent heterozygotes supporting the reference allele\">" << endl
        << "##INFO=<ID=ABA,Number=1,Type=Integer,Description=\"Alternate allele balance count: the number of sequence reads from apparent heterozygotes supporting the alternate allele\">" << endl
        << "##INFO=<ID=RUN,Number=1,Type=Integer,Description=\"Homopolymer run length: the number of consecutive nucleotides in the reference genome matching the alternate allele prior to the current position\">" << endl
        << "##INFO=<ID=BVAR,Number=0,Type=Flag,Description=\"The best genotype combination in the posterior is variant (non homozygous).\">" << endl
        << "##INFO=<ID=SNP,Number=0,Type=Flag,Description=\"SNP allele\">" << endl
        << "##INFO=<ID=TS,Number=0,Type=Flag,Description=\"transition SNP\">" << endl
        << "##INFO=<ID=TV,Number=0,Type=Flag,Description=\"transversion SNP\">" << endl
        << "##INFO=<ID=CpG,Number=0,Type=Flag,Description=\"CpG site (either CpG, TpG or CpA)\">" << endl
        << "##INFO=<ID=MNP,Number=0,Type=Integer,Description=\"Length of MNP allele, if present\">" << endl
        << "##INFO=<ID=INS,Number=1,Type=Integer,Description=\"Length of insertion allele, if present\">" << endl
        << "##INFO=<ID=DEL,Number=1,Type=Integer,Description=\"Length of deletion allele, if present\">" << endl
        << "##INFO=<ID=REPEAT,Number=1,Type=String,Description=\"Description of the local repeat structures flanking the current position\">" << endl
        << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl
        << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality, the PHRED-scaled marginal (or unconditional) probability of the called genotype\">" << endl
        << "##FORMAT=<ID=GL,Number=1,Type=Float,Description=\"Genotype Likelihood, log-scaled likeilhood of the data given the called genotype\">" << endl
        << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">" << endl
        << "##FORMAT=<ID=RA,Number=1,Type=Integer,Description=\"Reference allele observations\">" << endl
        << "##FORMAT=<ID=AA,Number=1,Type=Integer,Description=\"Alternate allele observations\">" << endl
        << "##FORMAT=<ID=SR,Number=1,Type=Integer,Description=\"Number of reference observations by strand, delimited by |: [forward]|[reverse]\">" << endl
        << "##FORMAT=<ID=SA,Number=1,Type=Integer,Description=\"Number of alternate observations by strand, delimited by |: [forward]|[reverse]\">" << endl
        << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (vector<string>::iterator s = samples.begin(); s != samples.end(); ++s) {
        out << "\t" << *s;
    }
    out << endl;

}

string vcf(
        long double pHom,
        //long double alleleSamplingProb,
        Samples& samples,
        string refbase,
        string altbase,
        Allele& altAllele,
        map<string, int> repeats,
        vector<string>& sampleNames,
        int coverage,
        GenotypeCombo& genotypeCombo,
        bool bestOverallComboIsHet,
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
    int altAlleleObservations = 0;
    int refAlleleObservations = 0;
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

            // TODO cleanup
            pair<pair<int,int>, pair<int,int> > baseCounts = sample.baseCount(refbase, altbase);
            baseCountsForwardBySample[*sampleName] = baseCounts.first;
            baseCountsReverseBySample[*sampleName] = baseCounts.second;
            baseCountsForwardTotal.first += baseCounts.first.first;
            baseCountsForwardTotal.second += baseCounts.first.second;
            baseCountsReverseTotal.first += baseCounts.second.first;
            baseCountsReverseTotal.second += baseCounts.second.second;
            refAlleleObservations += baseCounts.first.first + baseCounts.second.first;
            altAlleleObservations += baseCounts.first.second + baseCounts.second.second;
        }
    }

    int allObsCount = alternateObsCount + referenceObsCount;

    // 0-based variant position
    long unsigned int variantPosition = (long unsigned int) parser->currentPosition;

    string referenceSequence;
    string alternateSequence;
    switch (altAllele.type) {
        case ALLELE_SNP:
            referenceSequence = refbase;
            alternateSequence = altAllele.alternateSequence;
            break;
        case ALLELE_MNP:
            referenceSequence = parser->referenceSubstr(variantPosition, altAllele.length);
            alternateSequence = altAllele.alternateSequence;
            break;
        case ALLELE_DELETION:
            referenceSequence = parser->referenceSubstr(variantPosition - 1, altAllele.length + 1);
            // this decrement fixes deletion position reporting to match VCF
            // spec, in which we are reporting the position of the base prior
            // to the deletion
            --variantPosition;
            alternateSequence = referenceSequence.at(0);
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
    //out.setf(ios::fixed,ios::floatfield);
    out.precision(5);
    out << parser->currentSequenceName << "\t"
        << variantPosition + 1 << "\t"
        << "." << "\t"
        << referenceSequence << "\t"
        << alternateSequence << "\t"
        << float2phred(pHom) << "\t";
    out.precision(5);
    out << "." << "\t" // filter, no filter applied
        << "NS=" << samplesWithData << ";"
        << "DP=" << coverage << ";"
        << "AC=" << alternateCount << ";"
        << "AN=" << alleleCount << ";"
        << "AF=" << (double) alternateCount / (double) alleleCount << ";"
        << "RA=" << refAlleleObservations << ";"
        << "AA=" << altAlleleObservations << ";"
        << "SR=" << baseCountsForwardTotal.first << "|" << baseCountsReverseTotal.first << ";"
        << "SA=" << baseCountsForwardTotal.second << "|" << baseCountsReverseTotal.second << ";"
        << "SB=" << (double) baseCountsForwardTotal.second / (double) alternateObsCount << ";"
        << "AB=" << ((hetAllObsCount == 0) ? 0 : (double) hetReferenceObsCount / (double) hetAllObsCount ) << ";"
        << "ABR=" << hetReferenceObsCount <<  ";"
        << "ABA=" << hetAlternateObsCount <<  ";"
        << "RUN=" << parser->homopolymerRunLeft(altbase) + 1 + parser->homopolymerRunRight(altbase) << ";";

    if (bestOverallComboIsHet) {
        out << "BVAR;";
    }

    if (!repeats.empty()) {
        stringstream repeatsstr;
        for (map<string, int>::iterator c = repeats.begin(); c != repeats.end(); ++c) {
            repeatsstr << c->first << ":" << c->second << "|";
        }
        string repeatstr = repeatsstr.str();
        repeatstr = repeatstr.substr(0, repeatstr.size() - 1);
        out << "REPEAT=" << repeatstr << ";";
    }

    // allele class
    if (altAllele.type == ALLELE_DELETION) {
        out << "DEL=" << altAllele.length;
        // what is the class of deletion
        // microsatellite repeat?
        // "novel"?
        // how large is the repeat, if there is one?
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


    out << "\t" << "GT:GQ:GL:DP:RA:AA:SR:SA";
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
                << ":" << sample.genotypeLikelihood(genotype)
                << ":" << sample.observations->observationCount()
                << ":" << altAndRefCounts.second
                << ":" << altAndRefCounts.first
                << ":" << sample.observations->baseCount(refbase, STRAND_FORWARD) << "|" << sample.observations->baseCount(refbase, STRAND_REVERSE) 
                << ":" << sample.observations->baseCount(altbase, STRAND_FORWARD) << "|" << sample.observations->baseCount(altbase, STRAND_REVERSE) // TODO
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


