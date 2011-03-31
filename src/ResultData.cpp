#include "ResultData.h"
#include "TryCatch.h"

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
            if (g != sample.marginals.begin()) out << ",";
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

        // allele frequency metrics
        << "##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Total number of alternate alleles in called genotypes\">" << endl
        << "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">" << endl
        << "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Estimated allele frequency in the range (0,1]\">" << endl
        << "##INFO=<ID=HETAR,Number=1,Type=Integer,Description=\"Number of individuals heterozygous alternate / reference\">" << endl
        << "##INFO=<ID=HOMA,Number=1,Type=Integer,Description=\"Number of individuals homozygous for the alternate\">" << endl
        << "##INFO=<ID=HOMR,Number=1,Type=Integer,Description=\"Number of individuals homozygous for the reference\">" << endl

        // binomial balance metrics
        << "##INFO=<ID=RA,Number=1,Type=Integer,Description=\"Reference allele observations\">" << endl
        << "##INFO=<ID=AA,Number=1,Type=Integer,Description=\"Alternate allele observations\">" << endl
        << "##INFO=<ID=SRF,Number=1,Type=Integer,Description=\"Number of reference observations on the forward strand\">" << endl
        << "##INFO=<ID=SRR,Number=1,Type=Integer,Description=\"Number of reference observations on the reverse strand\">" << endl
        << "##INFO=<ID=SAF,Number=1,Type=Integer,Description=\"Number of alternate observations on the forward strand\">" << endl
        << "##INFO=<ID=SAR,Number=1,Type=Integer,Description=\"Number of alternate observations on the reverse strand\">" << endl
        << "##INFO=<ID=SRB,Number=1,Type=Float,Description=\"Strand bias for the reference allele: SRF / ( SRF + SRR )\">" << endl
        << "##INFO=<ID=SAB,Number=1,Type=Float,Description=\"Strand bias for the alternate allele: SAF / ( SAF + SAR )\">" << endl
        << "##INFO=<ID=SRP,Number=1,Type=Float,Description=\"Strand balance probability for the reference allele: Phred-scaled upper-bounds estimate of the probability of observing the deviation between SRF and SRR given E(SRF/SRR) ~ 0.5, derived using Hoeffding's inequality\">" << endl
        << "##INFO=<ID=SAP,Number=1,Type=Float,Description=\"Strand balance probability for the alternate allele: Phred-scaled upper-bounds estimate of the probability of observing the deviation between SAF and SAR given E(SAF/SAR) ~ 0.5, derived using Hoeffding's inequality\">" << endl
        << "##INFO=<ID=ABR,Number=1,Type=Integer,Description=\"Reference allele balance count: the number of sequence reads from apparent heterozygotes supporting the reference allele\">" << endl
        << "##INFO=<ID=ABA,Number=1,Type=Integer,Description=\"Alternate allele balance count: the number of sequence reads from apparent heterozygotes supporting the alternate allele\">" << endl
        << "##INFO=<ID=AB,Number=1,Type=Float,Description=\"Allele balance at heterozygous sites: a number between 0 and 1 representing the ratio of reads showing the reference allele to all reads, considering only reads from individuals called as heterozygous\">" << endl
        << "##INFO=<ID=ABP,Number=1,Type=Float,Description=\"Allele balance probability at heterozygous sites: Phred-scaled upper-bounds estimate of the probability of observing the deviation between ABR and ABA given E(ABR/ABA) ~ 0.5, derived using Hoeffding's inequality\">" << endl
        << "##INFO=<ID=RUN,Number=1,Type=Integer,Description=\"Homopolymer run length: the number of consecutive nucleotides in the reference genome matching the alternate allele prior to the current position\">" << endl
        << "##INFO=<ID=RL,Number=1,Type=Integer,Description=\"Reads Placed Left: number of reads supporting the alternate balanced to the left (5') of the alternate allele\">" << endl
        << "##INFO=<ID=RR,Number=1,Type=Integer,Description=\"Reads Placed Right: number of reads supporting the alternate balanced to the right (3') of the alternate allele\">" << endl
        << "##INFO=<ID=RPP,Number=1,Type=Float,Description=\"Read Placement Probability: Phred-scaled upper-bounds estimate of the probability of observing the deviation between RPL and RPR given E(RPL/RPR) ~ 0.5, derived using Hoeffding's inequality\">" << endl
        << "##INFO=<ID=EL,Number=1,Type=Integer,Description=\"Allele End Left: number of observations of the alternate where the alternate occurs in the left end of the read\">" << endl
        << "##INFO=<ID=ER,Number=1,Type=Integer,Description=\"Allele End Right: number of observations of the alternate where the alternate occurs in the right end of the read\">" << endl
        << "##INFO=<ID=EPP,Number=1,Type=Float,Description=\"End Placement Probability: Phred-scaled upper-bounds estimate of the probability of observing the deviation between EL and ER given E(EL/ER) ~ 0.5, derived using Hoeffding's inequality\">" << endl
        << "##INFO=<ID=BL,Number=1,Type=Integer,Description=\"Base Pairs Left: number of base pairs in reads supporting the alternate to the left (5') of the alternate allele\">" << endl
        << "##INFO=<ID=BR,Number=1,Type=Integer,Description=\"Base Pairs Right: number of base pairs in reads supporting the alternate to the right (3') of the alternate allele\">" << endl
        << "##INFO=<ID=LRB,Number=1,Type=Float,Description=\"((max(BR, BL) / (BR + BL)) - 0.5) * 2 : The proportion of base pairs in reads on one side of the alternate allele relative to total bases, scaled from [0.5,1] to [0,1]\">" << endl
        << "##INFO=<ID=LRBP,Number=1,Type=Float,Description=\"Left-Right Balance Probability: Phred-scaled upper-bounds estimate of the probability of observing the deviation between BL and BR given E(BR/BL) ~ 0.5, derived using Hoeffding's inequality\">" << endl

        // supplementary information about the site
        << "##INFO=<ID=BVAR,Number=0,Type=Flag,Description=\"The best genotype combination in the posterior is variant (non homozygous).\">" << endl
        << "##INFO=<ID=SNP,Number=0,Type=Flag,Description=\"SNP allele\">" << endl
        << "##INFO=<ID=TS,Number=0,Type=Flag,Description=\"transition SNP\">" << endl
        << "##INFO=<ID=TV,Number=0,Type=Flag,Description=\"transversion SNP\">" << endl
        << "##INFO=<ID=CpG,Number=0,Type=Flag,Description=\"CpG site (either CpG, TpG or CpA)\">" << endl
        << "##INFO=<ID=MNP,Number=0,Type=Flag,Description=\"MNP allele\">" << endl
        << "##INFO=<ID=INS,Number=0,Type=Flag,Description=\"insertion allele\">" << endl
        << "##INFO=<ID=DEL,Number=0,Type=Flag,Description=\"deletion allele\">" << endl
        << "##INFO=<ID=LEN,Number=1,Type=Integer,Description=\"allele length\">" << endl
        << "##INFO=<ID=MQM,Number=1,Type=Float,Description=\"Mean mapping quality of observed alternate alleles\">" << endl;
    if (parameters.showReferenceRepeats) {
        out << "##INFO=<ID=REPEAT,Number=1,Type=String,Description=\"Description of the local repeat structures flanking the current position\">" << endl;
    }

        // format fields for genotypes
    out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl
        << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype\">" << endl
        << "##FORMAT=<ID=GL,Number=1,Type=String,Description=\"Genotype Likelihood, log-scaled likeilhoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy\">" << endl
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
        map<string, vector<Allele*> >& alleleGroups,
        map<int, vector<Genotype> >& genotypesByPloidy,
        Results& results,
        AlleleParser* parser) {

    stringstream out;

    Parameters& parameters = parser->parameters;

    // TODO make it so you use the genotypeCombo... 

    GenotypeComboMap comboMap;
    genotypeCombo2Map(genotypeCombo, comboMap);

    map<int, vector<Genotype*> > presentGenotypesByPloidy;
    for (map<int, vector<Genotype> >::iterator gs = genotypesByPloidy.begin(); gs != genotypesByPloidy.end(); ++gs) {
        int ploidy = gs->first;
        vector<Genotype>& genotypes = gs->second;
        for (vector<Genotype>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
            Genotype* genotype = &*g;
            if (genotype->alleleFrequency(refbase) + genotype->alleleFrequency(altbase) == genotype->ploidy) {
                presentGenotypesByPloidy[ploidy].push_back(genotype);
            }
        }
    }

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
    int altAlleleObservations = 0;
    int refAlleleObservations = 0;
    int hetAltRefSamples = 0;
    int homAltSamples = 0;
    int homRefSamples = 0;

    pair<int, int> baseCountsForwardTotal = make_pair(0, 0);
    pair<int, int> baseCountsReverseTotal = make_pair(0, 0);
    map<string, pair<int, int> > baseCountsForwardBySample;
    map<string, pair<int, int> > baseCountsReverseBySample;
    for (vector<string>::iterator sampleName = sampleNames.begin(); sampleName != sampleNames.end(); ++sampleName) {
        GenotypeComboMap::iterator gc = comboMap.find(*sampleName);
        //cerr << "alternate count for " << altbase << " and " << *genotype << " is " << genotype->alleleCount(altbase) << endl;
        if (gc != comboMap.end()) {
            Genotype* genotype = gc->second->genotype;

            Sample& sample = *gc->second->sample;

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
                hetReferenceObsCount += sample.observationCount(refbase);
                hetAlternateObsCount += sample.observationCount(altbase);
                if (hetAlternateObsCount > 0) {
                    ++hetAltRefSamples;
                }
            } else {
                if (genotype->alleleFrequency(refbase) > 0) {
                    ++homRefSamples;
                } else {
                    ++homAltSamples;
                }
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
    int hetAllObsCount = hetReferenceObsCount + hetAlternateObsCount;

    unsigned int basesLeft = 0;
    unsigned int basesRight = 0;
    unsigned int readsLeft = 0;
    unsigned int readsRight = 0;
    unsigned int endLeft = 0;
    unsigned int endRight = 0;

    unsigned int mqsum = 0;
    vector<Allele*>& alternateAlleles = alleleGroups.at(altbase);
    for (vector<Allele*>::iterator app = alternateAlleles.begin(); app != alternateAlleles.end(); ++app) {
        Allele& allele = **app;
        basesLeft += allele.basesLeft;
        basesRight += allele.basesRight;
        if (allele.basesLeft >= allele.basesRight) {
            readsLeft += 1;
            if (allele.strand == STRAND_FORWARD) {
                endLeft += 1;
            } else {
                endRight += 1;
            }
        } else {
            readsRight += 1;
            if (allele.strand == STRAND_FORWARD) {
                endRight += 1;
            } else {
                endLeft += 1;
            }
        }
        mqsum += allele.mapQuality;
    }

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
            TRY { alternateSequence = refbase + altAllele.alternateSequence.substr(1); // strip leading "I"
            } CATCH;
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
        << "AF=" << ((alleleCount == 0) ? 0 : (double) alternateCount / (double) alleleCount) << ";"
        << "RA=" << refAlleleObservations << ";"
        << "AA=" << altAlleleObservations << ";"
        << "HETAR=" << hetAltRefSamples << ";"
        << "HOMA=" << homAltSamples << ";"
        << "HOMR=" << homRefSamples << ";"
        << "SRF=" << baseCountsForwardTotal.first << ";"
        << "SRR=" << baseCountsReverseTotal.first << ";"
        << "SAF=" << baseCountsForwardTotal.second << ";"
        << "SAR=" << baseCountsReverseTotal.second << ";"
        << "SRB=" << ((referenceObsCount == 0) ? 0 : (double) baseCountsForwardTotal.first / (double) referenceObsCount) << ";"
        << "SAB=" << ((alternateObsCount == 0) ? 0 : (double) baseCountsForwardTotal.second / (double) alternateObsCount) << ";"
        << "SRP=" << ((referenceObsCount == 0) ? 0 : ln2phred(hoeffdingln(baseCountsForwardTotal.first, referenceObsCount, 0.5))) << ";"
        << "SAP=" << ((alternateObsCount == 0) ? 0 : ln2phred(hoeffdingln(baseCountsForwardTotal.second, alternateObsCount, 0.5))) << ";"
        << "ABR=" << hetReferenceObsCount <<  ";"
        << "ABA=" << hetAlternateObsCount <<  ";"
        << "AB="  << ((hetAllObsCount == 0) ? 0 : (double) hetReferenceObsCount / (double) hetAllObsCount ) << ";"
        << "ABP=" << ((hetAllObsCount == 0) ? 0 : ln2phred(hoeffdingln(hetReferenceObsCount, hetAllObsCount, 0.5))) << ";"
        << "RUN=" << parser->homopolymerRunLeft(altbase) + 1 + parser->homopolymerRunRight(altbase) << ";"
        << "MQM=" << ((alternateAlleles.size() == 0) ? 0 : (double) mqsum / (double) alternateAlleles.size()) << ";"
        << "RL=" << readsLeft << ";"
        << "RR=" << readsRight << ";"
        << "RPP=" << ln2phred(hoeffdingln(readsLeft, readsRight + readsLeft, 0.5)) << ";"  // estimates upper bound for the lower tail of the binomial distribution
        << "EL=" << endLeft << ";"
        << "ER=" << endRight << ";"
        << "EPP=" << ((basesLeft + basesRight == 0) ? 0 : ln2phred(hoeffdingln(endLeft, endLeft + endRight, 0.5))) << ";"
        << "BL=" << basesLeft << ";"
        << "BR=" << basesRight << ";"
        << "LRB=" << ((double) max(basesLeft, basesRight) / (double) (basesRight + basesLeft) - 0.5) * 2 << ";"
        << "LRBP=" << ((basesLeft + basesRight == 0) ? 0 : ln2phred(hoeffdingln(basesLeft, basesLeft + basesRight, 0.5))) << ";";

    if (bestOverallComboIsHet) {
        out << "BVAR;";
    }

    if (parameters.showReferenceRepeats && !repeats.empty()) {
        stringstream repeatsstr;
        for (map<string, int>::iterator c = repeats.begin(); c != repeats.end(); ++c) {
            repeatsstr << c->first << ":" << c->second << "|";
        }
        string repeatstr = repeatsstr.str();
        TRY { repeatstr = repeatstr.substr(0, repeatstr.size() - 1); } CATCH;
        out << "REPEAT=" << repeatstr << ";";
    }

    // allele class
    if (altAllele.type == ALLELE_DELETION) {
        out << "DEL";
        // what is the class of deletion
        // microsatellite repeat?
        // "novel"?
        // how large is the repeat, if there is one?
    } else if (altAllele.type == ALLELE_INSERTION) {
        out << "INS";
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
        out << "MNP";
    }
    out << ";LEN=" << altAllele.length;


    out << "\t" << "GT:" << (parameters.calculateMarginals ? "GQ:" : "") << "GL:DP:RA:AA:SR:SA";

    // samples
    for (vector<string>::iterator sampleName = sampleNames.begin(); sampleName != sampleNames.end(); ++sampleName) {
        GenotypeComboMap::iterator gc = comboMap.find(*sampleName);
        Results::iterator s = results.find(*sampleName);
        if (gc != comboMap.end() && s != results.end()) {
            ResultData& sample = s->second;
            Genotype* genotype = gc->second->genotype;
            pair<int, int> altAndRefCounts = altAndRefCountsBySample[*sampleName]; // alternateAndReferenceCount(sample.observations, refbase, altbase);
            out << "\t"
                << genotype->relativeGenotype(refbase, altbase);
            if (parameters.calculateMarginals) {
                out << ":" << float2phred(1 - safe_exp(sample.marginals[genotype]));
            }

            // get data likelihoods for present genotypes, none if we have excluded genotypes from data likelihood calculations
            stringstream datalikelihoods;
            if (!parameters.excludeUnobservedGenotypes && !parameters.excludePartiallyObservedGenotypes) {
                vector<Genotype*>& presentGenotypes = presentGenotypesByPloidy[parser->currentSamplePloidy(*sampleName)];
                for (vector<Genotype*>::iterator g = presentGenotypes.begin(); g != presentGenotypes.end(); ++g) {
                    datalikelihoods << ((g == presentGenotypes.begin()) ? "" : ",") << (*g)->slashstr() << "=" << sample.genotypeLikelihood(*g);
                }
            }
            out << ":" << datalikelihoods.str()
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


