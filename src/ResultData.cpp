#include "ResultData.h"
#include "TryCatch.h"

using namespace std;



vcf::Variant& Results::vcf(
        vcf::Variant& var, // variant to update
        long double pHom,
        long double bestComboOddsRatio,
        //long double alleleSamplingProb,
        Samples& samples,
        string refbase,
        vector<Allele>& altAllelesIncludingNulls,
        map<string, int> repeats,
        vector<string>& sampleNames,
        int coverage,
        GenotypeCombo& genotypeCombo,
        bool bestOverallComboIsHet,
        map<string, vector<Allele*> >& alleleGroups,
        map<int, vector<Genotype> >& genotypesByPloidy,
        vector<string>& sequencingTechnologies,
        AlleleParser* parser) {

    Parameters& parameters = parser->parameters;

    GenotypeComboMap comboMap;
    genotypeCombo2Map(genotypeCombo, comboMap);

    // set up the reported reference allele
    long int referencePosition = (long int) parser->currentPosition; // 0-based

    // remove alt alleles
    vector<Allele> altAlleles;
    for (vector<Allele>::iterator aa = altAllelesIncludingNulls.begin(); aa != altAllelesIncludingNulls.end(); ++aa) {
        if (!aa->isNull()) {
            altAlleles.push_back(*aa);
        }
    }

    // adjust reference position, reference sequence, and alt alleles by
    // stripping invariant bases off the beginning and end of the alt alleles
    // first we find the minimum start and end matches
    vector<Allele> adjustedAltAlleles;
    int minStartMatch = 0;
    int minEndMatch = 0;
    for (vector<Allele>::iterator aa = altAlleles.begin(); aa != altAlleles.end(); ++aa) {
        vector<pair<int, string> > cigar = splitCigar(aa->cigar);
        int startmatch = 0;
        int endmatch = 0;
        if (cigar.front().second == "M") {
            startmatch = cigar.front().first;
        }
        if (cigar.back().second == "M") {
            endmatch = cigar.back().first;
        }
        if (cigar.size() > 1 && (cigar.at(1).second == "D" || cigar.at(1).second == "I")) {
            if (startmatch == 1) { // require at least one base flanking deletions
                startmatch = 0;
            }
        }
        if (aa == altAlleles.begin()) {
            minStartMatch = startmatch;
            minEndMatch = endmatch;
        } else {
            minStartMatch = min(minStartMatch, startmatch);
            minEndMatch = min(minEndMatch, endmatch);
        }
    }
    // if either is non-zero, we have to adjust cigars and alternate sequences to be printed
    // this is done solely for reporting, so the altAlleles structure is used
    // for stats generation out of the ML genotype combination
    map<string, string> adjustedCigar;
    if (minStartMatch || minEndMatch) {
        for (vector<Allele>::iterator aa = altAlleles.begin(); aa != altAlleles.end(); ++aa) {
            // subtract the minStartMatch and minEndMatch bases from the allele start and end
            adjustedAltAlleles.push_back(*aa);
            Allele& allele = adjustedAltAlleles.back();
            vector<pair<int, string> > cigar = splitCigar(allele.cigar);
            // TODO clean this up by writing a wrapper for Allele::subtract() (?)
            if (cigar.front().second == "M") {
                cigar.front().first -= minStartMatch;
                allele.alternateSequence = allele.alternateSequence.substr(minStartMatch);
            }
            if (cigar.back().second == "M") {
                cigar.back().first -= minEndMatch;
                allele.alternateSequence = allele.alternateSequence.substr(0, allele.alternateSequence.size() - minEndMatch);
            }
            allele.cigar = joinCigar(cigar);
            allele.position += minStartMatch;
            allele.referenceLength -= minStartMatch + minEndMatch;
            adjustedCigar[aa->base()] = allele.cigar;
        }
        referencePosition += minStartMatch;
    } else {
        adjustedAltAlleles = altAlleles;
        for (vector<Allele>::iterator aa = altAlleles.begin(); aa != altAlleles.end(); ++aa) {
            adjustedCigar[aa->base()] = aa->cigar;
        }
    }

    refbase = parser->referenceSubstr(referencePosition, 1);
    var.ref = refbase;

    // the reference sequence should be able to encompass all events at the site, +1bp on the left
    for (vector<Allele>::iterator aa = adjustedAltAlleles.begin(); aa != adjustedAltAlleles.end(); ++aa) {

        Allele& altAllele = *aa;
        switch (altAllele.type) {
            case ALLELE_SNP:
            case ALLELE_REFERENCE:
            case ALLELE_MNP:
                if (var.ref.size() < altAllele.referenceLength) {
                    var.ref = parser->referenceSubstr(referencePosition, altAllele.referenceLength);
                }
                break;
            case ALLELE_DELETION:
                // extend the reference sequence
                if (var.ref.size() < altAllele.referenceLength) {
                    var.ref = parser->referenceSubstr(referencePosition, altAllele.referenceLength);
                }
                break;
            case ALLELE_INSERTION:
                break;
            case ALLELE_COMPLEX:
                if (var.ref.size() < altAllele.referenceLength) {
                    var.ref = parser->referenceSubstr(referencePosition, altAllele.referenceLength);
                }
                break;
            default:
                cerr << "Unhandled allele type: " << altAllele.typeStr() << endl;
                break;
        }

    }

    for (vector<Allele>::iterator aa = adjustedAltAlleles.begin(); aa != adjustedAltAlleles.end(); ++aa) {
        Allele& altAllele = *aa;
        string altseq;
        switch (altAllele.type) {
            case ALLELE_REFERENCE:
                break;
            case ALLELE_SNP:
            case ALLELE_MNP:
                altseq = var.ref;
                altseq.replace(0, altAllele.alternateSequence.size(), altAllele.alternateSequence);
                var.alt.push_back(altseq);
                break;
            case ALLELE_DELETION:
            case ALLELE_INSERTION: // XXX is this correct???
            case ALLELE_COMPLEX:
                var.alt.push_back(altAllele.alternateSequence);
                break;
            default:
                cerr << "Unhandled allele type: " << altAllele.typeStr() << endl;
                break;
        }
    }


    // get the required size of the reference sequence


    // set up VCF record-wide variables

    var.sequenceName = parser->currentSequenceName;
    // XXX this should be the position of the matching reference haplotype
    var.position = referencePosition + 1;
    var.id = ".";
    var.filter = ".";
    // XXX this should be the size of the maximum deletion + 1bp on the left end
    var.quality = float2phred(pHom);


    // set up format string

    var.format.clear();
    var.format.push_back("GT");
    if (parameters.calculateMarginals) var.format.push_back("GQ");
    // XXX
    var.format.push_back("DP");
    var.format.push_back("RO");
    var.format.push_back("QR");
    var.format.push_back("AO");
    var.format.push_back("QA");
    if (!parameters.excludeUnobservedGenotypes) {
        var.format.push_back("GL");
    }
    //var.format.push_back("GLE");

    unsigned int refBasesLeft = 0;
    unsigned int refBasesRight = 0;
    unsigned int refReadsLeft = 0;
    unsigned int refReadsRight = 0;
    unsigned int refEndLeft = 0;
    unsigned int refEndRight = 0;
    unsigned int refmqsum = 0;
    unsigned int refProperPairs = 0;
    long double refReadMismatchSum = 0;
    long double refReadSNPSum = 0;
    long double refReadIndelSum = 0;
    unsigned int refObsCount = 0;
    map<string, int> refObsBySequencingTechnology;

    map<string, vector<Allele*> >::iterator f = alleleGroups.find(refbase);
    if (f != alleleGroups.end()) {
        vector<Allele*>& referenceAlleles = alleleGroups.at(refbase);
        refObsCount = referenceAlleles.size();
        for (vector<Allele*>::iterator app = referenceAlleles.begin(); app != referenceAlleles.end(); ++app) {
            Allele& allele = **app;
            refReadMismatchSum += allele.readMismatchRate;
            refReadSNPSum += allele.readSNPRate;
            refReadIndelSum += allele.readIndelRate;
            if (allele.isProperPair) {
                ++refProperPairs;
            }
            if (!allele.sequencingTechnology.empty()) {
                ++refObsBySequencingTechnology[allele.sequencingTechnology];
            }
            refBasesLeft += allele.basesLeft;
            refBasesRight += allele.basesRight;
            if (allele.basesLeft >= allele.basesRight) {
                refReadsLeft += 1;
                if (allele.strand == STRAND_FORWARD) {
                    refEndLeft += 1;
                } else {
                    refEndRight += 1;
                }
            } else {
                refReadsRight += 1;
                if (allele.strand == STRAND_FORWARD) {
                    refEndRight += 1;
                } else {
                    refEndLeft += 1;
                }
            }
            refmqsum += allele.mapQuality;
        }
    }

    long double refReadMismatchRate = (refObsCount == 0 ? 0 : refReadMismatchSum / (long double) refObsCount);
    long double refReadSNPRate = (refObsCount == 0 ? 0 : refReadSNPSum / (long double) refObsCount);
    long double refReadIndelRate = (refObsCount == 0 ? 0 : refReadIndelSum / (long double) refObsCount);

    var.info["XRM"].push_back(convert(refReadMismatchRate));
    var.info["XRS"].push_back(convert(refReadSNPRate));
    var.info["XRI"].push_back(convert(refReadIndelRate));

    var.info["MQMR"].push_back(convert((refObsCount == 0) ? 0 : (double) refmqsum / (double) refObsCount));
    var.info["RPPR"].push_back(convert((refObsCount == 0) ? 0 : ln2phred(hoeffdingln(refReadsLeft, refReadsRight + refReadsLeft, 0.5))));
    var.info["EPPR"].push_back(convert((refBasesLeft + refBasesRight == 0) ? 0 : ln2phred(hoeffdingln(refEndLeft, refEndLeft + refEndRight, 0.5))));
    var.info["PAIREDR"].push_back(convert((refObsCount == 0) ? 0 : (double) refProperPairs / (double) refObsCount));

    var.info["HWE"].push_back(convert(ln2phred(genotypeCombo.hweComboProb())));

    // loop over all alternate alleles
    for (vector<Allele>::iterator aa = altAlleles.begin(); aa != altAlleles.end(); ++aa) {

        Allele& altAllele = *aa;
        string altbase = altAllele.base();

        // count alternate alleles in the best genotyping
        unsigned int alternateCount = 0;
        unsigned int alleleCount = 0;
        // reference / alternate base counts by strand
        map<string, unsigned int> altCountBySample;
        map<string, unsigned int> altQualBySample;
        // het counts
        unsigned int hetReferenceObsCount = 0;
        unsigned int hetOtherObsCount = 0;
        unsigned int hetAlternateObsCount = 0;
        unsigned int hetAltSamples = 0;
        unsigned int homAltSamples = 0;
        unsigned int homRefSamples = 0;
        unsigned int refSampleObsCount = 0; // depth in hom-ref samples
        unsigned int altSampleObsCount = 0; // depth in samples with called alternates
        // unique alternate alleles / all alternate alleles in alt-associated samples
        unsigned int uniqueAllelesInAltSamples = 0;
        //unsigned int hetAllObsCount = hetOtherObsCount + hetAlternateObsCount + hetReferenceObsCount;
        unsigned int hetAllObsCount = 0;

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
                unsigned int observationCount = sample.observationCount();
                if (observationCount == 0) {
                    continue;
                }

                alternateCount += genotype->alleleCount(altbase);
                alleleCount += genotype->ploidy;

                unsigned int altCount = sample.observationCount(altbase);
                unsigned int refCount = sample.observationCount(refbase);
                if (!genotype->homozygous) {
                    // het case
                    if (altCount > 0) {
                        ++hetAltSamples;
                        hetAllObsCount += observationCount;
                        hetReferenceObsCount += refCount;
                        hetOtherObsCount += observationCount - altCount;
                        hetAlternateObsCount += altCount;
                        altSampleObsCount += observationCount;
                        uniqueAllelesInAltSamples += sample.size();
                        if (refCount > 0) {
                            --uniqueAllelesInAltSamples; // ignore reference allele
                        }
                    }
                } else {
                    if (altCount > 0) {
                        ++homAltSamples;
                        altSampleObsCount += observationCount;
                        uniqueAllelesInAltSamples += sample.size();
                        if (refCount > 0) {
                            --uniqueAllelesInAltSamples; // ignore reference allele
                        }
                    } else {
                        ++homRefSamples;
                        refSampleObsCount += observationCount;
                    }
                }
                altCountBySample[*sampleName] = altCount;

                altQualBySample[*sampleName] = sample.qualSum(altbase);

                pair<pair<int,int>, pair<int,int> > baseCounts = sample.baseCount(refbase, altbase);
                baseCountsForwardBySample[*sampleName] = baseCounts.first;
                baseCountsReverseBySample[*sampleName] = baseCounts.second;
                baseCountsForwardTotal.first += baseCounts.first.first;
                baseCountsForwardTotal.second += baseCounts.first.second;
                baseCountsReverseTotal.first += baseCounts.second.first;
                baseCountsReverseTotal.second += baseCounts.second.second;
            }
        }

        unsigned int altBasesLeft = 0;
        unsigned int altBasesRight = 0;
        unsigned int altReadsLeft = 0;
        unsigned int altReadsRight = 0;
        unsigned int altEndLeft = 0;
        unsigned int altEndRight = 0;
        unsigned int altmqsum = 0;
        unsigned int altproperPairs = 0;
        long double altReadMismatchSum = 0;
        long double altReadSNPSum = 0;
        long double altReadIndelSum = 0;
        unsigned int altObsCount = 0;
        map<string, int> altObsBySequencingTechnology;

        map<string, vector<Allele*> >::iterator f = alleleGroups.find(altbase);
        if (f != alleleGroups.end()) {
            vector<Allele*>& alternateAlleles = alleleGroups.at(altbase);
            altObsCount = alternateAlleles.size();
            for (vector<Allele*>::iterator app = alternateAlleles.begin(); app != alternateAlleles.end(); ++app) {
                Allele& allele = **app;
                altReadMismatchSum += allele.readMismatchRate;
                altReadSNPSum += allele.readSNPRate;
                altReadIndelSum += allele.readIndelRate;
                if (allele.isProperPair) {
                    ++altproperPairs;
                }
                if (!allele.sequencingTechnology.empty()) {
                    ++altObsBySequencingTechnology[allele.sequencingTechnology];
                }
                altBasesLeft += allele.basesLeft;
                altBasesRight += allele.basesRight;
                if (allele.basesLeft >= allele.basesRight) {
                    altReadsLeft += 1;
                    if (allele.strand == STRAND_FORWARD) {
                        altEndLeft += 1;
                    } else {
                        altEndRight += 1;
                    }
                } else {
                    altReadsRight += 1;
                    if (allele.strand == STRAND_FORWARD) {
                        altEndRight += 1;
                    } else {
                        altEndLeft += 1;
                    }
                }
                altmqsum += allele.mapQuality;
            }
        } /*else {
            cerr << "couldn't find altbase: " << altbase << " in allele groups" << endl;
            for (map<string, vector<Allele*> >::iterator a = alleleGroups.begin(); a != alleleGroups.end(); ++a) {
                cerr << a->first << " " << a->second.size() << endl;
            }
            assert(false);
        }
        */

        long double altReadMismatchRate = (altObsCount == 0 ? 0 : altReadMismatchSum / altObsCount);
        long double altReadSNPRate = (altObsCount == 0 ? 0 : altReadSNPSum / altObsCount);
        long double altReadIndelRate = (altObsCount == 0 ? 0 : altReadIndelSum / altObsCount);
        
        var.info["XAM"].push_back(convert(altReadMismatchRate));
        var.info["XAS"].push_back(convert(altReadSNPRate));
        var.info["XAI"].push_back(convert(altReadIndelRate));

        // alt/ref ratios
        //var.info["ARM"].push_back(convert(refReadMismatchRate == 0 ? 0 : altReadMismatchRate / refReadMismatchRate));
        //var.info["ARS"].push_back(convert(refReadSNPRate == 0 ? 0 : altReadSNPRate / refReadSNPRate));
        //var.info["ARI"].push_back(convert(refReadIndelRate == 0 ? 0 : altReadIndelRate / refReadIndelRate));

        //string refbase = parser->currentReferenceBase();
        // positional information
        // CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT
        //out.setf(ios::fixed,ios::floatfield);
        //out.precision(5);

        var.info["AC"].push_back(convert(alternateCount));
        var.info["AN"].clear(); var.info["AN"].push_back(convert(alleleCount)); // XXX hack...
        var.info["AF"].push_back(convert((alleleCount == 0) ? 0 : (double) alternateCount / (double) alleleCount));
        var.info["AO"].push_back(convert(altObsCount));
        if (homRefSamples > 0 && hetAltSamples + homAltSamples > 0) {
            double altSampleAverageDepth = (double) altSampleObsCount
                       / ( (double) hetAltSamples + (double) homAltSamples );
            double refSampleAverageDepth = (double) refSampleObsCount / (double) homRefSamples;
            var.info["DPRA"].push_back(convert(altSampleAverageDepth / refSampleAverageDepth));
        } else {
            var.info["DPRA"].push_back(convert(0));
        }
            //<< "HETAR=" << hetAltSamples << ";"
            //<< "HOMA=" << homAltSamples << ";"
            //<< "HOMR=" << homRefSamples << ";"
        var.info["SRP"].clear(); // XXX hack
        var.info["SRP"].push_back(convert((refObsCount == 0) ? 0 : ln2phred(hoeffdingln(baseCountsForwardTotal.first, refObsCount, 0.5))));
        var.info["SAP"].push_back(convert((altObsCount == 0) ? 0 : ln2phred(hoeffdingln(baseCountsForwardTotal.second, altObsCount, 0.5))));
            //<< "ABR=" << hetReferenceObsCount <<  ";"
            //<< "ABA=" << hetAlternateObsCount <<  ";"
        var.info["AB"].push_back(convert((hetAllObsCount == 0) ? 0 : (double) hetAlternateObsCount / (double) hetAllObsCount ));
        var.info["ABP"].push_back(convert((hetAllObsCount == 0) ? 0 : ln2phred(hoeffdingln(hetAlternateObsCount, hetAllObsCount, 0.5))));
        var.info["RUN"].push_back(convert(parser->homopolymerRunLeft(altbase) + 1 + parser->homopolymerRunRight(altbase)));
        var.info["MQM"].push_back(convert((altObsCount == 0) ? 0 : (double) altmqsum / (double) altObsCount));
        var.info["RPP"].push_back(convert((altObsCount == 0) ? 0 : ln2phred(hoeffdingln(altReadsLeft, altReadsRight + altReadsLeft, 0.5))));
        var.info["EPP"].push_back(convert((altBasesLeft + altBasesRight == 0) ? 0 : ln2phred(hoeffdingln(altEndLeft, altEndLeft + altEndRight, 0.5))));
        var.info["PAIRED"].push_back(convert((altObsCount == 0) ? 0 : (double) altproperPairs / (double) altObsCount));
        var.info["CIGAR"].push_back(adjustedCigar[altAllele.base()]);
        var.info["MEANALT"].push_back(convert((hetAltSamples + homAltSamples == 0) ? 0 : (double) uniqueAllelesInAltSamples / (double) (hetAltSamples + homAltSamples)));

        for (vector<string>::iterator st = sequencingTechnologies.begin();
                st != sequencingTechnologies.end(); ++st) { string& tech = *st;
            var.info["technology." + tech].push_back(convert((altObsCount == 0) ? 0
                        : (double) altObsBySequencingTechnology[tech] / (double) altObsCount ));
        }

        if (bestOverallComboIsHet) {
            var.infoFlags["BVAR"] = true;
        }

        // allele class
        if (altAllele.type == ALLELE_DELETION) {
            var.info["TYPE"].push_back("del");
            // what is the class of deletion
            // microsatellite repeat?
            // "novel"?
            // how large is the repeat, if there is one?
        } else if (altAllele.type == ALLELE_INSERTION) {
            var.info["TYPE"].push_back("ins");
        } else if (altAllele.type == ALLELE_COMPLEX) {
            var.info["TYPE"].push_back("complex");
        } else if (altAllele.type == ALLELE_SNP) {
            var.info["TYPE"].push_back("snp");

            // CpG
            if (parser->isCpG(altbase)) {
                var.infoFlags["CpG"] = true;
            }
        } else if (altAllele.type == ALLELE_MNP) {
            var.info["TYPE"].push_back("mnp");
        } else {
            cout << "What is this?" << endl;
            cout << altAllele.type << endl;
            cout << altAllele << endl;
        }
        var.info["LEN"].push_back(convert(altAllele.length));


        // samples
        for (vector<string>::iterator sn = sampleNames.begin(); sn != sampleNames.end(); ++sn) {
            string& sampleName = *sn;
            GenotypeComboMap::iterator gc = comboMap.find(sampleName);
            Results::iterator s = find(sampleName);
            if (gc != comboMap.end() && s != end()) {
                Result& sample = s->second;
                Genotype* genotype = gc->second->genotype;

                map<string, vector<string> >& sampleOutput = var.samples[sampleName];

                sampleOutput["AO"].push_back(convert(altCountBySample[sampleName]));
                sampleOutput["QA"].push_back(convert(altQualBySample[sampleName]));

            }
        }

        // TODO
        // mismatch rate of reads containing supporting observations
        // vs. mismatch rate of reads without the alternate, for each alternate

    }


    // set up site-wide INFO tags, non-multiple

    // info variables

    // site-wide coverage
    int samplesWithData = 0;
    int refAlleleObservations = 0;
    for (vector<string>::iterator sampleName = sampleNames.begin(); sampleName != sampleNames.end(); ++sampleName) {
        GenotypeComboMap::iterator gc = comboMap.find(*sampleName);
        //cerr << "alternate count for " << altbase << " and " << *genotype << " is " << genotype->alleleCount(altbase) << endl;
        if (gc != comboMap.end()) {
            Genotype* genotype = gc->second->genotype;

            Sample& sample = *gc->second->sample;
            refAlleleObservations += sample.observationCount(refbase);

            ++samplesWithData;
        }
    }

    var.info["NS"].push_back(convert(samplesWithData));
    var.info["DP"].push_back(convert(coverage));
    var.info["RO"].push_back(convert(refAlleleObservations));

    var.info["NUMALT"].push_back(convert(altAlleles.size()));

    if (parameters.showReferenceRepeats && !repeats.empty()) {
        stringstream repeatsstr;
        for (map<string, int>::iterator c = repeats.begin(); c != repeats.end(); ++c) {
            repeatsstr << c->first << ":" << c->second << "|";
        }
        string repeatstr = repeatsstr.str();
        TRY { repeatstr = repeatstr.substr(0, repeatstr.size() - 1); } CATCH;
        var.info["REPEAT"].clear();
        var.info["REPEAT"].push_back(repeatstr);
    }

    var.info["ODDS"].push_back(convert(bestComboOddsRatio));

    // samples

    // for ordering GLs
    // ordering is F(j/k) = (k*(k+1)/2)+j.
    map<int, map<string, int> > vcfGenotypeOrder;
    for (map<int, vector<Genotype> >::iterator gtg = genotypesByPloidy.begin(); gtg != genotypesByPloidy.end(); ++gtg) {

        int groupPloidy = gtg->first;
        vector<Genotype>& genotypes = gtg->second;
        for (vector<Genotype>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
            Genotype* genotypePtr = &*g;
            Genotype& genotype = *g;
            string genotypeStr = genotype.str();
            // only provide output for genotypes for which we have data
            bool fullySpecified = true;
            vector<int> gtspec;
            genotype.relativeGenotype(gtspec, refbase, altAlleles);
            // null allele case handled by the fact that we don't have any null alternate alleles
            for (vector<int>::iterator n = gtspec.begin(); n != gtspec.end(); ++n) {
                if (*n < 0) {
                    fullySpecified = false;
                    break;
                }
            }
            if (fullySpecified) {
                // XXX TODO ... EVIL HACKS
                if (groupPloidy == 2) {
                    int j = gtspec.front();
                    int k = gtspec.back();
                    vcfGenotypeOrder[groupPloidy][genotypeStr] = (k * (k + 1) / 2) + j;
                } else if (groupPloidy == 1) {
                    vcfGenotypeOrder[groupPloidy][genotypeStr] = gtspec.front();
                } else {
                    // XXX TODO ...
                }
            }
        }
    }

    // get the best genotypes from the combos, and set the output GTs and GQs using them
    for (vector<string>::iterator sn = sampleNames.begin(); sn != sampleNames.end(); ++sn) {
        string& sampleName = *sn;
        GenotypeComboMap::iterator gc = comboMap.find(sampleName);
        Results::iterator s = find(sampleName);
        map<string, vector<string> >& sampleOutput = var.samples[sampleName];
        if (gc != comboMap.end() && s != end()) {

            Sample& sample = *gc->second->sample;
            Result& sampleLikelihoods = s->second;
            Genotype* genotype = gc->second->genotype;
            sampleOutput["GT"].push_back(genotype->relativeGenotype(refbase, altAlleles));

            if (parameters.calculateMarginals) {
                sampleOutput["GQ"].push_back(convert(float2phred(1 - safe_exp(sampleLikelihoods.front().marginal))));
            }

            sampleOutput["DP"].push_back(convert(sample.observationCount()));
            sampleOutput["RO"].push_back(convert(sample.observationCount(refbase)));
            sampleOutput["QR"].push_back(convert(sample.qualSum(refbase)));


            // get data likelihoods for present genotypes, none if we have excluded genotypes from data likelihood calculations
            map<int, string> genotypeLikelihoods;
            if (!parameters.excludeUnobservedGenotypes && !parameters.excludePartiallyObservedGenotypes) {

                for (Result::iterator g = sampleLikelihoods.begin(); g != sampleLikelihoods.end(); ++g) {
                    if (g->genotype->hasNullAllele()) {
                        // if the genotype has null (unspecified) alleles, find
                        // the fully specified genotypes it can match with.
                        vector<Genotype*> nullmatchgts = g->genotype->nullMatchingGenotypes(genotypesByPloidy[g->genotype->ploidy]);
                        // the gls for these will be the same, so the gl for
                        // this genotype can be used for all of them.  these
                        // are the genotypes which the sample does not have,
                        // but for which one allele or no alleles match
                        for (vector<Genotype*>::iterator n = nullmatchgts.begin(); n != nullmatchgts.end(); ++n) {
                            map<string, int>::iterator o = vcfGenotypeOrder[(*n)->ploidy].find((*n)->str());
                            if (o != vcfGenotypeOrder[(*n)->ploidy].end()) {
                                genotypeLikelihoods[o->second] = convert(ln2log10(g->prob));
                            }
                        }
                    } else {
                        // otherwise, we are well-specified, and only one
                        // genotype should match
                        map<string, int>::iterator o = vcfGenotypeOrder[g->genotype->ploidy].find(g->genotype->str());
                        if (o != vcfGenotypeOrder[g->genotype->ploidy].end()) {
                            genotypeLikelihoods[o->second] = convert(ln2log10(g->prob));
                        }
                    }

                }

                //if (altAlleles.size() == 1) {
                //    if (!genotypeLikelihoods.size() == 3) {
                //        assert(false);
                //    }
                //}

                vector<string>& datalikelihoods = sampleOutput["GL"];
                // output is sorted by map
                for (map<int, string>::iterator gl = genotypeLikelihoods.begin(); gl != genotypeLikelihoods.end(); ++gl) {
                    datalikelihoods.push_back(gl->second);
                }


            }

        }
    }

    return var;

}
