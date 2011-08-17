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
        vector<Allele>& altAlleles,
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
    var.ref = refbase;

    // the reference sequence should be able to encompass all events at the site, +1bp on the left
    for (vector<Allele>::iterator aa = altAlleles.begin(); aa != altAlleles.end(); ++aa) {

        Allele& altAllele = *aa;

        switch (altAllele.type) {
            case ALLELE_SNP:
            case ALLELE_REFERENCE:
                break;
            case ALLELE_MNP:
                if (var.ref.size() < altAllele.length) {
                    var.ref = parser->referenceSubstr(referencePosition, altAllele.length);
                }
                break;
            case ALLELE_DELETION:
                // extend the reference sequence
                if (var.ref.size() - 1 <= altAllele.length) {
                    var.ref = parser->referenceSubstr(referencePosition, altAllele.length + 1);
                }
                break;
            case ALLELE_INSERTION:
                break;
            case ALLELE_COMPLEX:
                // if the complex event involves the deletion of some bp, then
                // we have to treat it like a deletion during reporting,
                // otherwise, as insertion
                if (var.ref.size() < altAllele.referenceLength + 1) {
                    var.ref = parser->referenceSubstr(referencePosition, altAllele.referenceLength + 1);
                }
                break;
            default:
                cerr << "Unhandled allele type: " << altAllele.typeStr() << endl;
                break;
        }

    }

    for (vector<Allele>::iterator aa = altAlleles.begin(); aa != altAlleles.end(); ++aa) {

        Allele& altAllele = *aa;
        string altSequence;

        switch (altAllele.type) {
            case ALLELE_REFERENCE:
                altSequence = var.ref;
                break;
            case ALLELE_SNP:
                altSequence = var.ref;
                // XXX hack, what happens when we are co-present with indels?
                // ... this would break.  in the current arrangement this is not an issue
                altSequence.replace(0, 1, altAllele.base());
                break;
            case ALLELE_MNP:
                // XXX also a hack, implicitly groups MNPs with SNPs
                altSequence = var.ref;
                altSequence.replace(0, altAllele.length, altAllele.base());
                break;
            case ALLELE_DELETION:
                altSequence = var.ref;
                altSequence.erase(1, altAllele.length);
                break;
            case ALLELE_INSERTION:
                altSequence = var.ref;
                // XXX hack... resolves the alt allele semantics issue
                // we have to strip off the leading 'I' from the allele name
                altSequence.insert(1, altAllele.base().substr(1));
                break;
            case ALLELE_COMPLEX:
                altSequence = var.ref;
                altSequence.erase(1, altAllele.referenceLength);
                altSequence.insert(1, altAllele.base().substr(1)); // chop off leading "C"
                break;
            default:
                cerr << "Unhandled allele type: " << altAllele.typeStr() << endl;
                break;
        }

        var.alt.push_back(altSequence);
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
    var.format.push_back("RA");
    var.format.push_back("QR");
    var.format.push_back("AA");
    var.format.push_back("QA");
    var.format.push_back("GL");
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
        unsigned int hetAlternateObsCount = 0;
        unsigned int homReferenceRefObsCount = 0;
        unsigned int homReferenceAltObsCount = 0;
        unsigned int homAlternateRefObsCount = 0;
        unsigned int homAlternateAltObsCount = 0;
        unsigned int hetAltRefSamples = 0;
        unsigned int homAltSamples = 0;
        unsigned int homRefSamples = 0;
        unsigned int refSampleObsCount = 0; // depth in hom-ref samples
        unsigned int altSampleObsCount = 0; // depth in samples with called alternates

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

                if (!genotype->homozygous) {
                    // het case
                    hetReferenceObsCount += sample.observationCount(refbase);
                    hetAlternateObsCount += sample.observationCount(altbase);
                    if (hetAlternateObsCount > 0) {
                        ++hetAltRefSamples;
                        altSampleObsCount += sample.observationCount();
                    }
                } else {
                    // homozygous cases
                    if (genotype->alleleCount(refbase) > 0) {
                        homReferenceRefObsCount += sample.observationCount(refbase);
                        homReferenceAltObsCount += sample.observationCount(altbase);
                        ++homRefSamples;
                        refSampleObsCount += sample.observationCount();
                    } else if (genotype->alleleCount(altbase) > 0) {
                        homAlternateAltObsCount += sample.observationCount(altbase);
                        homAlternateRefObsCount += sample.observationCount(refbase); // disagreeing obs
                        ++homAltSamples;
                        altSampleObsCount += sample.observationCount();
                    } // specifically exclude other alts from these sums
                }

                unsigned int altCount = sample.observationCount(altbase);
                altCountBySample[*sampleName] = altCount;
                //altObsCount += altCount;

                altQualBySample[*sampleName] = sample.qualSum(altbase);

                // TODO cleanup
                pair<pair<int,int>, pair<int,int> > baseCounts = sample.baseCount(refbase, altbase);
                baseCountsForwardBySample[*sampleName] = baseCounts.first;
                baseCountsReverseBySample[*sampleName] = baseCounts.second;
                baseCountsForwardTotal.first += baseCounts.first.first;
                baseCountsForwardTotal.second += baseCounts.first.second;
                baseCountsReverseTotal.first += baseCounts.second.first;
                baseCountsReverseTotal.second += baseCounts.second.second;
                // TODO, redundant...
                //refAlleleObservations += baseCounts.first.first + baseCounts.second.first;
                //altAlleleObservations += baseCounts.first.second + baseCounts.second.second;
            }
        }

        unsigned int hetAllObsCount = hetReferenceObsCount + hetAlternateObsCount;

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
        }

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
        var.info["AA"].push_back(convert(altObsCount));
        if (homRefSamples > 0 && hetAltRefSamples + homAltSamples > 0) {
            double altSampleAverageDepth = (double) altSampleObsCount
                       / ( (double) hetAltRefSamples + (double) homAltSamples );
            double refSampleAverageDepth = (double) refSampleObsCount / (double) homRefSamples;
            var.info["DPRA"].push_back(convert(altSampleAverageDepth / refSampleAverageDepth));
        } else {
            var.info["DPRA"].push_back(convert(0));
        }
            //<< "HETAR=" << hetAltRefSamples << ";"
            //<< "HOMA=" << homAltSamples << ";"
            //<< "HOMR=" << homRefSamples << ";"
        var.info["SRP"].clear(); // XXX hack
        var.info["SRP"].push_back(convert((refObsCount == 0) ? 0 : ln2phred(hoeffdingln(baseCountsForwardTotal.first, refObsCount, 0.5))));
        var.info["SAP"].push_back(convert((altObsCount == 0) ? 0 : ln2phred(hoeffdingln(baseCountsForwardTotal.second, altObsCount, 0.5))));
            //<< "ABR=" << hetReferenceObsCount <<  ";"
            //<< "ABA=" << hetAlternateObsCount <<  ";"
        var.info["AB"].push_back(convert((hetAllObsCount == 0) ? 0 : (double) hetReferenceObsCount / (double) hetAllObsCount ));
        var.info["ABP"].push_back(convert((hetAllObsCount == 0) ? 0 : ln2phred(hoeffdingln(hetReferenceObsCount, hetAllObsCount, 0.5))));
        var.info["RUN"].push_back(convert(parser->homopolymerRunLeft(altbase) + 1 + parser->homopolymerRunRight(altbase)));
        var.info["MQM"].push_back(convert((altObsCount == 0) ? 0 : (double) altmqsum / (double) altObsCount));
        var.info["RPP"].push_back(convert((altObsCount == 0) ? 0 : ln2phred(hoeffdingln(altReadsLeft, altReadsRight + altReadsLeft, 0.5))));
        var.info["EPP"].push_back(convert((altBasesLeft + altBasesRight == 0) ? 0 : ln2phred(hoeffdingln(altEndLeft, altEndLeft + altEndRight, 0.5))));
        var.info["PAIRED"].push_back(convert((altObsCount == 0) ? 0 : (double) altproperPairs / (double) altObsCount));

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
            // ts/tv
            if (isTransition(refbase, altbase)) {
                var.infoFlags["TS"] = true;
            } else {
                var.infoFlags["TV"] = true;
            }

            // CpG
            if (parser->isCpG(altbase)) {
                var.infoFlags["CpG"] = true;
            }
        } else if (altAllele.type == ALLELE_MNP) {
            var.info["TYPE"].push_back("mnp");
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

                sampleOutput["AA"].push_back(convert(altCountBySample[sampleName]));
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
    var.info["RA"].push_back(convert(refAlleleObservations));

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
    map<int, map<Genotype*, int> > vcfGenotypeOrder;
    for (map<int, vector<Genotype> >::iterator gtg = genotypesByPloidy.begin(); gtg != genotypesByPloidy.end(); ++gtg) {
        int groupPloidy = gtg->first;
        vector<Genotype>& genotypes = gtg->second;
        for (vector<Genotype>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
            Genotype* genotypePtr = &*g;
            Genotype& genotype = *g;
            vector<int> gtspec;
            genotype.relativeGenotype(gtspec, refbase, altAlleles);
            // XXX TODO ... EVIL HACKS
            if (groupPloidy == 2) {
                int j = gtspec.front();
                int k = gtspec.back();
                vcfGenotypeOrder[groupPloidy][genotypePtr] = (k * (k + 1) / 2) + j;
            } else if (groupPloidy == 1) {
                vcfGenotypeOrder[groupPloidy][genotypePtr] = gtspec.front();
            } else {
                // XXX TODO ...
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
            sampleOutput["RA"].push_back(convert(sample.observationCount(refbase)));
            sampleOutput["QR"].push_back(convert(sample.qualSum(refbase)));


            // get data likelihoods for present genotypes, none if we have excluded genotypes from data likelihood calculations
            map<int, string> genotypeLikelihoods;
            if (!parameters.excludeUnobservedGenotypes && !parameters.excludePartiallyObservedGenotypes) {

                // TODO
                //vector<string>& taggedDataLikelihoods = sampleOutput["GLE"];

                for (Result::iterator g = sampleLikelihoods.begin(); g != sampleLikelihoods.end(); ++g) {
                    //vector<Genotype>& genotypes = genotypesByPloidy[sample.ploidy];
                    genotypeLikelihoods[vcfGenotypeOrder[g->genotype->ploidy][g->genotype]] = convert(ln2log10(g->prob));

                    //string gle = g->genotype->relativeGenotype(refbase, altAlleles) + "=" + convert(ln2log10(g->prob));
                    //taggedDataLikelihoods.push_back(gle);
                }

                vector<string>& datalikelihoods = sampleOutput["GL"];
                //vector<string>& datalikelihoodsExplicit = sampleOutput["GLE"];
                // output is sorted by map
                for (map<int, string>::iterator gl = genotypeLikelihoods.begin(); gl != genotypeLikelihoods.end(); ++gl) {
                    datalikelihoods.push_back(gl->second);
                    //datalikelihoodsExplicit.push_back(gl->second);
                }


            }

        }
    }

    return var;

}
