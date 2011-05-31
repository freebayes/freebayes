#include "ResultData.h"
#include "TryCatch.h"

using namespace std;



vcf::Variant& Results::vcf(
        vcf::Variant& var, // variant to update
        long double pHom,
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
                break;
            case ALLELE_MNP:
                if (var.ref.size() < altAllele.length) {
                    var.ref = parser->referenceSubstr(referencePosition, altAllele.length);
                }
                break;
            case ALLELE_DELETION:
                // extend the reference sequence
                if (var.ref.size() - 1 <= altAllele.length) {
                    // XXX this a hack, due to fact that deletions are always intra-base, so this gets the previous base
                    // XXX this means things won't work if we have SNPs mixed with INDEL calls...
                    var.ref = parser->referenceSubstr(referencePosition, altAllele.length + 1);
                }
                break;
            case ALLELE_INSERTION:
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
                altSequence.insert(1, altAllele.alternateSequence.substr(1));
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


    // loop over all alternate alleles
    for (vector<Allele>::iterator aa = altAlleles.begin(); aa != altAlleles.end(); ++aa) {

        Allele& altAllele = *aa;
        string altbase = altAllele.base();

        // count alternate alleles in the best genotyping
        int alternateCount = 0;
        int alleleCount = 0;
        // reference / alternate base counts by strand
        map<string, int> altCountBySample;
        map<string, int> altQualBySample;
        int alternateObsCount = 0;
        // het counts
        int hetReferenceObsCount = 0;
        int hetAlternateObsCount = 0;
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
                }

                alternateCount += genotype->alleleCount(altbase);
                alleleCount += genotype->ploidy;

                if (!genotype->homozygous) {
                    hetReferenceObsCount += sample.observationCount(refbase);
                    hetAlternateObsCount += sample.observationCount(altbase);
                    if (hetAlternateObsCount > 0) {
                        ++hetAltRefSamples;
                    }
                } else {
                    if (genotype->alleleCount(refbase) > 0) {
                        ++homRefSamples;
                    } else {
                        ++homAltSamples;
                    }
                }

                int altCount = sample.observationCount(altbase);
                altCountBySample[*sampleName] = altCount;
                alternateObsCount += altCount;

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

        int hetAllObsCount = hetReferenceObsCount + hetAlternateObsCount;

        unsigned int basesLeft = 0;
        unsigned int basesRight = 0;
        unsigned int readsLeft = 0;
        unsigned int readsRight = 0;
        unsigned int endLeft = 0;
        unsigned int endRight = 0;

        unsigned int mqsum = 0;

        unsigned int properPairs = 0;

        map<string, int> obsBySequencingTechnology;

        vector<Allele*>& alternateAlleles = alleleGroups.at(altbase);
        for (vector<Allele*>::iterator app = alternateAlleles.begin(); app != alternateAlleles.end(); ++app) {
            Allele& allele = **app;
            if (allele.isProperPair) {
                ++properPairs;
            }
            if (!allele.sequencingTechnology.empty()) {
                ++obsBySequencingTechnology[allele.sequencingTechnology];
            }
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

        //string refbase = parser->currentReferenceBase();
        // positional information
        // CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT
        //out.setf(ios::fixed,ios::floatfield);
        //out.precision(5);

        var.info["AC"].push_back(convert(alternateCount));
        var.info["AN"].push_back(convert(alleleCount));
        var.info["AF"].push_back(convert((alleleCount == 0) ? 0 : (double) alternateCount / (double) alleleCount));
        var.info["AA"].push_back(convert(alternateObsCount));
            //<< "ADR=" << ((alleleCount == 0 || alternateCount == 0) ?
            //        0 : ( (double) altAlleleObservations / (double) alternateCount ) / ( (double) coverage / (double) samplesWithData )<< ";"
            //<< "HETAR=" << hetAltRefSamples << ";"
            //<< "HOMA=" << homAltSamples << ";"
            //<< "HOMR=" << homRefSamples << ";"
            //<< "SRF=" << baseCountsForwardTotal.first << ";"
            //<< "SRR=" << baseCountsReverseTotal.first << ";"
            //<< "SAF=" << baseCountsForwardTotal.second << ";"
            //<< "SAR=" << baseCountsReverseTotal.second << ";"
            //<< "SRB=" << ((referenceObsCount == 0) ? 0 : (double) baseCountsForwardTotal.first / (double) referenceObsCount) << ";"
            //<< "SAB=" << ((alternateObsCount == 0) ? 0 : (double) baseCountsForwardTotal.second / (double) alternateObsCount) << ";"
        //var.info["SRP"].push_back(convert((referenceObsCount == 0) ? 0 : ln2phred(hoeffdingln(baseCountsForwardTotal.first, referenceObsCount, 0.5))));
        var.info["SAP"].push_back(convert((alternateObsCount == 0) ? 0 : ln2phred(hoeffdingln(baseCountsForwardTotal.second, alternateObsCount, 0.5))));
            //<< "ABR=" << hetReferenceObsCount <<  ";"
            //<< "ABA=" << hetAlternateObsCount <<  ";"
        var.info["AB"].push_back(convert((hetAllObsCount == 0) ? 0 : (double) hetReferenceObsCount / (double) hetAllObsCount ));
        var.info["ABP"].push_back(convert((hetAllObsCount == 0) ? 0 : ln2phred(hoeffdingln(hetReferenceObsCount, hetAllObsCount, 0.5))));
        var.info["RUN"].push_back(convert(parser->homopolymerRunLeft(altbase) + 1 + parser->homopolymerRunRight(altbase)));
        var.info["MQM"].push_back(convert((alternateAlleles.size() == 0) ? 0 : (double) mqsum / (double) alternateAlleles.size()));
            //<< "RL=" << readsLeft << ";"
            //<< "RR=" << readsRight << ";"
        var.info["RPP"].push_back(convert(ln2phred(hoeffdingln(readsLeft, readsRight + readsLeft, 0.5))));
            //<< "EL=" << endLeft << ";"
            //<< "ER=" << endRight << ";"
        var.info["EPP"].push_back(convert((basesLeft + basesRight == 0) ? 0 : ln2phred(hoeffdingln(endLeft, endLeft + endRight, 0.5))));
            //<< "BL=" << basesLeft << ";"
            //<< "BR=" << basesRight << ";"
            //<< "LRB=" << ((double) max(basesLeft, basesRight) / (double) (basesRight + basesLeft) - 0.5) * 2 << ";"
            //<< "LRBP=" << ((basesLeft + basesRight == 0) ? 0 : ln2phred(hoeffdingln(basesLeft, basesLeft + basesRight, 0.5))) << ";"
        var.info["PAIRED"].push_back(convert((alternateAlleles.size() == 0) ? 0 : (double) properPairs / (double) alternateAlleles.size()));

        for (vector<string>::iterator st = sequencingTechnologies.begin();
                st != sequencingTechnologies.end(); ++st) { string& tech = *st;
            var.info["technology." + tech].push_back(convert((alternateAlleles.size() == 0) ? 0
                        : (double) obsBySequencingTechnology[tech] / (double) alternateAlleles.size() ));
        }

        if (bestOverallComboIsHet) {
            var.infoFlags["BVAR"] = true;
        }

        // allele class
        if (altAllele.type == ALLELE_DELETION) {
            var.infoFlags["DEL"] = true;
            // what is the class of deletion
            // microsatellite repeat?
            // "novel"?
            // how large is the repeat, if there is one?
        } else if (altAllele.type == ALLELE_INSERTION) {
            var.infoFlags["INS"] = true;
        } else if (altAllele.type == ALLELE_SNP) {
            var.infoFlags["SNP"] = true;
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
            var.infoFlags["MNP"] = true;
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
                // output is sorted by map
                for (map<int, string>::iterator gl = genotypeLikelihoods.begin(); gl != genotypeLikelihoods.end(); ++gl) {
                    datalikelihoods.push_back(gl->second);
                }


            }

        }
    }

    return var;

}
