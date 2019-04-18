#include "ResultData.h"
#include "TryCatch.h"

using namespace std;



vcflib::Variant& Results::vcf(
    vcflib::Variant& var, // variant to update
    BigFloat pHom,
    long double bestComboOddsRatio,
    //long double alleleSamplingProb,
    Samples& samples,
    string refbase,
    vector<Allele>& altAllelesIncludingNulls,
    map<string, int> repeats,
	int genotypingIterations,
    vector<string>& sampleNames,
    int coverage,
    GenotypeCombo& genotypeCombo,
    map<string, vector<Allele*> >& alleleGroups,
    map<string, vector<Allele*> >& partialObservationGroups,
    map<Allele*, set<Allele*> >& partialObservationSupport,
    map<int, vector<Genotype> >& genotypesByPloidy,
    vector<string>& sequencingTechnologies,
    AlleleParser* parser) {

    Parameters& parameters = parser->parameters;

    GenotypeComboMap comboMap;
    genotypeCombo2Map(genotypeCombo, comboMap);

    // set up the reported reference allele
    long int referencePosition = (long int) parser->currentPosition; // 0-based

    // remove NULL alt alleles
    vector<Allele> altAlleles;
    for (vector<Allele>::iterator aa = altAllelesIncludingNulls.begin(); aa != altAllelesIncludingNulls.end(); ++aa) {
        if (!aa->isNull()) {
            altAlleles.push_back(*aa);
        }
    }

    map<string, string> adjustedCigar;
    vector<Allele>& adjustedAltAlleles = altAlleles; // just an alias
    for (vector<Allele>::iterator aa = altAlleles.begin(); aa != altAlleles.end(); ++aa) {
        adjustedCigar[aa->base()] = aa->cigar;
        var.alt.push_back(aa->alternateSequence);
    }

    var.ref = refbase;
    assert(!var.ref.empty());

    // get the required size of the reference sequence
    // strip identical bases from start and/or end of alleles
    // if bases have been stripped from the beginning,

    // set up VCF record-wide variables

    var.sequenceName = parser->currentSequenceName;
    var.position = referencePosition + 1;
    var.id = ".";
    var.filter = ".";

    // note that we set QUAL to 0 at loci with no data
    var.quality = max((long double) 0, nan2zero(big2phred(pHom)));
    if (coverage == 0) {
        var.quality = 0;
    }


    // set up format string

    var.format.clear();
    var.format.push_back("GT");
    if (parameters.calculateMarginals) var.format.push_back("GQ");
    // XXX
    var.format.push_back("DP");
    var.format.push_back("AD");
    var.format.push_back("RO");
    var.format.push_back("QR");
    var.format.push_back("AO");
    var.format.push_back("QA");
    // add GL/GLE later, when we know if we need to use one or the other

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
    long double refReadSoftClipSum = 0;
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

    //var.info["XRM"].push_back(convert(refReadMismatchRate));
    //var.info["XRS"].push_back(convert(refReadSNPRate));
    //var.info["XRI"].push_back(convert(refReadIndelRate));

    var.info["MQMR"].push_back(convert((refObsCount == 0) ? 0 : (double) refmqsum / (double) refObsCount));
    var.info["RPPR"].push_back(convert((refObsCount == 0) ? 0 : nan2zero(ln2phred(hoeffdingln(refReadsLeft, refReadsRight + refReadsLeft, 0.5)))));
    var.info["EPPR"].push_back(convert((refBasesLeft + refBasesRight == 0) ? 0 : nan2zero(ln2phred(hoeffdingln(refEndLeft, refEndLeft + refEndRight, 0.5)))));
    var.info["PAIREDR"].push_back(convert((refObsCount == 0) ? 0 : (double) refProperPairs / (double) refObsCount));

    //var.info["HWE"].push_back(convert(nan2zero(ln2phred(genotypeCombo.hweComboProb()))));
    var.info["GTI"].push_back(convert(genotypingIterations));

    // loop over all alternate alleles
    for (vector<Allele>::iterator aa = altAlleles.begin(); aa != altAlleles.end(); ++aa) {

        Allele& altAllele = *aa;
        string altbase = altAllele.base();

        // count alternate alleles in the best genotyping
        unsigned int alternateCount = 0;
        unsigned int alleleCount = 0;
        double alternateQualitySum = 0;
        double partialObservationCount = 0;
        double partialObservationQualitySum;
        // reference / alternate base counts by strand
        //map<string, unsigned int> altCountBySample;
        //map<string, unsigned int> altQualBySample;
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

        StrandBaseCounts baseCountsTotal;
        map<string, StrandBaseCounts> baseCountsBySample;
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
                //altCountBySample[*sampleName] = altCount;

                //altQualBySample[*sampleName] = sample.qualSum(altbase);

                StrandBaseCounts baseCounts = sample.strandBaseCount(refbase, altbase);
                baseCountsBySample[*sampleName] = baseCounts;
                baseCountsTotal.forwardRef += baseCounts.forwardRef;
                baseCountsTotal.forwardAlt += baseCounts.forwardAlt;
                baseCountsTotal.reverseRef += baseCounts.reverseRef;
                baseCountsTotal.reverseAlt += baseCounts.reverseAlt;
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

        // TODO we need a partial obs structure to annotate partial obs
        map<string, vector<Allele*> >::iterator f = alleleGroups.find(altbase);
        if (f != alleleGroups.end()) {
            vector<Allele*>& alternateAlleles = alleleGroups.at(altbase);
            // TODO XXX XXX adjust to use partial observations
            altObsCount = alternateAlleles.size();
            for (vector<Allele*>::iterator app = alternateAlleles.begin(); app != alternateAlleles.end(); ++app) {
                Allele& allele = **app;
                altReadMismatchSum += allele.readMismatchRate;
                altReadSNPSum += allele.readSNPRate;
                altReadIndelSum += allele.readIndelRate;
				// TODO: add altReadSoftClipRate (avg)
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

        //var.info["XAM"].push_back(convert(altReadMismatchRate));
        //var.info["XAS"].push_back(convert(altReadSNPRate));
        //var.info["XAI"].push_back(convert(altReadIndelRate));

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
        var.info["PAO"].push_back(convert(samples.partialObservationCount(altbase)));
        var.info["QA"].push_back(convert(samples.qualSum(altbase)));
        var.info["PQA"].push_back(convert(samples.partialQualSum(altbase)));
        if (homRefSamples > 0 && hetAltSamples + homAltSamples > 0) {
            double altSampleAverageDepth = (double) altSampleObsCount
                / ( (double) hetAltSamples + (double) homAltSamples );
            double refSampleAverageDepth = (double) refSampleObsCount / (double) homRefSamples;
            var.info["DPRA"].push_back(convert(altSampleAverageDepth / refSampleAverageDepth));
        } else {
            var.info["DPRA"].push_back(convert(0));
        }

        var.info["SRP"].clear(); // XXX hack
        var.info["SRF"].clear();
        var.info["SRR"].clear();
        var.info["SRF"].push_back(convert(baseCountsTotal.forwardRef));
        var.info["SRR"].push_back(convert(baseCountsTotal.reverseRef));
        var.info["SRP"].push_back(convert((refObsCount == 0) ? 0 : nan2zero(ln2phred(hoeffdingln(baseCountsTotal.forwardRef, refObsCount, 0.5)))));
        var.info["SAF"].push_back(convert(baseCountsTotal.forwardAlt));
        var.info["SAR"].push_back(convert(baseCountsTotal.reverseAlt));
        var.info["SAP"].push_back(convert((altObsCount == 0) ? 0 : nan2zero(ln2phred(hoeffdingln(baseCountsTotal.forwardAlt, altObsCount, 0.5)))));
        var.info["AB"].push_back(convert((hetAllObsCount == 0) ? 0 : nan2zero((double) hetAlternateObsCount / (double) hetAllObsCount )));
        var.info["ABP"].push_back(convert((hetAllObsCount == 0) ? 0 : nan2zero(ln2phred(hoeffdingln(hetAlternateObsCount, hetAllObsCount, 0.5)))));
        var.info["RUN"].push_back(convert(parser->homopolymerRunLeft(altbase) + 1 + parser->homopolymerRunRight(altbase)));
        var.info["MQM"].push_back(convert((altObsCount == 0) ? 0 : nan2zero((double) altmqsum / (double) altObsCount)));
        var.info["RPP"].push_back(convert((altObsCount == 0) ? 0 : nan2zero(ln2phred(hoeffdingln(altReadsLeft, altReadsRight + altReadsLeft, 0.5)))));
        var.info["RPR"].push_back(convert(altReadsRight));
		var.info["RPL"].push_back(convert(altReadsLeft));
        var.info["EPP"].push_back(convert((altBasesLeft + altBasesRight == 0) ? 0 : nan2zero(ln2phred(hoeffdingln(altEndLeft, altEndLeft + altEndRight, 0.5)))));
        var.info["PAIRED"].push_back(convert((altObsCount == 0) ? 0 : nan2zero((double) altproperPairs / (double) altObsCount)));
        var.info["CIGAR"].push_back(adjustedCigar[altAllele.base()]);
        var.info["MEANALT"].push_back(convert((hetAltSamples + homAltSamples == 0) ? 0 : nan2zero((double) uniqueAllelesInAltSamples / (double) (hetAltSamples + homAltSamples))));

        for (vector<string>::iterator st = sequencingTechnologies.begin();
             st != sequencingTechnologies.end(); ++st) { string& tech = *st;
            var.info["technology." + tech].push_back(convert((altObsCount == 0) ? 0
                                                             : nan2zero((double) altObsBySequencingTechnology[tech] / (double) altObsCount )));
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

            /*
            // CpG
            if (parser->isCpG(altbase)) {
                var.infoFlags["CpG"] = true;
            }
            */
        } else if (altAllele.type == ALLELE_MNP) {
            var.info["TYPE"].push_back("mnp");
        } else {
            /*
            cerr << "What is this?"
                 << "type: " << altAllele.type << " "
                 << "allele: " << altAllele << endl;
            */
        }
        var.info["LEN"].push_back(convert(altAllele.length));

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
            //refAlleleObservations += sample.observationCount(refbase);
            refAlleleObservations += sample.observationCount(refbase);

            ++samplesWithData;
        }
    }

    var.info["NS"].push_back(convert(samplesWithData));
    var.info["DP"].push_back(convert(coverage));
    var.info["RO"].push_back(convert(refAlleleObservations));
    var.info["PRO"].push_back(convert(samples.partialObservationCount(refbase)));
    var.info["QR"].push_back(convert(samples.qualSum(refbase)));
    var.info["PQR"].push_back(convert(samples.partialQualSum(refbase)));

    // tally partial observations to get a mean coverage per bp of reference
    int haplotypeLength = refbase.size();
    int basesInObservations = 0;

    for (map<string, vector<Allele*> >::iterator g = alleleGroups.begin(); g != alleleGroups.end(); ++g) {
        for (vector<Allele*>::iterator a = g->second.begin(); a != g->second.end(); ++a) {
            basesInObservations += (*a)->alternateSequence.size();
        }
    }

    for (map<Allele*, set<Allele*> >::iterator p = partialObservationSupport.begin(); p != partialObservationSupport.end(); ++p) {
        basesInObservations += p->first->alternateSequence.size();
    }

    double depthPerBase = (double) basesInObservations / (double) haplotypeLength;
    var.info["DPB"].push_back(convert(depthPerBase));

    // number of alternate alleles
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

    bool outputExplicitGenotypeLikelihoods = false;
    bool outputAnyGenotypeLikelihoods = true;

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
                if (groupPloidy == 2) {
                    int j = gtspec.front();
                    int k = gtspec.back();
                    vcfGenotypeOrder[groupPloidy][genotypeStr] = (k * (k + 1) / 2) + j;
                } else if (groupPloidy == 1) {
                    vcfGenotypeOrder[groupPloidy][genotypeStr] = gtspec.front();
                } else {
                    outputAnyGenotypeLikelihoods = false; // XXX prevents output of GLs for polyploid data
                    outputExplicitGenotypeLikelihoods = true;
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
            if (sample.observationCount() == 0) {
                continue;
            }

            sampleOutput["GT"].push_back(genotype->relativeGenotype(refbase, altAlleles));

            if (parameters.calculateMarginals) {
                double val = nan2zero(big2phred((BigFloat)1 - big_exp(sampleLikelihoods.front().marginal)));
                if (parameters.strictVCF)
                    sampleOutput["GQ"].push_back(convert(int(round(val))));
                else
                    sampleOutput["GQ"].push_back(convert(val));
            }

            sampleOutput["DP"].push_back(convert(sample.observationCount()));
            sampleOutput["AD"].push_back(convert(sample.observationCount(refbase)));
            sampleOutput["RO"].push_back(convert(sample.observationCount(refbase)));
            sampleOutput["QR"].push_back(convert(sample.qualSum(refbase)));

            for (vector<Allele>::iterator aa = altAlleles.begin(); aa != altAlleles.end(); ++aa) {
                Allele& altAllele = *aa;
                string altbase = altAllele.base();
                sampleOutput["AO"].push_back(convert(sample.observationCount(altbase)));
                sampleOutput["AD"].push_back(convert(sample.observationCount(altbase)));
                sampleOutput["QA"].push_back(convert(sample.qualSum(altbase)));
            }

            if (outputAnyGenotypeLikelihoods && !parameters.excludeUnobservedGenotypes && !parameters.excludePartiallyObservedGenotypes) {

                // get data likelihoods for present genotypes, none if we have excluded genotypes from data likelihood calculations
                if (outputExplicitGenotypeLikelihoods) {

                    if (var.format.back() != "GLE") {
                        var.format.push_back("GLE");
                    }

                    map<string, string> genotypeLikelihoodsExplicit;

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
                                genotypeLikelihoodsExplicit[(*n)->relativeGenotype(refbase, altAlleles)] = convert(ln2log10(g->prob));
                            }
                        } else {
                            // otherwise, we are well-specified, and only one
                            // genotype should match
                            genotypeLikelihoodsExplicit[g->genotype->relativeGenotype(refbase, altAlleles)] = convert(ln2log10(g->prob));
                        }
                    }

                    vector<string> datalikelihoods;
                    for (map<string, string>::iterator gle = genotypeLikelihoodsExplicit.begin(); gle != genotypeLikelihoodsExplicit.end(); ++gle) {
                        datalikelihoods.push_back(gle->first + "^" + gle->second);
                    }
                    sampleOutput["GLE"].push_back(join(datalikelihoods, "|"));

                } else {

                    if (var.format.back() != "GL") {
                        var.format.push_back("GL");
                    }

                    map<int, double> genotypeLikelihoods;
                    map<int, string> genotypeLikelihoodsOutput;

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
                                    genotypeLikelihoods[o->second] = ln2log10(g->prob);
                                }
                            }
                        } else {
                            // otherwise, we are well-specified, and only one
                            // genotype should match
                            map<string, int>::iterator o = vcfGenotypeOrder[g->genotype->ploidy].find(g->genotype->str());
                            if (o != vcfGenotypeOrder[g->genotype->ploidy].end()) {
                                genotypeLikelihoods[o->second] = ln2log10(g->prob);
                            }
                        }
                    }

                    // normalize GLs to 0 max using division by max
                    long double minGL = 0;
                    for (map<int, double>::iterator g = genotypeLikelihoods.begin(); g != genotypeLikelihoods.end(); ++g) {
                        if (g->second < minGL) minGL = g->second;
                    }
                    long double maxGL = minGL;
                    for (map<int, double>::iterator g = genotypeLikelihoods.begin(); g != genotypeLikelihoods.end(); ++g) {
                        if (g->second > maxGL) maxGL = g->second;
                    }

                    if (parameters.limitGL == 0) {
                        for (map<int, double>::iterator g = genotypeLikelihoods.begin(); g != genotypeLikelihoods.end(); ++g) {
                            genotypeLikelihoodsOutput[g->first] = convert(g->second-maxGL);
                        }
                    } else {
                        for (map<int, double>::iterator g = genotypeLikelihoods.begin(); g != genotypeLikelihoods.end(); ++g) {
                            genotypeLikelihoodsOutput[g->first] = convert( max((long double) + parameters.limitGL, (g->second-maxGL)) );
                        }
                    }

                    vector<string>& datalikelihoods = sampleOutput["GL"];
                    // output is sorted by map
                    for (map<int, string>::iterator gl = genotypeLikelihoodsOutput.begin(); gl != genotypeLikelihoodsOutput.end(); ++gl) {
                        datalikelihoods.push_back(gl->second);
                    }

                }
            }
        }
    }

    return var;

}


vcflib::Variant& Results::gvcf(
    vcflib::Variant& var,
    NonCalls& nonCalls,
    AlleleParser* parser) {

    // what is the first position in the nonCalls?
    pair<string, long> start = nonCalls.firstPos();
    const string& startChrom = start.first;
    long startPos = start.second;
    // startPos and endPos are zero-based, half-open -- [startPos,endPos)

    // what is the current position? nb: can't be on a different chrom
    long endPos;
    if (startChrom != parser->currentSequenceName) {
        endPos = parser->reference.sequenceLength(startChrom);
    } else {
        endPos = parser->currentPosition;
    }
    long numSites = endPos - startPos;
    assert(numSites > 0);

    // set up site call
    var.ref = parser->referenceSubstr(startPos, 1);
    var.alt.push_back("<*>");
    var.sequenceName = parser->currentSequenceName;
    var.position = startPos + 1; // output text field is one-based
    var.id = ".";
    var.filter = ".";
    // TODO change to actual quality
    var.quality = 0;
    // set up format string
    var.format.clear();
    var.format.push_back("GQ");
    var.format.push_back("DP");
    var.format.push_back("MIN_DP");
    var.format.push_back("QR");
    var.format.push_back("QA");

    NonCall total = nonCalls.aggregateAll();

    /* This resets min depth to zero if nonCalls is less than numSites. */

    int minDepth = total.minDepth;

    if(numSites != total.nCount){
        minDepth = 0;
    }

    var.info["DP"].push_back(convert((total.refCount+total.altCount) / numSites));
    var.info["MIN_DP"].push_back(convert(minDepth));
    // The text END field is one-based, inclusive. We proudly conflate this
    // with our zero-based, exclusive endPos.
    var.info["END"].push_back(convert(endPos));

    // genotype quality is 1- p(polymorphic)

    map<string, NonCall> perSample;
    nonCalls.aggregatePerSample(perSample);

    // iterate through the samples and aggregate information about them
    for (vector<string>::const_iterator s = parser->sampleList.begin();
         s != parser->sampleList.end(); ++s) {
        const string& sampleName = *s;
        const NonCall& nc = perSample[sampleName];
        map<string, vector<string> >& sampleOutput = var.samples[sampleName];
        long double qual = nc.reflnQ - nc.altlnQ;
        sampleOutput["GQ"].push_back(convert(ln2phred(qual)));


      /* This resets min depth to zero if nonCalls is less than numSites. */

        int minDepth = nc.minDepth;

        if(numSites != nc.nCount){
            minDepth = 0;
        }

        sampleOutput["DP"].push_back(convert((nc.refCount+nc.altCount) / numSites));
        sampleOutput["MIN_DP"].push_back(convert(minDepth));
        sampleOutput["QR"].push_back(convert(ln2phred(nc.reflnQ)));
        sampleOutput["QA"].push_back(convert(ln2phred(nc.altlnQ)));
    }

    return var;
}
