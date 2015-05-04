
variant information is stored in:

    [line 730:]
    vcf::Variant var;

... the fields are set and printed to an output stream:

    out << results.vcf(
                var,
                pHom,
                bestComboOddsRatio,
                samples,
                referenceBase,
                alts,
                repeats,
                genotypingTotalIterations,
                parser->sampleList,
                coverage,
                bestCombo,
                alleleGroups,
                partialObservationGroups,
                partialObservationSupport,
                genotypesByPloidy,
                parser->sequencingTechnologies,
                parser)
                << endl;

GL order issue
==============

vcflib/src/Variant.cpp

l. 2054 

    list<list<int> > glorder(int ploidy, int alts) {}

l. 2028

    list<list<int> > _glorder(int ploidy, int alts) {}
    

