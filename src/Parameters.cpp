#include "Parameters.h"
// XXX A BUG in tclap 1.1.0 means we can only include this in this source file
// otherwise we'll get multiple definition errors
// see http://blog.csdn.net/happyashley/archive/2009/08/28/4492589.aspx (chinese)
#include <tclap/CmdLine.h>

using namespace TCLAP; 
using namespace std;


class MyOutput : public TCLAP::StdOutput {
public:

    string ProgramName;
    string ProgramDescription;
    string ProgramVersion;
    string ProgramDate;
    string ProgramDeveloper;
    string ProgramInstitution;
    string ProgramCopyrightDates;

    vector<ArgStruct> ArgList;

    virtual void failure(TCLAP::CmdLineInterface& c, TCLAP::ArgException& e);
    virtual void usage(TCLAP::CmdLineInterface& c);
    virtual void version(TCLAP::CmdLineInterface& c);

    MyOutput(void);
};


MyOutput::MyOutput(void)
    : ProgramName("freebayes")
    , ProgramDescription("Bayesian SNP and short-INDEL polymorphism discovery program.")
    , ProgramVersion(GIT_HEAD_COMMIT_ID)
    , ProgramDate(COMPILE_DATE)
    , ProgramDeveloper("Gabor T. Marth, Erik Garrison")
    , ProgramInstitution("Boston College")
    , ProgramCopyrightDates("2007, 2008, 2009, 2010.")
{ }

void MyOutput::failure(CmdLineInterface& c, ArgException& e)
{
    cerr << "################################################################################" << endl;
    cerr << "### Program: " << ProgramName << endl;
    cerr << "### Version: " <<  ProgramVersion << endl;
    cerr << "### Release date: " << ProgramDate << endl;
    cerr << "### " << ProgramDescription << endl;
    cerr << "### " << "Copyright: " << ProgramCopyrightDates << endl;
    cerr << "### " << ProgramDeveloper << ", " << ProgramInstitution << endl;
    cerr << "### All rights reserved." << endl;
    cerr << "###" << endl;
    cerr << "### Command line error:" << e.what() << endl;
    cerr << "### For usage please type: " << c.getProgramName() << " --help" << endl;
#ifdef VERBOSE_DEBUG
    cerr << "### (compiled with --debug2 support)" << endl;
#endif
    cerr << "################################################################################" << endl;
}

void MyOutput::usage(CmdLineInterface& c)
{

    cout << "################################################################################" << endl;
    cout << "### Program: " << ProgramName << " Version: " <<  ProgramVersion << " Release date: " << ProgramDate << endl;
    cout << "### " << ProgramDescription << endl;
    cout << "### " << "Copyright: " << ProgramCopyrightDates << " " << ProgramDeveloper << ", " << ProgramInstitution << endl;
    cout << "### All rights reserved." << endl;
    cout << "###" << endl;
    cout << "### Usage: " << c.getProgramName() << " [arguments], where:" << endl;

    for(vector<ArgStruct>::const_iterator it = ArgList.begin(); 
        it != ArgList.end(); it++) {
        ArgStruct arg = *it;

        string idString = "";
        if (arg.longId != "") {
            idString += "--" + arg.longId;
        }
        if (arg.shortId != "") {
            idString += ", -" + arg.shortId;
        }

        string multiString = "single-valued";
        if (arg.multi) {
            multiString = "multi-valued";
        }
        if (arg.required) {
            cout << "### " << idString << " [" << arg.type << ", required, no default, " << multiString << "]" << endl;
        } else {
            cout << "### " << idString << " [" << arg.type << ", optional, default=" << arg.defaultValueString << ", " << multiString << "]" << endl;
        }
        if (arg.constraint.size() > 0) {
            cout << "###     Permitted values: (";
            bool first = true;
            for (vector<string>::const_iterator iter = arg.constraint.begin();
             iter != arg.constraint.end(); iter++) {
                string value = *iter;
                if (! first) {
                    cout << "|";
                } 
                first = false;
                cout << value;
            }
            cout << ")" << endl;
        }
        cout << "###     Description: " << arg.description << endl;
    }
    cout << "################################################################################" << endl;
}

void MyOutput::version(CmdLineInterface& c)
{
    cerr << "################################################################################" << endl;
    cerr << "### Program: " << ProgramName << " Version: " <<  ProgramVersion << " Release date: " << ProgramDate << endl;
    cout << "### " << ProgramDescription << endl;
    cout << "### " << "Copyright: " << ProgramCopyrightDates << " " << ProgramDeveloper << ", " << ProgramInstitution << endl;
    cout << "### All rights reserved." << endl;
    cout << "###" << endl;
    cerr << "################################################################################" << endl;
}


// XXX this is out of date
ostream &operator<<(ostream &out, const Parameters &p) {

    out << "Complete list of parameter values:" << endl
         << "  --bam = " << p.bam << endl
         << "  --fasta = " << p.fasta << endl
         << "  --targets = " << p.targets << endl
         << "  --samples = " << p.samples << endl
         << "  --log = " << p.log << endl
         << "  --useRefAllele = " <<  (p.useRefAllele ? "true" : "false") << endl
         << "  --forceRefAllele = " <<  (p.forceRefAllele ? "true" : "false") << endl
         << "  --MQR = " << p.MQR << endl
         << "  --BQR = " << p.BQR << endl
         << "  --ploidy = " << p.ploidy << endl
         << "  --sampleNaming = " << p.sampleNaming << endl
         << "  --sampleDel = " << p.sampleDel << endl
         << "  --BQL0 = " << p.BQL0 << endl
         << "  --MQL0 = " << p.MQL0 << endl
         << "  --BQL1 = " << p.BQL1 << endl
         << "  --MQL1 = " << p.MQL1 << endl
         << "  --BQL2 = " << p.BQL2 << endl
         << "  --RMU = " << p.RMU << endl
         << "  --IDW = " << p.IDW << endl
         << "  --TH = " << p.TH << endl
         << "  --PVL = " << p.PVL << endl
         << "  --algorithm = " << p.algorithm << endl
         << "  --RDF = " << p.RDF << endl
         << "  --WB = " << p.WB << endl
         << "  --TB = " << p.TB << endl
         << "  --includeMonoB = " <<  (p.includeMonoB ? "true" : "false") << endl
         << "  --TR = " << p.TR << endl
         << "  --I = " << p.I << endl
         << "  --record = " << (p.record ? "true" : "false" ) << endl
         << "  --debug = " <<  (p.debug ? "true" : "false") << endl
         << "  --debug2 = " <<  (p.debug2 ? "true" : "false") << endl
         << "  --outputAlleles = " << p.outputAlleles << endl
         << "  --suppressOutput = " << p.suppressOutput << endl
         << endl;

    return out;

}

Parameters::Parameters (int argc, char** argv) {

    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------
    // constants
    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------

    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------
    // command line
    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------

    //----------------------------------------------------------------------------
    // Create new CmdLine object
    //----------------------------------------------------------------------------

    MyOutput my;
    CmdLine cmd("", ' ', my.ProgramVersion);
    cmd.setOutput(&my);

    //----------------------------------------------------------------------------
    // add program-specific command line arguments
    //----------------------------------------------------------------------------

    // initialize arg
    ArgStruct arg;

    // input file: BAM read alignments
    ArgStruct argBam;
    arg = argBam; 
    arg.shortId = ""; 
    arg.longId = "bam"; 
    arg.description = "Read alignment input file (indexed and sorted BAM format)";
    arg.required = false; 
    arg.defaultValueString = ""; 
    arg.type = "string"; 
    arg.multi = false; 
    my.ArgList.push_back(arg);
    ValueArg<string> cmd_bam(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

    // input file: MOSAIK binary reference sequence archive
    ArgStruct argMbr;
    arg = argMbr; 
    arg.shortId = ""; 
    arg.longId = "fasta"; 
    arg.description = "Reference sequence fasta file";
    arg.required = true; 
    arg.defaultValueString = ""; 
    arg.type = "string"; 
    arg.multi = false; 
    my.ArgList.push_back(arg);
    ValueArg<string> cmd_fasta(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

    // input file: target regions 
    ArgStruct argTargets;
    arg = argTargets; 
    arg.shortId = ""; 
    arg.longId = "targets"; 
    arg.description = "Target region input file (BED format)";
    arg.required = false; 
    arg.defaultValueString = ""; 
    arg.type = "string"; 
    arg.multi = false; 
    my.ArgList.push_back(arg);
    ValueArg<string> cmd_targets(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

    // input file: target regions 
    ArgStruct argRegion;
    arg = argRegion; 
    arg.shortId = ""; 
    arg.longId = "region"; 
    arg.description = "Target region in <sequence>:<start>..<end> format";
    arg.required = false; 
    arg.defaultValueString = ""; 
    arg.type = "string"; 
    arg.multi = false; 
    my.ArgList.push_back(arg);
    ValueArg<string> cmd_region(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

    // input file: list of samples to analyze
    ArgStruct argSamples;
    arg = argSamples; 
    arg.shortId = ""; 
    arg.longId = "samples"; 
    arg.description = "File containing list of samples to analyze";
    arg.required = false; 
    arg.defaultValueString = ""; 
    arg.type = "string"; 
    arg.multi = false; 
    my.ArgList.push_back(arg);
    ValueArg<string> cmd_samples(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

    // output file: VCF output file
    ArgStruct argOutputFile;
    arg = argOutputFile; 
    arg.shortId = ""; 
    arg.longId = "outputFile"; 
    arg.description = "output file (format selectable via --output option)";
    arg.required = false; 
    arg.defaultValueString = ""; 
    arg.type = "string"; 
    arg.multi = false; 
    my.ArgList.push_back(arg);
    ValueArg<string> cmd_outputFile(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

    // output option: output information about all alleles contributing to a given genotyping position
    ArgStruct argOutput;
    arg = argOutput; 
    arg.shortId = ""; 
    arg.longId = "output"; 
    arg.description = "Output option: stdout output format, json or vcf";
    arg.required = false; 
    arg.defaultValueString = "vcf"; 
    arg.type = "string"; 
    arg.multi = false;
    my.ArgList.push_back(arg);
    ValueArg<string> cmd_output(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

    // output option: output information about all alleles contributing to a given genotyping position
    ArgStruct argOutputAlleles;
    arg = argOutputAlleles; 
    arg.shortId = ""; 
    arg.longId = "outputAlleles"; 
    arg.description = "Output option: output information about all alleles contributing to a given genotyping position";
    arg.required = false; 
    arg.defaultValueString = "false"; 
    arg.type = "switch"; 
    arg.multi = false;
    my.ArgList.push_back(arg);
    SwitchArg cmd_outputAlleles(arg.shortId, arg.longId, arg.description, cmd, false);

    // algorithm, use gigabayes method for calculating data likelihoods
    ArgStruct argBamBayesDataLikelihoods;
    arg = argBamBayesDataLikelihoods; 
    arg.shortId = ""; 
    arg.longId = "bamBayesDataLikelihoods"; 
    arg.description = "algorithm, use GigaBayes method for calculating data likelihoods";
    arg.required = false; 
    arg.defaultValueString = "false"; 
    arg.type = "switch"; 
    arg.multi = false;
    my.ArgList.push_back(arg);
    SwitchArg cmd_bamBayesDataLikelihoods(arg.shortId, arg.longId, arg.description, cmd, false);

    // ignore duplicate reads
    ArgStruct argUseDuplicateReads;
    arg = argUseDuplicateReads; 
    arg.shortId = ""; 
    arg.longId = "useDuplicateReads"; 
    arg.description = "Use alignments which are marked as duplicates in the analysis, default behavior is to drop them.";
    arg.required = false; 
    arg.defaultValueString = "false"; 
    arg.type = "switch"; 
    arg.multi = false;
    my.ArgList.push_back(arg);
    SwitchArg cmd_useDuplicateReads(arg.shortId, arg.longId, arg.description, cmd, false);

    // output option: output information about all alleles contributing to a given genotyping position
    ArgStruct argSuppressOutput;
    arg = argSuppressOutput; 
    arg.shortId = ""; 
    arg.longId = "suppressOutput"; 
    arg.description = "Output option: suppress output stream (for debugging work)";
    arg.required = false; 
    arg.defaultValueString = "false"; 
    arg.type = "switch"; 
    arg.multi = false;
    my.ArgList.push_back(arg);
    SwitchArg cmd_suppressOutput(arg.shortId, arg.longId, arg.description, cmd, false);

    // log file
    ArgStruct argLog;
    arg = argLog; 
    arg.shortId = ""; 
    arg.longId = "log"; 
    arg.description = "Log file";
    arg.required = false; 
    arg.defaultValueString = ""; 
    arg.type = "string"; 
    arg.multi = false;
    my.ArgList.push_back(arg);
    ValueArg<string> cmd_log(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

    // log file
    ArgStruct argTrace;
    arg = argTrace; 
    arg.shortId = ""; 
    arg.longId = "trace"; 
    arg.description = "Trace file output algorithmic traces to this file";
    arg.required = false; 
    arg.defaultValueString = ""; 
    arg.type = "string"; 
    arg.multi = false;
    my.ArgList.push_back(arg);
    ValueArg<string> cmd_traceFile(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

    // useRefAllele
    ArgStruct argUseRefAllele;
    arg = argUseRefAllele;
    arg.shortId = "";
    arg.longId = "useRefAllele";
    arg.description = "Use reference sequence allele in polymorphism calling?";
    arg.required = false;
    arg.defaultValueString = "false";
    arg.type = "switch";
    arg.multi = false;
    my.ArgList.push_back(arg);
    SwitchArg cmd_useRefAllele(arg.shortId, arg.longId, arg.description, cmd, false);

    // forceRefAllele
    ArgStruct argForceRefAllele;
    arg = argForceRefAllele;
    arg.shortId = "";
    arg.longId = "forceRefAllele";
    arg.description = "Force reference sequence allele to be always considered?";
    arg.required = false;
    arg.defaultValueString = "false";
    arg.type = "switch";
    arg.multi = false;
    my.ArgList.push_back(arg);
    SwitchArg cmd_forceRefAllele(arg.shortId, arg.longId, arg.description, cmd, false);

    // best N alleles
    ArgStruct argUseBestNAlleles;
    arg = argUseBestNAlleles;
    arg.shortId = "";
    arg.longId = "useBestNAlleles";
    arg.description = "Consider this many best possible alleles genotypes in analysis, even those without evidence.";
    arg.required = false;
    arg.defaultValueString = "0";
    arg.type = "int";
    arg.multi = false;
    my.ArgList.push_back(arg);
    ValueArg<int> cmd_useBestNAlleles(arg.shortId, arg.longId, arg.description, arg.required, 0, arg.type, cmd);

    // MQR: reference sequence mapping quality value
    ArgStruct argMQR;
    arg = argMQR;
    arg.shortId = "";
    arg.longId = "MQR";
    arg.description = "Reference sequence mapping quality value";
    arg.required = false;
    arg.defaultValueString = "100";
    arg.type = "int";
    arg.multi = false;
    my.ArgList.push_back(arg);
    ValueArg<int> cmd_MQR(arg.shortId, arg.longId, arg.description, arg.required, 100, arg.type, cmd);

    // BQR: reference base quality value
    ArgStruct argBQR;
    arg = argBQR;
    arg.shortId = "";
    arg.longId = "BQR";
    arg.description = "Reference sequence base quality value";
    arg.required = false;
    arg.defaultValueString = "60";
    arg.type = "int";
    arg.multi = false;
    my.ArgList.push_back(arg);
    ValueArg<int> cmd_BQR(arg.shortId, arg.longId, arg.description, arg.required, 60, arg.type, cmd);

    // ploidy: sample ploidy
    ArgStruct argPloidy;
    arg = argPloidy;
    arg.shortId = "";
    arg.longId = "ploidy";
    arg.description = "Sample ploidy";
    arg.required = false;
    arg.defaultValueString = "2";
    arg.type = "int";
    arg.multi = false;
    my.ArgList.push_back(arg);
    ValueArg<int> cmd_ploidy(arg.shortId, arg.longId, arg.description, arg.required, 2, arg.type, cmd);

    // sample: naming scheme for matching reads to samples
    ArgStruct argSampleNaming;
    arg = argSampleNaming;
    arg.shortId = "";
    arg.longId = "sampleNaming";
    arg.description = "Naming scheme for matching reads to samples";
    arg.required = false;
    arg.defaultValueString = "multiple";
    arg.type = "string";
    arg.multi = false;
    vector<string> allowedSampleNaming;
    allowedSampleNaming.push_back("single");
    allowedSampleNaming.push_back("multiple");
    allowedSampleNaming.push_back("trio");
    allowedSampleNaming.push_back("unknown");
    arg.constraint = allowedSampleNaming;
    ValuesConstraint<string> allowedSampleNamingVals(allowedSampleNaming); 
    my.ArgList.push_back(arg);
    ValueArg<string> cmd_sampleNaming(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, &allowedSampleNamingVals, cmd);

    // sampleDel
    ArgStruct argSampleDel;
    arg = argSampleDel;
    arg.shortId = "";
    arg.longId = "sampleDel";
    arg.description = "Delimiter string separating sample name and read name";
    arg.required = false;
    arg.defaultValueString = "-";
    arg.type = "string";
    arg.multi = false;
    my.ArgList.push_back(arg);
    ValueArg<string> cmd_sampleDel(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

    // MQL0: minimum mapping quality value to consider read
    ArgStruct argMQL0;
    arg = argMQL0;
    arg.shortId = "";
    arg.longId = "MQL0";
    arg.description = "Minimum mapping quality value to consider read";
    arg.required = false;
    arg.defaultValueString = "40";
    arg.type = "int";
    arg.multi = false;
    my.ArgList.push_back(arg);
    ValueArg<int> cmd_MQL0(arg.shortId, arg.longId, arg.description, arg.required, 40, arg.type, cmd);

    // BQL0: minimum read base quality value
    ArgStruct argBQL0;
    arg = argBQL0;
    arg.shortId = "";
    arg.longId = "BQL0";
    arg.description = "Minimum read base quality value";
    arg.required = false;
    arg.defaultValueString = "10";
    arg.type = "int";
    arg.multi = false;
    my.ArgList.push_back(arg);
    ValueArg<int> cmd_BQL0(arg.shortId, arg.longId, arg.description, arg.required, 10, arg.type, cmd);

    // MQL1: minimum mapping quality value required for at least one read for each allele
    ArgStruct argMQL1;
    arg = argMQL1;
    arg.shortId = "";
    arg.longId = "MQL1";
    arg.description = "Minimum mapping quality value required for at least one read for each allele";
    arg.required = false;
    arg.defaultValueString = "40";
    arg.type = "int";
    arg.multi = false;
    my.ArgList.push_back(arg);
    ValueArg<int> cmd_MQL1(arg.shortId, arg.longId, arg.description, arg.required, 40, arg.type, cmd);

    // BQL1: minimum base quality value required for at least one base for each allele
    ArgStruct argBQL1;
    arg = argBQL1;
    arg.shortId = "";
    arg.longId = "BQL1";
    arg.description = "Minimum read base quality value required for at least one base for each allele";
    arg.required = false;
    arg.defaultValueString = "30";
    arg.type = "int";
    arg.multi = false;
    my.ArgList.push_back(arg);
    ValueArg<int> cmd_BQL1(arg.shortId, arg.longId, arg.description, arg.required, 30, arg.type, cmd);

    // BQL2: minimum mismatch base quality value for read mismatch filtering
    ArgStruct argBQL2;
    arg = argBQL2;
    arg.shortId = "";
    arg.longId = "BQL2";
    arg.description = "Minimum mismatch base quality value for read mismatch filtering";
    arg.required = false;
    arg.defaultValueString = "10";
    arg.type = "int";
    arg.multi = false;
    my.ArgList.push_back(arg);
    ValueArg<int> cmd_BQL2(arg.shortId, arg.longId, arg.description, arg.required, 10, arg.type, cmd);

    // RMU: maximum number of mismatches between read and reference sequence
    ArgStruct argRMU;
    arg = argRMU;
    arg.shortId = "";
    arg.longId = "RMU";
    arg.description = "Maximum number of mismatches between read and refrence sequence";
    arg.required = false;
    arg.defaultValueString = "10000000";
    arg.type = "int";
    arg.multi = false;
    my.ArgList.push_back(arg);
    ValueArg<int> cmd_RMU(arg.shortId, arg.longId, arg.description, arg.required, 10000000, arg.type, cmd);

    // IDW: Window of base exclusion around INDELs in reads
    ArgStruct argIDW;
    arg = argIDW;
    arg.shortId = "";
    arg.longId = "IDW";
    arg.description = "Window of base exclusion around INDELs in reads";
    arg.required = false;
    arg.defaultValueString = "-1";  // this means 'no exclusion', whilst 0 means 'exclude indels'
    arg.type = "int";
    arg.multi = false;
    my.ArgList.push_back(arg);
    ValueArg<int> cmd_IDW(arg.shortId, arg.longId, arg.description, arg.required, -1, arg.type, cmd);

    // RDF: read dependence factor
    ArgStruct argRDF;
    arg = argRDF;
    arg.shortId = "";
    arg.longId = "RDF";
    arg.description = "Read dependence factor";
    arg.required = false;
    arg.defaultValueString = "0.9";
    arg.type = "double";
    arg.multi = false;
    my.ArgList.push_back(arg);
    ValueArg<double> cmd_RDF(arg.shortId, arg.longId, arg.description, arg.required, 0.9, arg.type, cmd);

    // minAltFraction: filter against positions, don't process
    // if no sample has less than this fraction of alternate alleles
    ArgStruct argMinAltFraction;
    arg = argMinAltFraction;
    arg.shortId = "";
    arg.longId = "minAltFraction";
    arg.description = "skip positions if no sample has >= than this fraction of alternate alleles";
    arg.required = false;
    arg.defaultValueString = "0.0";
    arg.type = "double";
    arg.multi = false;
    my.ArgList.push_back(arg);
    ValueArg<double> cmd_minAltFraction(arg.shortId, arg.longId, arg.description, arg.required, 0.0, arg.type, cmd);

    // minAltCount: filter against positions, don't process if no sample has
    // less than this count of alternate alleles
    ArgStruct argMinAltCount;
    arg = argMinAltCount;
    arg.shortId = "";
    arg.longId = "minAltCount";
    arg.description = "skip positions if no sample has >= this count of alternate alleles";
    arg.required = false;
    arg.defaultValueString = "1";
    arg.type = "int";
    arg.multi = false;
    my.ArgList.push_back(arg);
    ValueArg<int> cmd_minAltCount(arg.shortId, arg.longId, arg.description, arg.required, 1, arg.type, cmd);

    // TH: pairwise nucleotide diversity (theta)
    ArgStruct argTH;
    arg = argTH;
    arg.shortId = "";
    arg.longId = "TH";
    arg.description = "Pairwise nucleotide diversity (theta)";
    arg.required = false;
    arg.defaultValueString = "10E-3";
    arg.type = "double";
    arg.multi = false;
    my.ArgList.push_back(arg);
    ValueArg<double> cmd_TH(arg.shortId, arg.longId, arg.description, arg.required, 1E-3, arg.type, cmd);

    // PVL: minimum P(VAR) value for site to be reported in output
    ArgStruct argPVL;
    arg = argPVL;
    arg.shortId = "";
    arg.longId = "PVL";
    arg.description = "minimum P(VAR) value for site to be reported in output";
    arg.required = false;
    arg.defaultValueString = "0.0";
    arg.type = "double";
    arg.multi = false;
    my.ArgList.push_back(arg);
    ValueArg<double> cmd_PVL(arg.shortId, arg.longId, arg.description, arg.required, 0.0, arg.type, cmd);

    // algorithm
    ArgStruct argAlgorithm;
    arg = argAlgorithm;
    arg.shortId = "";
    arg.longId = "algorithm";
    arg.description = "P(SNP) calculation algorithm";
    arg.required = false;
    arg.defaultValueString = "banded";
    arg.type = "string";
    arg.multi = false;
    vector<string> allowedAlgorithm;
    allowedAlgorithm.push_back("banded");
    allowedAlgorithm.push_back("recursive");
    arg.constraint = allowedAlgorithm;
    ValuesConstraint<string> allowedAlgorithmVals(allowedAlgorithm); 
    my.ArgList.push_back(arg);
    ValueArg<string> cmd_algorithm(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, &allowedAlgorithmVals, cmd);

    // WB: bandwidth for banded P(SNP) calculation algorithm (genotype combination level)
    ArgStruct argWB;
    arg = argWB;
    arg.shortId = "";
    arg.longId = "WB";
    arg.description = "Bandwidth for banded P(SNP) calculation algorithm (genotype combination level)";
    arg.required = false;
    arg.defaultValueString = "2";
    arg.type = "int";
    arg.multi = false;
    my.ArgList.push_back(arg);
    ValueArg<int> cmd_WB(arg.shortId, arg.longId, arg.description, arg.required, 2, arg.type, cmd);

    // TB: number of terms to keep in banded P(SNP) calculation algorithm
    ArgStruct argTB;
    arg = argTB;
    arg.shortId = "";
    arg.longId = "TB";
    arg.description = "Number of terms to keep per cycle in banded P(SNP) calculation algorithm (0: keep all)";
    arg.required = false;
    arg.defaultValueString = "10";
    arg.type = "int";
    arg.multi = false;
    my.ArgList.push_back(arg);
    ValueArg<int> cmd_TB(arg.shortId, arg.longId, arg.description, arg.required, 10, arg.type, cmd);

    // includeMonoB
    ArgStruct argIncludeMonoB;
    arg = argIncludeMonoB;
    arg.shortId = "";
    arg.longId = "IncludeMonoB";
    arg.description = "Include monomorphic genotype combinations in initial list in banded SNP prob calculation?";
    arg.required = false;
    arg.defaultValueString = "false";
    arg.type = "switch";
    arg.multi = false;
    my.ArgList.push_back(arg);
    SwitchArg cmd_includeMonoB(arg.shortId, arg.longId, arg.description, cmd, false);

    // TR: number of terms to keep in recursive P(SNP) calculation algorithm
    ArgStruct argTR;
    arg = argTR;
    arg.shortId = "";
    arg.longId = "TR";
    arg.description = "Number of terms to keep per cycle of recursive P(SNP) calculation algorithm (0: keep all)";
    arg.required = false;
    arg.defaultValueString = "10";
    arg.type = "int";
    arg.multi = false;
    my.ArgList.push_back(arg);
    ValueArg<int> cmd_TR(arg.shortId, arg.longId, arg.description, arg.required, 10, arg.type, cmd);

    // I: report interval for debugging printouts
    ArgStruct argI;
    arg = argI;
    arg.shortId = "";
    arg.longId = "I";
    arg.description = "Report interval for debugging printouts";
    arg.required = false;
    arg.defaultValueString = "10000";
    arg.type = "int";
    arg.multi = false;
    my.ArgList.push_back(arg);
    ValueArg<int> cmd_I(arg.shortId, arg.longId, arg.description, arg.required, 10000, arg.type, cmd);

    // debug
    ArgStruct argRecord;
    arg = argRecord;
    arg.shortId = "";
    arg.longId = "record";
    arg.description = "Record messages to logfile?";
    arg.required = false;
    arg.defaultValueString = "false";
    arg.type = "switch";
    arg.multi = false;
    my.ArgList.push_back(arg);
    SwitchArg cmd_record(arg.shortId, arg.longId, arg.description, cmd, false);

    // debug
    ArgStruct argDebug;
    arg = argDebug;
    arg.shortId = "";
    arg.longId = "debug";
    arg.description = "Print debugging messages?";
    arg.required = false;
    arg.defaultValueString = "false";
    arg.type = "switch";
    arg.multi = false;
    my.ArgList.push_back(arg);
    SwitchArg cmd_debug(arg.shortId, arg.longId, arg.description, cmd, false);

    // debug2
    ArgStruct argDebug2;
    arg = argDebug2;
    arg.shortId = "";
    arg.longId = "debug2";
    arg.description = "Print very detailed debugging messages?";
    arg.required = false;
    arg.defaultValueString = "false";
    arg.type = "switch";
    arg.multi = false;
    my.ArgList.push_back(arg);
    SwitchArg cmd_debug2(arg.shortId, arg.longId, arg.description, cmd, false);

    //----------------------------------------------------------------------------
    // register (but not add to cmd) special arguments that are automatically
    // added to cmd
    //----------------------------------------------------------------------------

    // help
    ArgStruct argHelp;
    arg = argHelp;
    arg.shortId = "h";
    arg.longId = "help";
    arg.description = "Print usage statement?";
    arg.required = false;
    arg.defaultValueString = "false";
    arg.type = "switch";
    arg.multi = false;
    my.ArgList.push_back(arg);

    // version
    ArgStruct argVersion;
    arg = argVersion;
    arg.shortId = "";
    arg.longId = "version";
    arg.description = "Print program version?";
    arg.required = false;
    arg.defaultValueString = "false";
    arg.type = "switch";
    arg.multi = false;
    my.ArgList.push_back(arg);

    //----------------------------------------------------------------------------
    // parse command line and catch possible errors
    //----------------------------------------------------------------------------
    try {
        cmd.parse(argc,argv);
    } 
    catch ( ArgException& e ) { 
        cerr << "ERROR: " << e.error() << " " << e.argId() << endl; 
    }
  
    //----------------------------------------------------------------------------
    // assign command line parameters
    //----------------------------------------------------------------------------

    bam = cmd_bam.getValue();
    boost::split(bams, bam, boost::is_any_of(" \n"));

    fasta = cmd_fasta.getValue();
    targets = cmd_targets.getValue();
    region = cmd_region.getValue();
    samples = cmd_samples.getValue();
    outputFile = cmd_outputFile.getValue();
    log = cmd_log.getValue();
    output = cmd_output.getValue();
    outputAlleles = cmd_outputAlleles.getValue();
    bamBayesDataLikelihoods = cmd_bamBayesDataLikelihoods.getValue();
    traceFile = cmd_traceFile.getValue();
    useDuplicateReads = cmd_useDuplicateReads.getValue();
    suppressOutput = cmd_suppressOutput.getValue();

    useRefAllele = cmd_useRefAllele.getValue();
    forceRefAllele = cmd_forceRefAllele.getValue();
    useBestNAlleles = cmd_useBestNAlleles.getValue();
    MQR = cmd_MQR.getValue();
    BQR = cmd_BQR.getValue();
    ploidy = cmd_ploidy.getValue();
    sampleNaming = cmd_sampleNaming.getValue();
    sampleDel = cmd_sampleDel.getValue();
    MQL0 = cmd_MQL0.getValue();
    BQL0 = cmd_BQL0.getValue();
    MQL1 = cmd_MQL1.getValue();
    BQL1 = cmd_BQL1.getValue();
    BQL2 = cmd_BQL2.getValue();
    RMU = cmd_RMU.getValue();
    IDW = cmd_IDW.getValue();
    TH = (long double)cmd_TH.getValue();
    PVL = (long double)cmd_PVL.getValue();
    algorithm = cmd_algorithm.getValue();
    RDF = (long double)cmd_RDF.getValue();
    minAltFraction = (long double) cmd_minAltFraction.getValue();
    minAltCount = cmd_minAltCount.getValue();
    WB = cmd_WB.getValue();
    TB = cmd_TB.getValue();
    includeMonoB = cmd_includeMonoB.getValue();
    TR = cmd_TR.getValue();
    I = cmd_I.getValue();
    record = cmd_record.getValue();
    debug = cmd_debug.getValue();
    debug2 = cmd_debug2.getValue();

    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------
    // check and fix command line options
    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------

  
    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------
    // derived variables
    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------

    record = false;
    if (log != "")
        record = true;

    trace = false;
    if (traceFile != "")
        trace = true;


}
