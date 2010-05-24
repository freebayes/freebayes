//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// bamBayes
// PolyBayes for high-throughput next-generation sequences
// Copyright 2007, 2008, 2009 Gabor T. Marth, Boston College
// All rights reserved.
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//
//
//  Overview:  (man page)
//
//  --- this a snp caller
//  --- it operates on (an) alignment file, a reference, and a target file (bed
//  file) and a set of samples (optionally defined external file)
//
//  operation:
//
//  for every position in the target file that corresponds to regions of the
//  reference sequence,
//
//  take the set of observed bases from the alignments and estimate the
//  probability that the set of observations represents a snp,
//
//  specifically we estimate the probability that there are at least two or
//  more reference alleles present at a specific location (more general
//  solution would look at more than 2, but we use only two).
//
//  input reads are filtered from analysis on a variety of parameters:
//
//  - mapping quality value (if low, the read may not map to the region on the
//  ref genome)  --- output by mosaik
//
//  - number of differences between the read and the reference.
//
//  base (in read) level filters:
//
//  - if the quality value is too low, we might not want to use that base in
//  that read at all.
//
//
//  algorithmic outline:
//
//  question is, is the number of true alleles > 2 (psnp)
//
//  we evaluate this question using bayes theorem.  in the context of this
//  problem the theorem states that the probability of a snp event in our data
//  given the reference should approximate the 'data likelihood' (the
//  probability that our data reflects an event if one has occurred) multiplied
//  by the probability of the event occurring (our prior probability a
//  probability we estimate using data from population genetics).
//
//  
//
//  caveats:
//
//  the algorithmic complexity is 10^n for n individuals.
//  to resolve this issue we estimate the most likely data likelihood.
//  guess which are the highest terms (the dominant terms)
//  
//  keep the term which gives the highest data likelihood, but then vary the
//  genotype by some bandwidth (say 1 or 2).
//

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// includes
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// standard includes
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <time.h>

// "tclap" commandline parsing library
#include <tclap/CmdLine.h>

// "boost" regular expression library
#include <boost/regex.hpp>

// "boost" string manipulation
#include <boost/algorithm/string/join.hpp>

// "hash_map" true hashes
#include <ext/hash_map>

// private libraries
#include "Class-GigReader.h"
#include "Function-Sequence.h"
#include "Function-Generic.h"
#include "Function-Math.h"
#include "Class-BedReader.h"
#include "Class-FastaReader.h"
//#include "BamReader.h"
#include "BamMultiReader.h"
#include "ReferenceSequenceReader.h"
#include "Fasta.h"
#include "TryCatch.h"

// uses
using namespace std; 
using namespace TCLAP; 
using namespace __gnu_cxx;
using namespace BamTools;

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// typedefs
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

struct ArgStruct {
  string shortId;
  string longId;
  string description;
  bool required;
  string defaultValueString;
  int defaultValueInt;
  double defaultValueDouble;
  bool defaultValueBool;
  string type;
  bool multi;
  vector<string> constraint;
};

struct AlleleData {
  string base;
  short qual;
};

struct CigarElement {
  string type;
  unsigned int length;
};

typedef vector<CigarElement> Cigar;

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// template headers
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

template< typename keyType, typename valueType >
vector<keyType> sortKeysByValue(map<keyType, valueType, std::less<keyType> > hash, bool descend);

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// subroutine headers
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Prints a record into GFF format output file.
//------------------------------------------------------------------------------
bool printVar2Glt(
		  ofstream & gltFile,
		  string conName,
		  int posBegin,
		  int posEnd,
		  Variation & var
		  );

//------------------------------------------------------------------------------
// Prints a variation record into GFF format output file.
//------------------------------------------------------------------------------
bool printVar2Gff(
		  ofstream & gffFile,
		  int O,
		  string conName,
		  string ProgramName,
		  string allele1,
		  string allele2,
		  int posPadded,
		  int posBegin,
		  int posEnd,
		  Variation & var,
		  map<int, map<string, vector<Basecall> > > & individualBasecalls
		  );

bool printVar(
	      int O,
	      string conName,
	      string ProgramName,
	      string allele1,
	      string allele2,
	      int posPadded,
	      int posBegin,
	      int posEnd,
	      Variation & var,
	      map<int, map<string, vector<Basecall> > > & individualBasecalls
	      );

//------------------------------------------------------------------------------
// register an alignment into coverage data structures
//------------------------------------------------------------------------------
bool registerAlignment(
		       BamAlignment, 
		       map<int, string > &, 
		       map<int, string > &,
		       const int leftBound,
		       const int rightBound
		       );

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// static variables
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

static string ProgramName("bamBayes");
static string ProgramDescription("Bayesian SNP and short-INDEL polymorphism discovery program.");
static string ProgramVersion("0.7.2");
static string ProgramDate("2009-08-27");
static string ProgramDeveloper("Gabor T. Marth");
static string ProgramInstitution("Boston College");
static string ProgramCopyrightDates("2007, 2008, 2009.");

static vector<ArgStruct> ArgList;

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// classes
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

class MyOutput : public StdOutput {
public:
  
  virtual void failure(CmdLineInterface& c, ArgException& e)
  {
    cerr << "################################################################################" << endl;
    cerr << "### Program: " << ProgramName << " Version: " <<  ProgramVersion << " Release date: " << ProgramDate << endl;
    cerr << "### " << ProgramDescription << endl;
    cerr << "### " << "Copyright: " << ProgramCopyrightDates << " " << ProgramDeveloper << ", " << ProgramInstitution << endl;
    cerr << "### All rights reserved." << endl;
    cerr << "###" << endl;
    cerr << "### Command line error:" << e.what() << endl;
    cerr << "### For usage please type: " << c.getProgramName() << " --help" << endl;
    cerr << "################################################################################" << endl;
  }
  
  virtual void usage(CmdLineInterface& c)
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
      if (arg.longId != "") {idString += "--" + arg.longId;}
      if (arg.shortId != "") {idString += ", -" + arg.shortId;}

      string multiString = "single-valued";
      if (arg.multi) {multiString = "multi-valued";}

      if (arg.required) {
	cout << "### " << idString << " [" << arg.type << ", required, no default, " << multiString << "]" << endl;
      }
      else {
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
  
  virtual void version(CmdLineInterface& c)
  {
    cerr << "################################################################################" << endl;
    cerr << "### Program: " << ProgramName << " Version: " <<  ProgramVersion << " Release date: " << ProgramDate << endl;
    cout << "### " << ProgramDescription << endl;
    cout << "### " << "Copyright: " << ProgramCopyrightDates << " " << ProgramDeveloper << ", " << ProgramInstitution << endl;
    cout << "### All rights reserved." << endl;
    cout << "###" << endl;
    cerr << "################################################################################" << endl;
  }
};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// main
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

int main (int argc, char *argv[]) {

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
  CmdLine cmd("", ' ', ProgramVersion);
    
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
  arg.description = "Read alignment input file or files (indexed and sorted BAM format, multiple files must refer to the same reference)";
  arg.required = false; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false; 
  ArgList.push_back(arg);
  ValueArg<string> cmd_bam(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  //  TODO change this to fasta input
  //  ... try to find a .fai file corresponding to the fasta file
  //  unless one is supplied via another command line option
  //  if one doesn't exist, raise an error and say 'hey, convert it using samtools'

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
  ArgList.push_back(arg);
  ValueArg<string> cmd_fasta(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // TODO verify that Gabor is ok with this removal prior to release
  // input file: MOSAIK binary reference sequence archive
  /*
  ArgStruct argMbr;
  arg = argMbr; 
  arg.shortId = ""; 
  arg.longId = "mbr"; 
  arg.description = "Reference sequence archive file (MOSAIK binary format)";
  arg.required = true; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false; 
  ArgList.push_back(arg);
  ValueArg<string> cmd_mbr(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);
  */

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
  ArgList.push_back(arg);
  ValueArg<string> cmd_targets(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

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
  ArgList.push_back(arg);
  ValueArg<string> cmd_samples(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // output file: tab-delimited alignment and SNP report text file
  ArgStruct argRpt;
  arg = argRpt; 
  arg.shortId = ""; 
  arg.longId = "rpt"; 
  arg.description = "Output file: tab-delimited alignment and SNP report text file";
  arg.required = false; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<string> cmd_rpt(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // output file: VCF (variant call format) file
  ArgStruct argVcf;
  arg = argVcf; 
  arg.shortId = ""; 
  arg.longId = "vcf"; 
  arg.description = "Output file: VCF (variant call format) file";
  arg.required = true; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<string> cmd_vcf(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

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
  ArgList.push_back(arg);
  ValueArg<string> cmd_log(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

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
  ArgList.push_back(arg);
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
  ArgList.push_back(arg);
  SwitchArg cmd_forceRefAllele(arg.shortId, arg.longId, arg.description, cmd, false);

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
  ArgList.push_back(arg);
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
  ArgList.push_back(arg);
  ValueArg<int> cmd_BQR(arg.shortId, arg.longId, arg.description, arg.required, 60, arg.type, cmd);

  // ploidy: sample ploidy
  ArgStruct argPloidy;
  arg = argPloidy;
  arg.shortId = "";
  arg.longId = "ploidy";
  arg.description = "Sample ploidy";
  arg.required = false;
  arg.defaultValueString = "diploid";
  arg.type = "string";
  arg.multi = false;
  vector<string> allowedPloidy;
  allowedPloidy.push_back("haploid");
  allowedPloidy.push_back("diploid");
  arg.constraint = allowedPloidy;
  ValuesConstraint<string> allowedPloidyVals(allowedPloidy); 
  ArgList.push_back(arg);
  ValueArg<string> cmd_ploidy(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, &allowedPloidyVals, cmd);

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
  ArgList.push_back(arg);
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
  ArgList.push_back(arg);
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
  ArgList.push_back(arg);
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
  ArgList.push_back(arg);
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
  ArgList.push_back(arg);
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
  ArgList.push_back(arg);
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
  ArgList.push_back(arg);
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
  ArgList.push_back(arg);
  ValueArg<int> cmd_RMU(arg.shortId, arg.longId, arg.description, arg.required, 10000000, arg.type, cmd);

  // IDW: Window of base exclusion around INDELs in reads
  ArgStruct argIDW;
  arg = argIDW;
  arg.shortId = "";
  arg.longId = "IDW";
  arg.description = "Window of base exclusion around INDELs in reads";
  arg.required = false;
  arg.defaultValueString = "0";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_IDW(arg.shortId, arg.longId, arg.description, arg.required, 0, arg.type, cmd);

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
  ArgList.push_back(arg);
  ValueArg<double> cmd_RDF(arg.shortId, arg.longId, arg.description, arg.required, 0.9, arg.type, cmd);

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
  ArgList.push_back(arg);
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
  ArgList.push_back(arg);
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
  ArgList.push_back(arg);
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
  ArgList.push_back(arg);
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
  ArgList.push_back(arg);
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
  ArgList.push_back(arg);
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
  ArgList.push_back(arg);
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
  ArgList.push_back(arg);
  ValueArg<int> cmd_I(arg.shortId, arg.longId, arg.description, arg.required, 10000, arg.type, cmd);

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
  ArgList.push_back(arg);
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
  ArgList.push_back(arg);
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
  ArgList.push_back(arg);

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
  ArgList.push_back(arg);

  //----------------------------------------------------------------------------
  // add custom output handler
  //----------------------------------------------------------------------------
  MyOutput my;
  cmd.setOutput(&my);

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

  string bam = cmd_bam.getValue();
  //string mbr = cmd_mbr.getValue();
  string fasta = cmd_fasta.getValue();
  string targets = cmd_targets.getValue();
  string samples = cmd_samples.getValue();
  string rpt = cmd_rpt.getValue();
  string vcf = cmd_vcf.getValue();
  string log = cmd_log.getValue();

  bool useRefAllele = cmd_useRefAllele.getValue();
  bool forceRefAllele = cmd_forceRefAllele.getValue();
  int MQR = cmd_MQR.getValue();
  int BQR = cmd_BQR.getValue();
  string ploidy = cmd_ploidy.getValue();
  string sampleNaming = cmd_sampleNaming.getValue();
  string sampleDel = cmd_sampleDel.getValue();
  int MQL0 = cmd_MQL0.getValue();
  int BQL0 = cmd_BQL0.getValue();
  int MQL1 = cmd_MQL1.getValue();
  int BQL1 = cmd_BQL1.getValue();
  int BQL2 = cmd_BQL2.getValue();
  int RMU = cmd_RMU.getValue();
  int IDW = cmd_IDW.getValue();
  long double TH = (long double)cmd_TH.getValue();
  long double PVL = (long double)cmd_PVL.getValue();
  string algorithm = cmd_algorithm.getValue();
  long double RDF = (long double)cmd_RDF.getValue();
  int WB = cmd_WB.getValue();
  int TB = cmd_TB.getValue();
  bool includeMonoB = cmd_includeMonoB.getValue();
  int TR = cmd_TR.getValue();
  int I = cmd_I.getValue();
  bool debug = cmd_debug.getValue();
  bool debug2 = cmd_debug2.getValue();

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // check and fix command line options
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  
  // check that we have one output file specified
  if (rpt == "" && vcf == "") {
      cerr << "You must specify at least one output file (via --rpt or --vcf)" << endl;
      exit(1);
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // derived variables
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  bool diploid = false;
  if (ploidy == "diploid") {
    diploid = true;
  }

  bool banded = false;
  if (algorithm == "banded") {
    banded = true;
  }

  bool record = false;
  if (log != "") {
    record = true;
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // open log file if required
  //----------------------------------------------------------------------------
  //---------------------------------------------------------------------------- 
  ofstream logFile(log.c_str(), ios::out);
  if (record) {
    // report
    if (debug) {cerr << "Opening log file: " << log << " ...";}

    // open
    if (!logFile) {
      cerr << " unable to open file: " << log << endl;
      exit(1);
    }
    if (debug) {cerr << " done." << endl;}
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // report command line and parameters used
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  map<bool, string, less<bool> > bool2String;
  bool2String[false] = "false";
  bool2String[true] = "true";

  if (record) {
    logFile << "Command line:";
    for (int i=0; i<argc; i++) {
      logFile << " " << argv[i];
    }
    logFile << endl;
    logFile << endl;
    logFile << "Complete list of parameter values:" << endl;
    logFile << "  --bam = " << bam << endl;
    logFile << "  --fasta = " << fasta << endl;
    logFile << "  --targets = " << targets << endl;
    logFile << "  --samples = " << samples << endl;
    logFile << "  --rpt = " << rpt << endl;
    logFile << "  --log = " << log << endl;
    logFile << "  --useRefAllele = " <<  bool2String[useRefAllele] << endl;
    logFile << "  --forceRefAllele = " <<  bool2String[forceRefAllele] << endl;
    logFile << "  --MQR = " << MQR << endl;
    logFile << "  --BQR = " << BQR << endl;
    logFile << "  --ploidy = " << ploidy << endl;
    logFile << "  --sampleNaming = " << sampleNaming << endl;
    logFile << "  --sampleDel = " << sampleDel << endl;
    logFile << "  --MQL0 = " << MQL0 << endl;
    logFile << "  --BQL0 = " << BQL0 << endl;
    logFile << "  --MQL1 = " << MQL1 << endl;
    logFile << "  --BQL1 = " << BQL1 << endl;
    logFile << "  --BQL2 = " << BQL2 << endl;
    logFile << "  --RMU = " << RMU << endl;
    logFile << "  --IDW = " << IDW << endl;
    logFile << "  --TH = " << TH << endl;
    logFile << "  --PVL = " << PVL << endl;
    logFile << "  --algorithm = " << algorithm << endl;
    logFile << "  --RDF = " << RDF << endl;
    logFile << "  --WB = " << WB << endl;
    logFile << "  --TB = " << TB << endl;
    logFile << "  --includeMonoB = " <<  bool2String[includeMonoB] << endl;
    logFile << "  --TR = " << TR << endl;
    logFile << "  --I = " << I << endl;
    logFile << "  --debug = " <<  bool2String[debug] << endl;
    logFile << "  --debug2 = " <<  bool2String[debug2] << endl;
    logFile << endl;
  }
  if (debug) {
    cerr << "Command line:";
    for (int i=0; i<argc; i++) {
      cerr << " " << argv[i];
    }
    cerr << endl;
    cerr << endl;
    cerr << "Complete list of parameter values:" << endl;
    cerr << "  --bam = " << bam << endl;
    cerr << "  --fasta = " << fasta << endl;
    cerr << "  --targets = " << targets << endl;
    cerr << "  --samples = " << samples << endl;
    cerr << "  --rpt = " << rpt << endl;
    cerr << "  --log = " << log << endl;
    cerr << "  --useRefAllele = " <<  bool2String[useRefAllele] << endl;
    cerr << "  --forceRefAllele = " <<  bool2String[forceRefAllele] << endl;
    cerr << "  --MQR = " << MQR << endl;
    cerr << "  --BQR = " << BQR << endl;
    cerr << "  --ploidy = " << ploidy << endl;
    cerr << "  --sampleNaming = " << sampleNaming << endl;
    cerr << "  --sampleDel = " << sampleDel << endl;
    cerr << "  --BQL0 = " << BQL0 << endl;
    cerr << "  --MQL0 = " << MQL0 << endl;
    cerr << "  --BQL1 = " << BQL1 << endl;
    cerr << "  --MQL1 = " << MQL1 << endl;
    cerr << "  --BQL2 = " << BQL2 << endl;
    cerr << "  --RMU = " << RMU << endl;
    cerr << "  --IDW = " << IDW << endl;
    cerr << "  --TH = " << TH << endl;
    cerr << "  --PVL = " << PVL << endl;
    cerr << "  --algorithm = " << algorithm << endl;
    cerr << "  --RDF = " << RDF << endl;
    cerr << "  --WB = " << WB << endl;
    cerr << "  --TB = " << TB << endl;
    cerr << "  --includeMonoB = " <<  bool2String[includeMonoB] << endl;
    cerr << "  --TR = " << TR << endl;
    cerr << "  --I = " << I << endl;
    cerr << "  --debug = " <<  bool2String[debug] << endl;
    cerr << "  --debug2 = " <<  bool2String[debug2] << endl;
    cerr << endl;
  }


  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  // open BAM input file
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  
  // report
  
  // open
  //string  bamFileNameC = bam.c_str();
  BamMultiReader bReader;
  vector<string> bams;
  boost::split(bams, bam, boost::is_any_of(" \n"));
  if (bams.size() == 1) {
      if (record) {logFile << "Opening BAM fomat alignment input file: " << bam << " ...";}
      if (debug) {cerr << "Opening BAM format alignment input file: " << bam << " ...";}
  } else if (bams.size() > 1) {
      if (record) {logFile << "Opening BAM fomat alignment input files: " << endl;}
      if (debug) {cerr << "Opening BAM format alignment input files: " << endl;}
      for (vector<string>::const_iterator b = bams.begin(); b != bams.end(); ++b) {
          if (record) {logFile << *b << endl;}
          if (debug) {cerr << *b << endl;}
      }
  }
  bReader.Open(bams);
  // report
  if (record) {logFile << " done." << endl;}
  if (debug) {cerr << " done." << endl;}
  
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  // read sample list file or get sample names from bam file header
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  //
  // If a sample file is given, use it.  But otherwise process the bam file
  // header to get the sample names.
  //
  vector<string> sampleList;
  if (samples != "") {
    ifstream sampleFile(samples.c_str(), ios::in);
    if (! sampleFile) {
      cerr << "unable to open file: " << samples << endl;
      exit(1);
    }
    boost::regex patternSample("^(\\S+)\\s*(.*)$");
    boost::regex re("\\s+");
    boost::regex pr("^(\\S+):(\\S+)$");
    boost::smatch match;
    string line;
    while (getline(sampleFile, line)) {
      
      // if proper line
      if (boost::regex_search(line, match, patternSample)) {
	
	// assign content
	string s = match[1];
    if (debug) cerr << "found sample " << s << endl;
	sampleList.push_back(s);
      }
    }
  } 
      // retrieve header information

  vector<string> sampleListFromBam;

  string bamHeader = bReader.GetHeaderText();

  vector<string> headerLines;
  boost::split(headerLines, bamHeader, boost::is_any_of("\n"));

  for (vector<string>::const_iterator it = headerLines.begin(); it != headerLines.end(); ++it) {

      // get next line from header, skip if empty
      string headerLine = *it;
      if ( headerLine.empty() ) { continue; }

      // lines of the header look like:
      // "@RG     ID:-    SM:NA11832      CN:BCM  PL:454"
      //                     ^^^^^^^\ is our sample name
      if ( headerLine.find("@RG") == 0 ) {
          vector<string> readGroupParts;
          boost::split(readGroupParts, headerLine, boost::is_any_of("\t "));
          vector<string> nameParts;
          boost::split(nameParts, readGroupParts.at(2), boost::is_any_of(":"));
          string name = nameParts.back();
          //mergedHeader.append(1, '\n');
          if (debug) cerr << "found sample " << name << endl;
          sampleListFromBam.push_back(name);
      }
  }
   // no samples file given, read from BAM file header for sample names
  if (sampleList.size() == 0) {
      if (debug) cerr << "no sample list file given, reading sample names from bam file" << endl;
      for (vector<string>::const_iterator s = sampleListFromBam.begin(); s != sampleListFromBam.end(); ++s) {
          if (debug) cerr << "found sample " << *s << endl;
          sampleList.push_back(*s);
      }
  } else {
      // verify that the samples in the sample list are present in the bam,
      // and raise an error and exit if not
      for (vector<string>::const_iterator s = sampleList.begin(); s != sampleList.end(); ++s) {
          bool in = false;
          for (vector<string>::const_iterator b = sampleListFromBam.begin(); b != sampleListFromBam.end(); ++b) {
              if (*s == *b) { in = true; break; }
          }
          if (!in) {
              cerr << "ERROR: sample " << *s << " listed in sample file "
                  << samples.c_str() << " is not listed in the header of BAM file "
                  << bam << endl;
              exit(1);
          }
      }
  }
  
  //--------------------------------------------------------------------------
  // read reference sequences from input file
  //--------------------------------------------------------------------------
  
  // store the names of all the reference sequences in the BAM file
  RefVector refDatas = bReader.GetReferenceData();
  
  // data structs
  map<string, int > RefLength, RefId;
  for (RefVector::const_iterator refIter = refDatas.begin();
       refIter != refDatas.end(); refIter++) {
    RefData refData = *refIter;
    string refName = refData.RefName;
    RefLength[refName] = refData.RefLength;
    RefId[refName] = bReader.GetReferenceID(refName);
  }
  
  if (record) {
    logFile << "Number of ref seqs: " << bReader.GetReferenceCount() << endl;
  }
  if (debug) {
    cerr << "Number of ref seqs: " << bReader.GetReferenceCount() << endl;
  }

  //--------------------------------------------------------------------------
  // process input fasta file
  //--------------------------------------------------------------------------
  // This call loads the reference and reads any index file it can find.
  // If it can't find an index file for the reference, it will attempt to
  // generate one alongside it.
  FastaReference* fastaReference = new FastaReference(fasta);
  map<string, unsigned int > refseqLength;
  map<string, string > refseqDna;

  //--------------------------------------------------------------------------
  // load ref seq names into hash
  //--------------------------------------------------------------------------
  for(vector<FastaIndexEntry>::const_iterator it = fastaReference->index->begin(); 
          it != fastaReference->index->end(); ++it) {

    FastaIndexEntry entry = *it;

    if (debug) cerr << "reading sequence... " << entry << endl;

    // we split out the first word of the full sequence name for keying our sequences
    vector<string> sequenceNameParts;
    boost::split(sequenceNameParts, entry.name, boost::is_any_of(" "));
    string name = sequenceNameParts.front();

    refseqLength[name] = entry.length;
    // get the reference base sequence
    refseqDna[name] = fastaReference->getSequence(entry.name);

  }

  // close the reference sequence reader
  delete fastaReference;
  
  // report success
  if (record) {logFile << " done." << endl;}
  if (debug) {logFile << " done." << endl;}
  
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  // process target region input file if required
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  
  // reference-specific target lists
  // XXX this need not be a map between a string and a vector; it could go
  // vector to vector, as the reference sequences are referred to by number in
  // BamTools anyway
  map<string, vector<BedData> > RefTargetList;
  
  // if target file specified use targets from file
  int tc = 0;
  if (targets != "") {
    
    //------------------------------------------------------------------------
    // open input BED file if required
    //------------------------------------------------------------------------
    
    // report
    if (record) {logFile << "Making BedReader object for target file: " << targets << " ...";}
    if (debug) {cerr << "Making BedReader object for target file: " << targets << " ...";}
    
    // make
    BedReader * br = new BedReader(targets);
    
    if (! br->isOpen()) {
      cerr << "Unable to open target file: " << targets << "... terminating." << endl;
      exit(1);
    }
    
    //------------------------------------------------------------------------
    // iterate through entries
    //------------------------------------------------------------------------
    BedData bd;
    while (br->getNextEntry(bd)) {
      if (debug2) {
	cerr << bd.seq << "\t" << bd.left << "\t" << bd.right << "\t" << bd.desc << endl;
      }
      if (RefId.count(bd.seq) > 0) {
	if (bd.left < 1 || bd.right > RefLength[bd.seq]) {
	  cerr << "Target region coordinate outside of reference sequence bounds... terminating." << endl;
	  exit(1);
	}
	RefTargetList[bd.seq].push_back(bd);
	tc++;
      }
    }
    
    //------------------------------------------------------------------------
    // close
    //------------------------------------------------------------------------
    br->close();

    if (debug) {cerr << "done" << endl;}
    if (record) {logFile << "done" << endl;}
  }

  // otherwise analyze all reference sequences from BAM file
  else {
    RefVector::iterator refIter = refDatas.begin();
    RefVector::iterator refEnd  = refDatas.end();
    for( ; refIter != refEnd; ++refIter) {
      RefData refData = *refIter;
      string refName = refData.RefName;
      int refId = bReader.GetReferenceID(refName);
      int refLength = refData.RefLength;
      BedData bd;
      bd.seq = refName;
      bd.left = 1;
      bd.right = refLength;
      RefTargetList[bd.seq].push_back(bd);
      tc++;
    }
  }
  if (debug) {cerr << "Number of target regions: " << tc << endl;}
  if (record) {logFile << "Number of target regions: " << tc << endl;}
  

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // open output file(s)
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // open report output file
  //----------------------------------------------------------------------------

  // report
  if (record) {logFile << "opening report output file for writing: " << rpt << "...";}
  if (debug) {cerr << "opening report output file for writing: " << rpt << "...";}

  // open output streams
  ofstream rptFile, vcfFile;
  bool outputRPT, outputVCF; // for legibility

  if (rpt != "") {
      outputRPT = true;
      rptFile.open(rpt.c_str());
      if (!rptFile) {
        if (record) {logFile << " unable to open file: " << rpt << endl;}
        cerr << " unable to open file: " << rpt << endl;
        exit(1);
      }
  } else { outputRPT = false; }

  if (vcf != "") {
      outputVCF = true;
      vcfFile.open(vcf.c_str());
      if (!vcfFile) {
        if (record) {logFile << " unable to open file: " << vcf << endl;}
        cerr << " unable to open file: " << vcf << endl;
        exit(1);
      }
  } else { outputVCF = false; }
  if (record) {logFile << " done." << endl;}
  if (debug) {cerr << " done." << endl;}

  //----------------------------------------------------------------------------
  // write header information
  //----------------------------------------------------------------------------
  if (outputRPT) {
      rptFile << "# Command line that generated this output:";
      for (int i=0; i<argc; i++) {
        rptFile << " " << argv[i];
      }
      rptFile << endl;
      rptFile << "#" << endl;
      rptFile << "# Complete list of parameter values:" << endl;
      rptFile << "#   --bam = " << bam << endl;
      rptFile << "#   --fasta = " << fasta << endl;
      rptFile << "#   --targets = " << targets << endl;
      rptFile << "#   --samples = " << samples << endl;
      rptFile << "#   --rpt = " << rpt << endl;
      rptFile << "#   --log = " << log << endl;
      rptFile << "#   --useRefAllele = " <<  bool2String[useRefAllele] << endl;
      rptFile << "#   --forceRefAllele = " <<  bool2String[forceRefAllele] << endl;
      rptFile << "#   --MQR = " << MQR << endl;
      rptFile << "#   --BQR = " << BQR << endl;
      rptFile << "#   --ploidy = " << ploidy << endl;
      rptFile << "#   --sampleNaming = " << sampleNaming << endl;
      rptFile << "#   --sampleDel = " << sampleDel << endl;
      rptFile << "#   --BQL0 = " << BQL0 << endl;
      rptFile << "#   --MQL0 = " << MQL0 << endl;
      rptFile << "#   --BQL1 = " << BQL1 << endl;
      rptFile << "#   --MQL1 = " << MQL1 << endl;
      rptFile << "#   --BQL2 = " << BQL2 << endl;
      rptFile << "#   --RMU = " << RMU << endl;
      rptFile << "#   --IDW = " << IDW << endl;
      rptFile << "#   --TH = " << TH << endl;
      rptFile << "#   --PVL = " << PVL << endl;
      rptFile << "#   --algorithm = " << algorithm << endl;
      rptFile << "#   --RDF = " << RDF << endl;
      rptFile << "#   --WB = " << WB << endl;
      rptFile << "#   --TB = " << TB << endl;
      rptFile << "#   --includeMonoB = " <<  bool2String[includeMonoB] << endl;
      rptFile << "#   --TR = " << TR << endl;
      rptFile << "#   --I = " << I << endl;
      rptFile << "#   --debug = " <<  bool2String[debug] << endl;
      rptFile << "#   --debug2 = " <<  bool2String[debug2] << endl;
      rptFile << "#" << endl;
  }

  
  if (outputVCF) {
      time_t rawtime;
      struct tm * timeinfo;
      char datestr [80];

      time(&rawtime);
      timeinfo = localtime(&rawtime);

      strftime(datestr, 80, "%Y%m%d %X", timeinfo);

      vcfFile << "##format=VCFv3.3" << endl
              << "##fileDate=" << datestr << endl
              << "##source=gigabayes" << endl
              << "##reference=1000GenomesPilot-NCBI36" << endl
              << "##phasing=none" << endl
              << "##notes=\"All FORMAT fields matching *i* (e.g. NiBAll, NiA) refer to individuals.\"" << endl
             
              << "##INFO=NS,1,Integer,\"total number of samples\"" << endl
              << "##INFO=ND,1,Integer,\"total number of non-duplicate samples\"" << endl
              << "##INFO=DP,1,Integer,\"total read depth at this base\"" << endl
              << "##INFO=AC,1,Integer,\"total number of alternate alleles in called genotypes\"" << endl
              //<< "##INFO=AN,1,Integer,\"total number of alleles in called genotypes\"" << endl

              // these are req'd
              << "##FORMAT=GT,1,String,\"Genotype\"" << endl // g
              << "##FORMAT=GQ,1,Integer,\"Genotype Quality\"" << endl // phred prob of genotype
              << "##FORMAT=DP,1,Integer,\"Read Depth\"" << endl // NiBAll[ind]
              << "##FORMAT=HQ,2,Integer,\"Haplotype Quality\"" << endl
              << "##FORMAT=QiB,1,Integer,\"Total base quality\"" << endl
              << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" 
              << boost::algorithm::join(sampleList, "\t")
              << endl;

  }
  //-------------------------------------------------------------------------
  //-------------------------------------------------------------------------
  // process target regions
  //-------------------------------------------------------------------------
  //-------------------------------------------------------------------------

  //-------------------------------------------------------------------------
  // global coverage stat variables
  //-------------------------------------------------------------------------
  long int gT = 0;
  long int gL = 0;
  long int gAAll = 0;
  long int gANondup = 0;
  long int gA = 0;
  long int gDAll = 0;
  long int gDNondup = 0;
  long int gD = 0;

  // outer loop reference sequences in list
  for (map<string, vector<BedData> >::const_iterator refIter = RefTargetList.begin(); 
       refIter != RefTargetList.end(); refIter++) {
    string refName = refIter->first;
    vector<BedData> targetRegions = refIter->second;
    for(vector<BedData>::const_iterator targetIter = targetRegions.begin();
	targetIter != targetRegions.end(); targetIter++) {
      
      // retrieve target region
      BedData target = * targetIter;
      
      // find target ref seq in list of ref seqs from BAM file
      int refId = bReader.GetReferenceID(target.seq);
      RefData refData = refDatas[refId];
      int refLength = refData.RefLength;
      
      // report
      if (record) {logFile << "  Processing target: " << target.seq << ":" << target.left << "-" << target.right << endl;}
      if (debug) {cerr << "  Processing target: " << target.seq << ":" << target.left << "-" << target.right << endl;}

      //----------------------------------------------------------------------
      // coverage structures
      //----------------------------------------------------------------------
      long int tL = target.right - target.left + 1;
      long int tAAll = 0;
      long int tANondup = 0;
      long int tA = 0;
      long int tDAll = 0;
      long int tDNondup = 0;
      long int tD = 0;

      //----------------------------------------------------------------------
      // data structure for aligned position specific allele layout
      //----------------------------------------------------------------------
      
      // contig-position specific coverage quantities
      // XXX all these maps can be replaced by one data structure
      //     mapping a range of ints into a vector of structures which hold these
      // XXX these are unused??
      map<int, int > readCount, readCountPlus, readCountMinus;
      map<int, map_string_int > alleleCount, alleleCountPlus, alleleCountMinus;
      map<int, map_string_longDouble > alleleQual, alleleQualPlus, alleleQualMinus;
      
      // contig-position and individual specific coverage quantities
      // XXX this map as well
      map<int, map<string, vector<Basecall> > >
          individualBasecalls, 
          individualBasecallsAll,
          individualBasecallsNondup;
      // it should be like:
      // Contig contig(reference, start_position, end_position);
      // initialize it with our current position,
      //  and access would be like:
      // contig[position][individual]
      
      // collect targets for which there is no read data
      bool badTarget = false;
      int numReadsInTarget = 0;    
      BamAlignment ba;

      // skip target if invalid reference
      if (!(refId < bReader.GetReferenceCount())) {
        badTarget = true;
        if (record) {logFile << "      WARNING: Target ref ID invalid: " << refId << endl;}
        if (debug) {cerr << "      WARNING: Target ref ID invalid: " << refId << endl;}	
      }
      else {
        if (record) {logFile << "      Target ref ID OK: " << refId << endl;}
        if (debug) {cerr << "      Target ref ID OK: " << refId << endl;}	
      }
      
      // skip target if ref seq does not have any alignments
      if (! refData.RefHasAlignments) {
        badTarget = true;
        if (record) {logFile << "      WARNING: Target has no alignments..." << endl;}
        if (debug) {cerr << "      WARNING: Target has no alignments..." << endl;}	
      }
      else {
        if (record) {logFile << "      Target has alignments..." << endl;}
        if (debug) {cerr << "      Target has alignments..." << endl;}	
      }
      
      if (record) {logFile << "      Jumping to target start in BAM file. RefId: " << refId << " pos: " << target.left << endl;}
      if (debug) {cerr << "      Jumping to target start in BAM file. RefId: " << refId << " pos: " << target.left << endl;}	
      if (! bReader.Jump(refId, target.left)) {
        badTarget = true;
        if (record) {logFile << "      WARNING: Cannot jump to target start in BAM file. RefId: " << refId << " pos: " << target.left << endl;}
        if (debug) {cerr << "      WARNING: Cannot jump to target start in BAM file. REfId: " << refId << " pos: " << target.left << endl;}	
      }
      else {
        if (record) {logFile << "      Jumped to target start in BAM file. RefId: " << refId << " pos: " << target.left << endl;}
        if (debug) {cerr << "      Jumped to target start in BAM file. RefId: " << refId << " pos: " << target.left << endl;}	

        /* 
         * Get the number of alignments in this target space, and jump back to
         * the start.  We do this so that we can properly step through the
         * positions at the end of the region that do not have reads following
         * them without having to loop through them one by one.
         */
        while (bReader.GetNextAlignment(ba) && (ba.Position + 1 <= target.right)) {
            ++numReadsInTarget;
        }
        if (debug) cerr << "numReadsInTarget " << numReadsInTarget << endl;
        bReader.Jump(refId, target.left); // jump back
      }

      //------------------------------------------------------------------
      // if target is bad, warn
      //------------------------------------------------------------------
      if (badTarget) {
        if (record) {logFile << "    Target not OK... adding default stats" << endl;}
        if (debug) {cerr << "    Target not OK... adding default stats" << endl;}

        cerr << "    Warning: Bad target starting at position " << target.left << endl;

      }
 
      // otherwise analyze target
      else {
	if (record) {logFile << "    Target OK... processing..." << endl;}
	if (debug) {cerr << "    Target OK... processing..." << endl;}
      
    int processedReadsInTarget = 0;
	int refPosLast = target.left - 1;
	while (bReader.GetNextAlignment(ba) && (ba.Position + 1 <= target.right)) {
        ++processedReadsInTarget;
	  
	  // only process if mapped
	  if (! ba.IsMapped()) {
	    continue;
	  }
	  
	  // mark read if duplicate
	  bool dupRead = false;
	  if (ba.IsDuplicate()) {
	    dupRead = true;
	  }
	  
	  // mark read if mapping quality is low
	  bool lowMapQualRead = false;
	  if (ba.MapQuality < MQL0) {
	    lowMapQualRead = true;
	  }
	  

      // extract sample name and information
      string readName = ba.Name;
      string sampleName;
      if (! ba.GetReadGroup(sampleName)) {
          cerr << "ERROR: Couldn't find read group id (@RG tag) for BAM Alignment " << ba.Name 
               << " EXITING!" << endl;
          exit(1);
               //<< " ... attempting to read from read name" << endl;
          //SampleInfo sampleInfo = extractSampleInfo(readName, sampleNaming, sampleDel);
          //string sampleName = sampleInfo.sampleId;
      }
	  
	  //----------------------------------------------------------------
	  // get mapping quality
	  //----------------------------------------------------------------
	  short readMq = ba.MapQuality;

	  //----------------------------------------------------------------
	  // parse cigar
	  //----------------------------------------------------------------
	  string rDna = ba.QueryBases;
	  string rQual = ba.Qualities;
	  
	  //----------------------------------------------------------------
	  // 1. filter on high quality mismatches between read and reference
	  // 2. register neighborhood of high quality INDELs in read
	  // 3. make stats on duplicate reads
	  //----------------------------------------------------------------

	  // number of mismatches between read and ref seq
	  int mru = 0;
	  
	  // indicator function of ref seq positions where this read
	  //   has an INDEL
      //   XXX this one can also get mapped into the above data structure
	  map<int, bool > readINDEL;

	  // initialize reference sequence position and read position
	  int sp = ba.Position + 1;
	  int rp = 1;
	  vector<CigarOp>::const_iterator cigarIter = ba.CigarData.begin();
	  vector<CigarOp>::const_iterator cigarEnd  = ba.CigarData.end();
	  for ( ; cigarIter != cigarEnd; ++cigarIter ) {
	    unsigned int l = (*cigarIter).Length;
	    char t = (*cigarIter).Type;
	    
	    if (t == 'S') { // soft clip
	      rp += l;
	    }
	    else if (t == 'M') { // match or mismatch
	      for (int i=1; i<=l; i++) {
		
		// extract aligned base
        string b;
		TRY { b = rDna.substr(rp-1, 1); } CATCH;
		char Q = rQual[rp-1];
		
		// convert base quality value into short int
		short q = static_cast<short>(Q) - 33;
		
		// get reference allele
        string sb;
		TRY { sb = refseqDna[refName].substr(sp-1, 1); } CATCH;

		// register mismatch
		if (b != sb && q >= BQL2) {mru++;}
		
		// update positions
		sp++;
		rp++;
	      }
	    }
	    else if (t == 'D' || t == 'N') { // deletion or skipped region
	      if (t == 'D') {
		// extract base quality of left and right flanking non-deleted bases
		char QL = rQual[rp-1];
		char QR = rQual[rp];
		
		// convert base quality values into short int
		short qL = static_cast<short>(QL) - 33;
		short qR = static_cast<short>(QR) - 33;
		
		// calculate maximum of two values
		short q = max(qL, qR);
		
		// iterate through each deleted base
		for (int i=1; i<=l; i++) {

		  // register deletion + base quality with reference sequence
		  if (q >= BQL2) {
		    readINDEL[sp] = true;
		  }

		  // update refseq position
		  sp++;
		}
	      }
	    }
	    else if (t == 'I') { // insertion
	      for (int i=1; i<=l; i++) {

		// extract base quality of inserted base
		char Q = rQual[rp-1];
		
		// convert base quality value into short int
		short q = static_cast<short>(Q) - 33;
		
		// register insertion + base quality with reference sequence
		if (q >= BQL2) {
		  readINDEL[sp] = true;
		  readINDEL[sp+1] = true;
		}

		// update read position
		rp += 1;
	      }	    
	    }
	  } // end cigar iter loop
	  
	  //----------------------------------------------------------------
	  // mark read if too many mismatches
	  //----------------------------------------------------------------
	  bool filtRead = false;
	  if (mru > RMU) {
	    filtRead = true;
	  }
	  
	  //----------------------------------------------------------------
	  // mask out refseq bases in a neighborhood of read INDELs for this
	  //   read
	  //----------------------------------------------------------------
	  
	  // mask data structure
      // XXX this can also be encapsulated
	  map<int, bool > readMask;

	  // !!! too simplistic! Replace with more sophisticated method!
	  for (map<int, bool >::const_iterator posIter = readINDEL.begin();
	       posIter != readINDEL.end(); ++posIter) {
	    int pos = posIter->first;
	    for(int i=pos-IDW; i<=pos+IDW; i++) {
	      readMask[i] = true;
	    }
	  }
	  
	  //----------------------------------------------------------------
	  // update target alignment counts
	  //----------------------------------------------------------------
	  tAAll++;
	  if (! dupRead) {
	    tANondup++;
	    if (!filtRead && !lowMapQualRead) {
	      tA++;
	    }
	  }
	  
	  //----------------------------------------------------------------
	  // register alignment in coverage structures
	  //----------------------------------------------------------------
	  
	  // initialize reference sequence position and read position
	  sp = ba.Position + 1;
	  rp = 1;
	  
	  // iterate through Cigar operations and register read alignment
	  cigarIter = ba.CigarData.begin();
	  cigarEnd  = ba.CigarData.end();
	  for ( ; cigarIter != cigarEnd; ++cigarIter ) {
	    unsigned int l = (*cigarIter).Length;
	    char t = (*cigarIter).Type;
	    
	    if (t == 'S') { // soft clip
	      rp += l;
	    }
	    else if (t == 'M') { // match or mismatch
	      for (int i=1; i<=l; i++) {
		
		// extract aligned base
        string b;
		TRY { b = rDna.substr(rp-1, 1); } CATCH;
		char Q = rQual[rp-1];
		
		// convert base quality value into short int
		short q = static_cast<short>(Q) - 33;
		
		// make base call data struct
		Basecall bc;
		bc.seqName = readName;
		bc.map = readMq;
		bc.base = b;
		bc.qual = q;
		bc.strand = "+";
		if (ba.IsReverseStrand()) {
		  bc.strand = "-";
		}
		
		// register base call
		individualBasecallsAll[sp][sampleName].push_back(bc);
		if (!dupRead) {
		  individualBasecallsNondup[sp][sampleName].push_back(bc);
		  if (!filtRead && !lowMapQualRead && q >= BQL0 && !readMask[sp]) {
		    individualBasecalls[sp][sampleName].push_back(bc);
		  }
		}
		
		// update positions
		sp++;
		rp++;
	      }
	    }
	    else if (t == 'D' || t == 'N') { // deletion or skipped region
	      // update refseq position
	      sp += l;
	    }
	    else if (t == 'I') { // insertion
	      // update read position
	      rp += l;	    
	    }
	  } // end read alignment registration loop
	  
	  //----------------------------------------------------------------
	  // analyze positions up to the position before the beginning of this alignment
      // including positions after the alignment if it is the last one in the target region
	  //----------------------------------------------------------------
	  int refPosLeft = refPosLast + 1;
	  if (refPosLeft < target.left) {refPosLeft = target.left;}

      if (debug) {
          if (processedReadsInTarget >= numReadsInTarget) { cerr << "at last read in target region" << endl; }
      }
      // check if we are at the last read
      // if so set the refPosRight to the target right
	  int refPosRight = (processedReadsInTarget < numReadsInTarget) ? min((int)ba.Position, target.right) : target.right;
	  for (int p = refPosLeft; p <= refPosRight; p++) {
	    
	    // report progress if necessary
	    if (record && (p % I == 0)) {
	      logFile << "    Processing refseq:" << refName << " position: " << p << endl;
	    }
	    if (debug && (p % I == 0)) {
	      cerr << "    Processing ref seq:" << refName << " position: " << p << endl;
	    }
	    
	    //---------------------------------------------------------------
	    // calculate coverage and determine if allele coverage condition
	    //   is met for each allele
	    //---------------------------------------------------------------

	    // allele coverage condition
	    bool okA = false;
	    bool okC = false;
	    bool okG = false;
	    bool okT = false;

	    // multi-individual coverage data
	    int ntBAll = 0, ntBNondup = 0;
	    int ntB = 0, qtB = 0, ntBp = 0, ntBm = 0;
	    int ntA = 0, qtA = 0, ntAp = 0, ntAm = 0; 
	    int ntC = 0, qtC = 0, ntCp = 0, ntCm = 0; 
	    int ntG = 0, qtG = 0, ntGp = 0, ntGm = 0; 
	    int ntT = 0, qtT = 0, ntTp = 0, ntTm = 0;
	    
	    // single-individual coverage maps
        // XXX this should be a unified data structure
        //   ... not 22 maps
	    map<string, int > NiBAll, NiBNondup;
	    map<string, int > NiB, QiB, NiBp, NiBm;
	    map<string, int > NiA, QiA, NiAp, NiAm;
	    map<string, int > NiC, QiC, NiCp, NiCm;
	    map<string, int > NiG, QiG, NiGp, NiGm;
	    map<string, int > NiT, QiT, NiTp, NiTm;
	    
	    for (vector<string>::const_iterator sampleIter = sampleList.begin();
		 sampleIter != sampleList.end(); sampleIter++) {
	      string ind = *sampleIter;
	      
	      // read coverage data
	      int niBAll = 0, niBNondup = 0;
	      int niB = 0, qiB = 0, niBp = 0, niBm = 0;
	      int niA = 0, qiA = 0, niAp = 0, niAm = 0; 
	      int niC = 0, qiC = 0, niCp = 0, niCm = 0; 
	      int niG = 0, qiG = 0, niGp = 0, niGm = 0; 
	      int niT = 0, qiT = 0, niTp = 0, niTm = 0;

	      // read coverage by ALL reads
	      vector<Basecall> basecallsAll = individualBasecallsAll[p][ind];
	      niBAll += basecallsAll.size();
	      
	      // read coverage by non-dup reads
	      vector<Basecall> basecallsNondup = individualBasecallsNondup[p][ind];
	      niBNondup += basecallsNondup.size();
	      
	      // read coverage by non-dup, non-filtered reads
	      vector<Basecall> basecalls = individualBasecalls[p][ind];
	      for (vector<Basecall>::const_iterator bcIter = basecalls.begin();
		   bcIter != basecalls.end(); bcIter++) {
		Basecall bc = *bcIter;
		short m = bc.map;
		string b = bc.base;
		short q = bc.qual;
		string s = bc.strand;
		
		// allele coverage condition check
		if (m >= MQL1 && q >= BQL1) {
		  if (b == "A") {
		    okA = true;
		  }
		  else if (b == "C") {
		    okC = true;
		  }
		  else if (b == "G") {
		    okG = true;
		  }
		  else if (b == "T") {
		    okT = true;
		  }
		}

		niB++;
		qiB += q;
		if (s == "+") {
		  niBp++;
		  if (b == "A") {
		    niA++; qiA += q; niAp++;
		  }
		  else if (b == "C") {
		    niC++; qiC += q; niCp++;
		  }
		  else if (b == "G") {
		    niG++; qiG += q; niGp++;
		  }
		  else if (b == "T") {
		    niT++; qiT += q; niTp++;
		  }
		}
		else if (s == "-") {
		  niBm++;
		  if (b == "A") {
		    niA++; qiA += q; niAm++;
		  }
		  else if (b == "C") {
		    niC++; qiC += q; niCm++;
		  }
		  else if (b == "G") {
		    niG++; qiG += q; niGm++;
		  }
		  else if (b == "T") {
		    niT++; qiT += q; niTm++;
		  }
		}
	      }

	      //-------------------------------------------------------------
	      // register individual coverage at this position
	      //-------------------------------------------------------------
	      NiBAll[ind] = niBAll; NiBNondup[ind] = niBNondup;
	      NiB[ind] = niB; QiB[ind] = qiB; NiBp[ind] = niBp; NiBm[ind] = niBm;
	      NiA[ind] = niA; QiA[ind] = qiA; NiAp[ind] = niAp; NiAm[ind] = niAm; 
	      NiC[ind] = niC; QiC[ind] = qiC; NiCp[ind] = niCp; NiCm[ind] = niCm; 
	      NiG[ind] = niG; QiG[ind] = qiG; NiGp[ind] = niGp; NiGm[ind] = niGm; 
	      NiT[ind] = niT; QiT[ind] = qiT; NiTp[ind] = niTp; NiTm[ind] = niTm;


	      //-------------------------------------------------------------
	      // update total sample coverage at this position
	      //-------------------------------------------------------------
	      ntBAll += niBAll; ntBNondup += niBNondup;
	      ntB += niB; qtB += qiB; ntBp += niBp; ntBm += niBm;
	      ntA += niA; qtA += qiA; ntAp += niAp; ntAm += niAm; 
	      ntC += niC; qtC += qiC; ntCp += niCp; ntCm += niCm; 
	      ntG += niG; qtG += qiG; ntGp += niGp; ntGm += niGm; 
	      ntT += niT; qtT += qiT; ntTp += niTp; ntTm += niTm;

	    }
	    
	    //---------------------------------------------------------------
	    // update target base coverage stats
	    //---------------------------------------------------------------
	    tDAll += ntBAll;
	    tDNondup += ntBNondup;
	    tD += ntB;
	    
	    //---------------------------------------------------------------
	    // add reference allele if needed
	    //---------------------------------------------------------------
	    
	    // set defaults
	    string sb = "?";
	    short sm = MQR;
	    short sq = BQR;
	    string sbPrev = "?";
	    string sbNext = "?";
	    
        // XXX very bad, why are we appending the ref seq to the end of this vector?
	    vector<string> sampleListPlusRef = sampleList;
	    
	    if (useRefAllele) {
	      
	      // fill real values
	      TRY { sb = refseqDna[refName].substr(p-1, 1); } CATCH;

	      //	      sb = RefFastaData[refName].sequence.substr(p-1, 1);	      
	      if (p > 1 ) {
            TRY { sbPrev = refseqDna[refName].substr(p-2, 1); } CATCH;
		//		sbPrev = RefFastaData[refName].sequence.substr(p-2, 1);
	      }
	      if (p < refLength) {
		//		sbNext = RefFastaData[refName].sequence.substr(p, 1);
            TRY { sbNext = refseqDna[refName].substr(p, 1); } CATCH;
	      }
	      
	      // only use ref allele if real base
	      if (sb == "A" || sb == "C" || sb == "G" || sb == "T") {

		// make reference base call
		Basecall bc;
		bc.seqName = target.seq;
		bc.map = sm;
		bc.base = sb;
		bc.qual = sq;
		bc.strand = "?";
	      
        // XXX
		sampleListPlusRef.push_back(target.seq);
		vector<Basecall> basecalls;
		basecalls.push_back(bc);
		individualBasecalls[p][target.seq] = basecalls;
	      }
	    }
	    
	    //---------------------------------------------------------------
	    // calculate genotype data likelihoods and data sufficiency
	    //---------------------------------------------------------------

	    // define per-sample genotype data likelihood
	    map<string, map<string, long double > > 
	      LnProbBiGivenGi;			

	    // define per-sample data sufficiency indicator
	    map<string, bool > DataSufficientGi;
	    
	    // calculate genotype data likelihoods and data sufficiency
	    for (vector<string>::const_iterator sampleIter = sampleListPlusRef.begin();
		 sampleIter != sampleListPlusRef.end(); sampleIter++) {
	      string ind = *sampleIter;
	      vector<Basecall> basecalls = individualBasecalls[p][ind];
	      
          // DATA LIKELIHOOD function
	      // calculate genotype log likelihoods
          //
          // maps genotype to probability for this individual
	      map<string, long double > logGl 
		= logGenotypeLikelihoods(  // this gets loaded into another map, indexed by individual
					 basecalls,
					 diploid,
					 RDF,
					 debug2
					 );
	      
	      // assign genotype log likelihoods in multi-individual data structure
	      LnProbBiGivenGi[ind] = logGl;

	      // determine data sufficiency
	      DataSufficientGi[ind] = false;
	      if (basecalls.size() > 0) {
		DataSufficientGi[ind] = true;
	      }
	    }	      
	    
	    //---------------------------------------------------------------
	    // find best two alleles
	    //---------------------------------------------------------------
	    
	    // load sorting hash
        map<string, int> AQ;
	    AQ["A"] = qtA; AQ["C"] = qtC; AQ["G"] = qtG; AQ["T"] = qtT;
	    
	    // sort eligible alleles according to descending allele quality
	    vector<string> sortedAlleles = sortKeysByValue(AQ, true);
	    
	    // find best alleles
	    string allele1 = sortedAlleles[0];
	    string allele2 = sortedAlleles[1];
	    
	    // if ref allele is to be used make sure that is one of the alleles
	    if (useRefAllele && forceRefAllele) {
	      if (sb != allele1 && sb != allele2) {
            allele2 = allele1;
            allele1 = sb;
	      }
	    }
	    

	    //---------------------------------------------------------------
	    // run posterior calculation for two best alleles
	    //---------------------------------------------------------------
	    Variation var = posteriorProb2(
					  LnProbBiGivenGi,
					  DataSufficientGi,
					  sampleListPlusRef,
					  diploid,
					  TH,
					  banded,
					  WB,
					  TB,
					  includeMonoB,
					  TR,
					  allele1,
					  allele2,
					  debug2
					  );

	    
	    //---------------------------------------------------------------
	    // report if required
	    //---------------------------------------------------------------
	    if (var.pSnp >= PVL) {

            if (outputRPT) {
              rptFile << "POSITION" << "\t" << refName << "\t" <<  p  
                  << "\t" << ntBAll << "\t" << ntBNondup
                  << "\t" << ntB << "\t" << qtB << "\t" << ntBp << "\t" << ntBm
                  << "\t" << ntA << "\t" << qtA << "\t" << ntAp << "\t" << ntAm 
                  << "\t" << ntC << "\t" << qtC << "\t" << ntCp << "\t" << ntCm 
                  << "\t" << ntG << "\t" << qtG << "\t" << ntGp << "\t" << ntGm 
                  << "\t" << ntT << "\t" << qtT << "\t" << ntTp << "\t" << ntTm;	    
              rptFile << "\t" << sb << "\t" << sq << "\t" << sbPrev << "\t" << sbNext;
            
              //-------------------------------------------------------------
              // print two best alleles and P(SNP) value
              //-------------------------------------------------------------
              rptFile << "\t" << allele1 << "\t" << allele2 << "\t" << var.pSnp;
              
              //-------------------------------------------------------------
              // print whether allele coverage conditions are met for each allele
              //-------------------------------------------------------------
              rptFile << "\t" << okA << "\t" << okC << "\t" << okG << "\t" << okT;

              //-------------------------------------------------------------
              // print number of individuals
              //-------------------------------------------------------------
              rptFile << "\t" << sampleList.size();
            
              //-------------------------------------------------------------
              // print individual coverage & posterior genotype probabilities
              //-------------------------------------------------------------
              for (vector<string>::const_iterator sampleIter = sampleList.begin();
               sampleIter != sampleList.end(); sampleIter++) {
                string ind = *sampleIter;
            
                // print individual id
                rptFile << "\t" << ind;
                
                // print coverage values
                rptFile  << "\t" << NiBAll[ind] << "\t" << NiBNondup[ind]
                     << "\t" << NiB[ind] << "\t" << QiB[ind] << "\t" << NiBp[ind] << "\t" << NiBm[ind]
                     << "\t" << NiA[ind] << "\t" << QiA[ind] << "\t" << NiAp[ind] << "\t" << NiAm[ind] 
                     << "\t" << NiC[ind] << "\t" << QiC[ind] << "\t" << NiCp[ind] << "\t" << NiCm[ind] 
                     << "\t" << NiG[ind] << "\t" << QiG[ind] << "\t" << NiGp[ind] << "\t" << NiGm[ind]
                     << "\t" << NiT[ind] << "\t" << QiT[ind] << "\t" << NiTp[ind] << "\t" << NiTm[ind];
            
                for (map<string, long double >::const_iterator 
                       gIter = var.individualGenotypeProbability[ind].begin();
                     gIter != var.individualGenotypeProbability[ind].end();
                     gIter++) {
                  string g = gIter->first;
                  long double p = gIter->second;
                  rptFile << "\t" << g << "\t" << p;
                }
              }
              rptFile << endl;
            }

            if (outputVCF) {

                vcfFile << refName << "\t" 
                    <<  p + 1 << "\t"
                    << "." << "\t"  // . is the default signifier for "no id"
                    << sb << "\t"; // reference base at this position
                //vcfFile << allele1 << "," << allele2 << "\t";
                // the VCF spec suggests that these shouldn't be the reference; but when we use forceRefAllele they can be
                if (allele1 == sb) {
                    vcfFile << allele2 << "\t";
                } else if (allele2 == sb) {
                    vcfFile << allele1 << "\t";
                } else {
                    vcfFile << allele1 << "," << allele2 << "\t";
                }

                vcfFile << phred(var.pSnp) << "\t" // quality of the snp call
                    << "." << "\t" // filters, . means "no data"
                    << "NS=" << sampleList.size() << ";" // number of samples
                    << "ND=" << ntBNondup << ";" // 
                    << "DP=" << tDAll
                    << "\t"
                    // format string
                    //<< "NiBAll:NiBNondup:NiB:QiB:NiBp:NiBm:NiA:QiA:NiAp:NiAm:NiC:QiC:NiCp:NiCm:NiG:QiG:NiGp:NiGm:NiT:QiT:NiTp:NiTm:GiP"
                    << "GT:GQ:DP" // FIXME this is incomplete
                    << "\t";
                // samples
                bool firstSampleEntry = true; // flag to help with tab insertion between fields
                for (vector<string>::const_iterator sampleIter = sampleList.begin();
                 sampleIter != sampleList.end(); sampleIter++) {
                  string ind = *sampleIter;
                
                  if (firstSampleEntry)
                      firstSampleEntry = false;
                  else
                      vcfFile << "\t";  // join fields on tab

                  // find the most likely genotype for this individual
                  long double max = 0;
                  string bestGt = "";
                  for (map<string, long double >::const_iterator 
                         gIter = var.individualGenotypeProbability[ind].begin();
                       gIter != var.individualGenotypeProbability[ind].end();
                       gIter++) {
                    string g = gIter->first;
                    long double p = gIter->second;
                    if (p > max) {
                        max = p;
                        bestGt = gIter->first;
                    }
                  }
                  // then print the best genotype score and coverage info
                  long double p = var.individualGenotypeProbability[ind][bestGt];
                  // format is diploid genotype=phred, or N/N=phred, e.g. 0/1=56
                  vcfFile << ((bestGt.at(0) == sb.at(0)) ? "0" : "1")
                      << "/" << ((bestGt.at(1) == sb.at(0)) ? "0" : "1")
                      << ":" << phred(p)
                      << ":" << NiBAll[ind];
                }
                vcfFile << endl;
            }
	    }

	    //--------------------------------------------------------------
	    // erase allele data for this position
	    //--------------------------------------------------------------
	    individualBasecallsAll.erase(p);
	    individualBasecallsNondup.erase(p);
	    individualBasecalls.erase(p);
	  }
	  // update last alignment position analyzed
	  // remember: ba.Position is 0 base
	  refPosLast = refPosRight;
	}
	
      }

      //----------------------------------------------------------------
      // update global alignment counts
      //----------------------------------------------------------------
      gAAll += tAAll;
      gANondup += tANondup;
      gA += tA;
      gDAll += tDAll;
      gDNondup += tDNondup;
      gD += tD;
      gT++;
      gL += tL;
    } // end bed targets loop
  } // end reference targets loop

  //--------------------------------------------------------------------
  // report global stats
  //--------------------------------------------------------------------
  if (debug) {
    cerr << "Total numnber of targets: " << gT << endl;
    cerr << "Total numnber of reads overlapping targets: " << gA << endl;
    cerr << "Total length of targets: " << gL << endl;
    cerr << "Total depth in targets: " << gD << endl;
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // finish up
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // report
  if (record) {logFile << "Closing files...";}
  if (debug) {cerr << "Closing files...";}
  
  
  // close BamReader
  bReader.Close();
	
  // close rpt file
  if (outputRPT) rptFile.close();
  if (outputVCF) vcfFile.close();
  
  // report
  if (record) {
    logFile << " done." << endl;
    logFile << "RPT output file: " << rpt << endl;
    logFile << "Log file: " << log << endl;
    logFile << "Program completed." << endl;
  }
  if (debug) {
    cerr << " done." << endl;
    cerr << "RPT output file: " << rpt << endl;
    if (record) {
      cerr << "Log file: " << log << endl;
    }
    cerr << "Program completed." << endl;
  }

  // close log file last
  if (record) {logFile.close();}
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// templates
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// sorts keys of a hash in order of associated value
//------------------------------------------------------------------------------
template< typename keyType, typename valueType >
vector<keyType> sortKeysByValue(map<keyType, valueType, std::less<keyType> > hash, bool descend) {

  // instantiate inverse hash as a multimap
  multimap<valueType, keyType, std::less<valueType> > inverseHash;

  // load elements of hash into inverseHash
  for (typename map<keyType, valueType, std::less<keyType> >::const_iterator iter = hash.begin(); iter != hash.end(); iter++) {
    keyType key = iter->first;
    valueType value = iter->second;
    inverseHash.insert(typename multimap<valueType, keyType, std::less<valueType> >::value_type(value, key));
  }

  // compose vector of original keys sorted by original values
  vector<keyType> sortedKeys;
  for(typename multimap<valueType, keyType, std::less<valueType> >::const_iterator iter = inverseHash.begin(); 
      iter != inverseHash.end(); iter++) {
    keyType key = iter->second;
    sortedKeys.push_back(key);
  }

  // reverse if descending order was required
  if (descend) {
    reverse(sortedKeys.begin(), sortedKeys.end());
  }

  // return
  return sortedKeys;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// subroutines
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Prints a record into GFF format output file.
//------------------------------------------------------------------------------
bool printVar2Glt(
		  ofstream & gltFile,
		  string conName,
		  int posBegin,
		  int posEnd,
		  Variation & var
		  ) {

  gltFile << conName 
	  << "\t" << posBegin 
	  << "\t" << posEnd;


  for (map<string, map<string, long double > >::const_iterator 
	 indIter = var.individualGenotypeLogDataLikelihood.begin(); 
       indIter != var.individualGenotypeLogDataLikelihood.end();
       indIter++) {
      
      string ind = indIter->first;
      map<string, long double > gl = indIter->second;

      gltFile << "\t" << ind << ":";
      
      bool first = true;
      for (map<string, long double >::const_iterator 
	     glIter = gl.begin();
	   glIter != gl.end();
	   glIter++) {
	string g = glIter->first;
	long double l = glIter->second;

	if (! first) {
	  gltFile << "&";
	}
	first = false;

	gltFile << g << "|" << l;
      }
    }
  gltFile << endl;
}

//------------------------------------------------------------------------------
// Prints a variation record into GFF format output file.
//------------------------------------------------------------------------------
bool printVar2Gff(
		  ofstream & gffFile,
		  int O,
		  string conName,
		  string ProgramName,
		  string allele1,
		  string allele2,
		  int posPadded,
		  int posBegin,
		  int posEnd,
		  Variation & var,
		  map<int, map<string, vector<Basecall> > > & individualBasecalls
		  ) {

  // write variation entry in output file
  gffFile << conName << "\t"
	  << ProgramName << "\t";
  if (allele1 == "-" || allele1 == "*" || allele2 == "-" || allele2 == "*") {
    gffFile << "INDEL" << "\t";
  }
  else {
    gffFile << "SNP" << "\t";
  }
  gffFile << posBegin << "\t"
	  << posEnd << "\t"
	  << var.pSnp << "\t"
	  << "." << "\t"
	  << "." << "\t"
	  << "alleles=" << allele1 << "," << allele2;
  if (O >= 1) {
    gffFile << ";individualGenotypes=";
    bool firstInd = true;
    for (map<string, map<string, long double > >::const_iterator 
	   indIter = var.individualGenotypeProbability.begin(); 
	 indIter != var.individualGenotypeProbability.end();
	 indIter++) {
      
      string ind = indIter->first;
      if (! firstInd) {
	gffFile << ",";
      }
      firstInd = false;
      vector<string> genotypes = sortKeysByValue(var.individualGenotypeProbability[ind], true);
      string g = genotypes[0];
      long double p = var.individualGenotypeProbability[ind][g];
      gffFile << ind << ":" << g << "|" << p;
    }
  }
  if (O >= 2) {
    gffFile << ";individualGenotypeProbabilities=";
    bool firstInd = true;
    for (map<string, map<string, long double > >::const_iterator 
	   indIter = var.individualGenotypeProbability.begin(); 
	 indIter != var.individualGenotypeProbability.end();
	 indIter++) {
      
      string ind = indIter->first;
      if (! firstInd) {
	gffFile << ",";
      }
      firstInd = false;
      gffFile << ind << ":";
      bool firstG = true;
      for (map<string, long double >::const_iterator 
	     gIter = var.individualGenotypeProbability[ind].begin();
	   gIter != var.individualGenotypeProbability[ind].end();
	   gIter++) {
	string g = gIter->first;
	long double p = gIter->second;
	if (! firstG) {
	  gffFile << "&";
	}
	firstG = false;
	gffFile << g << "|" << p;
      }
    }
  }
  if (O >= 3) {
    gffFile << ";individualAlleleCounts=";
    bool firstInd = true;
    for (map<string, vector<Basecall> >::const_iterator 
	   indIter = individualBasecalls[posPadded].begin();
	 indIter != individualBasecalls[posPadded].end(); indIter++) {
      
      // retreive individual
      string ind = indIter->first;
      vector<Basecall> baseCalls = indIter->second;
      if (! firstInd) {
	gffFile << ",";
      }
      firstInd = false;
      gffFile << ind << ":";
      
      int count1 = 0; int count2 = 0; int count1Pos = 0; int count1Min = 0; int count2Pos = 0; int count2Min = 0; 
      int qual1 = 0; int qual2 = 0; int qual1Pos = 0; int qual1Min = 0; int qual2Pos = 0; int qual2Min = 0;
      for(vector<Basecall>::const_iterator iter = baseCalls.begin();
	  iter != baseCalls.end(); iter++) {
	Basecall bc = *iter;
	
	// filter irrelevant base calls
	if ((bc.base != allele1) && (bc.base != allele2)) {
	  continue;
	}
	
	// update allele counts
	if (bc.base == allele1) {
	  count1++;
	  qual1 += bc.qual;
	  if (bc.strand == "+") {
	    count1Pos++;
	    qual1Pos += bc.qual;
	  }
	  else if (bc.strand == "-") {
	    count1Min++;
	    qual1Min += bc.qual;
	  }
	}
	else if (bc.base == allele2) {
	  count2++;
	  qual2 += bc.qual;
	  if (bc.strand == "+") {
	    count2Pos++;
	    qual2Pos += bc.qual;
	  }
	  else if (bc.strand == "-") {
	    count2Min++;
	    qual2Min += bc.qual;
	  }
	}
      }
      gffFile << allele1 << "|" << count1 << "|" << qual1 << "|+|" << count1Pos << "|" << qual1Pos << "|-|" << count1Min << "|" << qual1Min << "&" << allele2 << "|" << count2 << "|" << qual2 << "|+|" << count2Pos << "|" << qual2Pos << "|-|" << count2Min << "|" << qual2Min;
    }
  }
  if (O >= 4) {
    gffFile << ";individualReadAlleles=";
    bool firstInd = true;
    for (map<string, vector<Basecall> >::const_iterator 
	   indIter = individualBasecalls[posPadded].begin();
	 indIter != individualBasecalls[posPadded].end(); indIter++) {
      
      // retreive individual
      string ind = indIter->first;
      vector<Basecall> baseCalls = indIter->second;
      if (! firstInd) {
	gffFile << ",";
      }
      firstInd = false;
      gffFile << ind << ":";
      bool firstB = true;
      for(vector<Basecall>::const_iterator iter = baseCalls.begin();
	  iter != baseCalls.end(); iter++) {
	Basecall bc = *iter;
	
	string seqName = bc.seqName;
	string b = bc.base;
	short q = bc.qual;
	string s = bc.strand;
	if (! firstB) {
	  gffFile << "&";
	}
	firstB = false;
	gffFile << seqName << "|" << s << "|" << b << "|" << q;
      }
    }
  }
  gffFile << endl;

  //----------------------------------------------------------------------------
  // return
  //----------------------------------------------------------------------------
  return true;
}

//------------------------------------------------------------------------------
// Prints a variation record into GFF format output file.
//------------------------------------------------------------------------------
bool printVar(
		  int O,
		  string conName,
		  string ProgramName,
		  string allele1,
		  string allele2,
		  int posPadded,
		  int posBegin,
		  int posEnd,
		  Variation & var,
		  map<int, map<string, vector<Basecall> > > & individualBasecalls
		  ) {

  // write variation entry in output file
  cout << conName << "\t"
	  << ProgramName << "\t";
  if (allele1 == "-" || allele1 == "*" || allele2 == "-" || allele2 == "*") {
    cout << "INDEL" << "\t";
  }
  else {
    cout << "SNP" << "\t";
  }
  cout << posBegin << "\t"
	  << posEnd << "\t"
	  << var.pSnp << "\t"
	  << "." << "\t"
	  << "." << "\t"
	  << "alleles=" << allele1 << "," << allele2;
  if (O >= 1) {
    cout << ";individualGenotypes=";
    bool firstInd = true;
    for (map<string, map<string, long double > >::const_iterator 
	   indIter = var.individualGenotypeProbability.begin(); 
	 indIter != var.individualGenotypeProbability.end();
	 indIter++) {
      
      string ind = indIter->first;
      if (! firstInd) {
	cout << ",";
      }
      firstInd = false;
      vector<string> genotypes = sortKeysByValue(var.individualGenotypeProbability[ind], true);
      string g = genotypes[0];
      long double p = var.individualGenotypeProbability[ind][g];
      cout << ind << ":" << g << "|" << p;
    }
  }
  if (O >= 2) {
    cout << ";individualGenotypeProbabilities=";
    bool firstInd = true;
    for (map<string, map<string, long double > >::const_iterator 
	   indIter = var.individualGenotypeProbability.begin(); 
	 indIter != var.individualGenotypeProbability.end();
	 indIter++) {
      
      string ind = indIter->first;
      if (! firstInd) {
	cout << ",";
      }
      firstInd = false;
      cout << ind << ":";
      bool firstG = true;
      for (map<string, long double >::const_iterator 
	     gIter = var.individualGenotypeProbability[ind].begin();
	   gIter != var.individualGenotypeProbability[ind].end();
	   gIter++) {
	string g = gIter->first;
	long double p = gIter->second;
	if (! firstG) {
	  cout << "&";
	}
	firstG = false;
	cout << g << "|" << p;
      }
    }
  }
  if (O >= 3) {
    cout << ";individualAlleleCounts=";
    bool firstInd = true;
    for (map<string, vector<Basecall> >::const_iterator 
	   indIter = individualBasecalls[posPadded].begin();
	 indIter != individualBasecalls[posPadded].end(); indIter++) {
      
      // retreive individual
      string ind = indIter->first;
      vector<Basecall> baseCalls = indIter->second;
      if (! firstInd) {
	cout << ",";
      }
      firstInd = false;
      cout << ind << ":";
      
      int count1 = 0; int count2 = 0; int count1Pos = 0; int count1Min = 0; int count2Pos = 0; int count2Min = 0; 
      int qual1 = 0; int qual2 = 0; int qual1Pos = 0; int qual1Min = 0; int qual2Pos = 0; int qual2Min = 0;
      for(vector<Basecall>::const_iterator iter = baseCalls.begin();
	  iter != baseCalls.end(); iter++) {
	Basecall bc = *iter;
	
	// filter irrelevant base calls
	if ((bc.base != allele1) && (bc.base != allele2)) {
	  continue;
	}
	
	// update allele counts
	if (bc.base == allele1) {
	  count1++;
	  qual1 += bc.qual;
	  if (bc.strand == "+") {
	    count1Pos++;
	    qual1Pos += bc.qual;
	  }
	  else if (bc.strand == "-") {
	    count1Min++;
	    qual1Min += bc.qual;
	  }
	}
	else if (bc.base == allele2) {
	  count2++;
	  qual2 += bc.qual;
	  if (bc.strand == "+") {
	    count2Pos++;
	    qual2Pos += bc.qual;
	  }
	  else if (bc.strand == "-") {
	    count2Min++;
	    qual2Min += bc.qual;
	  }
	}
      }
      cout << allele1 << "|" << count1 << "|" << qual1 << "|+|" << count1Pos << "|" << qual1Pos << "|-|" << count1Min << "|" << qual1Min << "&" << allele2 << "|" << count2 << "|" << qual2 << "|+|" << count2Pos << "|" << qual2Pos << "|-|" << count2Min << "|" << qual2Min;
    }
  }
  if (O >= 4) {
    cout << ";individualReadAlleles=";
    bool firstInd = true;
    for (map<string, vector<Basecall> >::const_iterator 
	   indIter = individualBasecalls[posPadded].begin();
	 indIter != individualBasecalls[posPadded].end(); indIter++) {
      
      // retreive individual
      string ind = indIter->first;
      vector<Basecall> baseCalls = indIter->second;
      if (! firstInd) {
	cout << ",";
      }
      firstInd = false;
      cout << ind << ":";
      bool firstB = true;
      for(vector<Basecall>::const_iterator iter = baseCalls.begin();
	  iter != baseCalls.end(); iter++) {
	Basecall bc = *iter;
	
	string seqName = bc.seqName;
	string b = bc.base;
	short q = bc.qual;
	string s = bc.strand;
	if (! firstB) {
	  cout << "&";
	}
	firstB = false;
	cout << seqName << "|" << s << "|" << b << "|" << q;
      }
    }
  }
  cout << endl;

  //----------------------------------------------------------------------------
  // return
  //----------------------------------------------------------------------------
  return true;
}

//------------------------------------------------------------------------------
// register an alignment in coverage structures
//------------------------------------------------------------------------------
bool registerAlignment(
		       BamAlignment a, 
		       map<int, string > & pileDna, 
		       map<int, string > & pileQual,
		       const int leftBound,
		       const int rightBound
		       ) {

  //----------------------------------------------------------------------------
  // parse cigar and register alignment in coverage data structures
  //----------------------------------------------------------------------------
  string rDna = a.QueryBases;
  string rQual = a.Qualities;
  
  // initialize reference sequence position and read position
  int sp = a.Position;
  int rp = 1;
  
  // iterate through Cigar operations and register read alignment
  vector<CigarOp>::const_iterator cigarIter = a.CigarData.begin();
  vector<CigarOp>::const_iterator cigarEnd  = a.CigarData.end();
  for ( ; cigarIter != cigarEnd; ++cigarIter ) {
    unsigned int l = (*cigarIter).Length;
    char t = (*cigarIter).Type;
    
    if (t == 'S') { // soft clip
      rp += l;
    }
    else if (t == 'M') { // match or mismatch
      for (int i=1; i<=l; i++) {
	
	if (sp >= leftBound && sp <= rightBound) {
	  // extract aligned base
      string b, q;
	  TRY { b = rDna.substr(rp-1, 1); } CATCH;
	  TRY { q = rQual.substr(rp-1, 1); } CATCH;
	
	  // add aligned base to pile
	  pileDna[sp] += b;
	  pileQual[sp] += q;
	}
	
	// update positions
	sp++;
	rp++;
      }	    
    }
    else if (t == 'D' || t == 'N') { // deletion or skipped region
      for (int i=1; i<=l; i++) {
	
	if (sp >= leftBound && sp <= rightBound) {
	  // assign aligned base as gap
	  string b = "-";
	  string q = "!";
	  
	  // add gap base to pile
	  pileDna[sp] += b;
	  pileQual[sp] += q;
	}

	// update refseq position
	sp++;
      }
    }
    else if (t == 'I') { // insertion
      rp += l;	    
    }
  }
  //----------------------------------------------------------------------------
  // return status
  //----------------------------------------------------------------------------
  return true;
}
