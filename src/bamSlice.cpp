//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// bamSlice
// Writes out aligned base info for a given ref pos in a BAM file
// Copyright 2009 Gabor T. Marth, Boston College
// All rights reserved.
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

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

// "tclap" commandline parsing library
#include <tclap/CmdLine.h>

// "boost" regular expression library
#include <boost/regex.hpp>

// "hash_map" true hashes
#include <ext/hash_map>

// private libraries
#include "Function-Generic.h"
#include "Function-Math.h"
#include "BamReader.h"
#include "ReferenceSequenceReader.h"

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
//------------------------------------------------------------------------------
// static variables
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

static string ProgramName("bamSlice");
static string ProgramDescription("Reports base coverage data at a ref pos in a BAM file.");
static string ProgramVersion("0.0.1");
static string ProgramDate("2009-06-04");
static string ProgramDeveloper("Gabor T. Marth");
static string ProgramInstitution("Boston College");
static string ProgramCopyrightDates("2009.");

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
  arg.description = "Read alignment input file (indexed binary BAM format)";
  arg.required = true; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false; 
  ArgList.push_back(arg);
  ValueArg<string> cmd_bam(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // input file: MOSAIK binary reference sequence archive
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

  // ref
  ArgStruct argRef;
  arg = argRef; 
  arg.shortId = ""; 
  arg.longId = "ref"; 
  arg.description = "Name of reference sequence";
  arg.required = false; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false; 
  ArgList.push_back(arg);
  ValueArg<string> cmd_ref(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // pos
  ArgStruct argPos;
  arg = argPos;
  arg.shortId = "";
  arg.longId = "pos";
  arg.description = "Reference position";
  arg.required = false;
  arg.defaultValueString = "1";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_pos(arg.shortId, arg.longId, arg.description, arg.required, 1, arg.type, cmd);

  // sample: naming scheme for matching reads to samples
  ArgStruct argSample;
  arg = argSample;
  arg.shortId = "";
  arg.longId = "sample";
  arg.description = "Naming scheme for matching reads to samples";
  arg.required = false;
  arg.defaultValueString = "multiple";
  arg.type = "string";
  arg.multi = false;
  vector<string> allowedSample;
  allowedSample.push_back("single");
  allowedSample.push_back("multiple");
  allowedSample.push_back("trio");
  allowedSample.push_back("unknown");
  arg.constraint = allowedSample;
  ValuesConstraint<string> allowedSampleVals(allowedSample); 
  ArgList.push_back(arg);
  ValueArg<string> cmd_sample(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, &allowedSampleVals, cmd);

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

  // an
  ArgStruct argAn;
  arg = argAn;
  arg.shortId = "";
  arg.longId = "an";
  arg.description = "Non-reference allele";
  arg.required = false;
  arg.defaultValueString = "A";
  arg.type = "string";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<string> cmd_an(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // prob
  ArgStruct argProb;
  arg = argProb;
  arg.shortId = "";
  arg.longId = "prob";
  arg.description = "Calculate SNP probability?";
  arg.required = false;
  arg.defaultValueString = "false";
  arg.type = "switch";
  arg.multi = false;
  ArgList.push_back(arg);
  SwitchArg cmd_prob(arg.shortId, arg.longId, arg.description, cmd, false);

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

  // QR: reference allele quality value
  ArgStruct argQR;
  arg = argQR;
  arg.shortId = "";
  arg.longId = "QR";
  arg.description = "Reference allele quality value (if = 0 value from GIG file used)";
  arg.required = false;
  arg.defaultValueString = "60";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_QR(arg.shortId, arg.longId, arg.description, arg.required, 60, arg.type, cmd);

  // ploidy: sample ploidy
  ArgStruct argPloidy;
  arg = argPloidy;
  arg.shortId = "";
  arg.longId = "ploidy";
  arg.description = "Sample ploidy";
  arg.required = false;
  arg.defaultValueString = "haploid";
  arg.type = "string";
  arg.multi = false;
  vector<string> allowedPloidy;
  allowedPloidy.push_back("haploid");
  allowedPloidy.push_back("diploid");
  arg.constraint = allowedPloidy;
  ValuesConstraint<string> allowedPloidyVals(allowedPloidy); 
  ArgList.push_back(arg);
  ValueArg<string> cmd_ploidy(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, &allowedPloidyVals, cmd);

  // QRL: minimum read base quality value
  ArgStruct argQRL;
  arg = argQRL;
  arg.shortId = "";
  arg.longId = "QRL";
  arg.description = "Minimum read base quality value";
  arg.required = false;
  arg.defaultValueString = "10";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_QRL(arg.shortId, arg.longId, arg.description, arg.required, 10, arg.type, cmd);

  // MRU: maximum number of mismatches between read and reference sequence
  ArgStruct argMRU;
  arg = argMRU;
  arg.shortId = "";
  arg.longId = "MRU";
  arg.description = "Maximum number of mismatches between read and refrence sequence";
  arg.required = false;
  arg.defaultValueString = "10000000";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_MRU(arg.shortId, arg.longId, arg.description, arg.required, 10000000, arg.type, cmd);

  // QML: minimum mismatch base quality value
  ArgStruct argQML;
  arg = argQML;
  arg.shortId = "";
  arg.longId = "QML";
  arg.description = "Minimum mismatch base quality value";
  arg.required = false;
  arg.defaultValueString = "10";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_QML(arg.shortId, arg.longId, arg.description, arg.required, 10, arg.type, cmd);

  // D: pairwise nucleotide diversity (theta)
  ArgStruct argD;
  arg = argD;
  arg.shortId = "";
  arg.longId = "D";
  arg.description = "Pairwise nucleotide diversity (theta)";
  arg.required = false;
  arg.defaultValueString = "10E-3";
  arg.type = "double";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<double> cmd_D(arg.shortId, arg.longId, arg.description, arg.required, 1E-3, arg.type, cmd);

  // M: per-generation somatic mutation rate (mu)
  ArgStruct argM;
  arg = argM;
  arg.shortId = "";
  arg.longId = "M";
  arg.description = "Per-generation somatic mutation rate (mu)";
  arg.required = false;
  arg.defaultValueString = "1E-8";
  arg.type = "double";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<double> cmd_M(arg.shortId, arg.longId, arg.description, arg.required, 1E-8, arg.type, cmd);

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
  string mbr = cmd_mbr.getValue();
  string ref = cmd_ref.getValue();
  int pos = cmd_pos.getValue();
  string sample = cmd_sample.getValue();
  string sampleDel = cmd_sampleDel.getValue();

  bool prob = cmd_prob.getValue();
  bool useRefAllele = cmd_useRefAllele.getValue();
  string an = cmd_an.getValue();
  int QR = cmd_QR.getValue();
  string ploidy = cmd_ploidy.getValue();
  long double D = (long double)cmd_D.getValue();
  long double M = (long double)cmd_M.getValue();
  int QRL = cmd_QRL.getValue();
  int MRU = cmd_MRU.getValue();
  int QML = cmd_QML.getValue();
  string algorithm = cmd_algorithm.getValue();
  int WB = cmd_WB.getValue();
  int TB = cmd_TB.getValue();
  int TR = cmd_TR.getValue();

  bool debug = cmd_debug.getValue();
  bool debug2 = cmd_debug2.getValue();

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

  bool diploid = false;
  if (ploidy == "diploid") {
    diploid = true;
  }

  bool banded = false;
  if (algorithm == "banded") {
    banded = true;
  }

  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  // constants
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  
  map<string, string, less<string> > IUPAC;
  IUPAC["A"] = "A"; 
  IUPAC["C"] = "C"; 
  IUPAC["G"] = "G"; 
  IUPAC["T"] = "T"; 
  IUPAC["AA"] = "A"; 
  IUPAC["AC"] = "M"; 
  IUPAC["AG"] = "R"; 
  IUPAC["AT"] = "W"; 
  IUPAC["CA"] = "M"; 
  IUPAC["CC"] = "C"; 
  IUPAC["CG"] = "S"; 
  IUPAC["CT"] = "Y"; 
  IUPAC["GA"] = "R"; 
  IUPAC["GC"] = "S"; 
  IUPAC["GG"] = "G"; 
  IUPAC["GT"] = "K"; 
  IUPAC["TA"] = "W"; 
  IUPAC["TC"] = "Y"; 
  IUPAC["TG"] = "K"; 
  IUPAC["TT"] = "T"; 

  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  // open MBR input file
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  
  // report
  if (debug) {cerr << "Opening MOSAIK binary format reference sequence archive input file: " << mbr << " ...";}
  
  // if this is not a valid MOSAIK reference file, complain and bomb
  if (! Mosaik::CReferenceSequenceReader::CheckFile(mbr, false)) {
    cerr << "ERROR: MOSAIK reference file not valid: " << mbr << ". Exiting..." << endl;
    exit(1);
  }

  // open the MOSAIK alignments file
  Mosaik::CReferenceSequenceReader rsr;
  rsr.Open(mbr);

  // retrieve all of the reference sequence metadata
  vector<Mosaik::ReferenceSequence> referenceSequences;
  rsr.GetReferenceSequences(referenceSequences);

  // load ref seq names into hash
  map<string, unsigned int, less<string> > refseqLength;
  for (vector<Mosaik::ReferenceSequence>::const_iterator rsIter = referenceSequences.begin();
       rsIter != referenceSequences.end(); ++rsIter) {
    refseqLength[rsIter->Name] = rsIter->NumBases;
  }
  
  // if ref is not valid bomb
  if (refseqLength.count(ref) <= 0) {
    cerr << "ERROR: ref does not exist: " << ref << endl;
    exit(1);
  }

  // if pos is not valid, bomb
  if (pos < 0 || pos > refseqLength[ref]) {
    cerr << "ERROR: pos out of bounds: " << pos << endl;
    exit(1);
  }

  // get the reference base sequence
  string bases;
  rsr.GetReferenceSequence(ref, bases);

  // close the reference sequence reader
  rsr.Close();
  
  // report success
  if (debug) {cerr << " done." << endl;}
  
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  // open BAM input file
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  
  // report
  if (debug) {cerr << "Opening BAM format alignment input file: " << bam << " ...";}
  
  // open
  const char * bamFileNameC = bam.c_str();
  string bamIndexFileName = bam + ".bai";
  const char * bamIndexFileNameC = bamIndexFileName.c_str();
 
  BamReader bReader;

  bReader.Open(bam, bamIndexFileName);
  
  // report
  if (debug) {cerr << " done." << endl;}
  
  //--------------------------------------------------------------------------
  // read reference sequences from input file
  //--------------------------------------------------------------------------
  
  // retrieve header information
  string bamHeader = bReader.GetHeaderText();
  
  // store the names of all the reference sequences in the BAM file
  RefVector refDatas = bReader.GetReferenceData();
  
  // data structs
  map<string, int, less<string> > RefLength, RefId;
  for (RefVector::const_iterator refIter = refDatas.begin();
       refIter != refDatas.end(); refIter++) {
    RefData refData = *refIter;
    string refName = refData.RefName;
    RefLength[refName] = refData.RefLength;
    RefId[refName] = bReader.GetReferenceID(refName);
  }
  
  if (debug) {
    cerr << "Number of ref seqs: " << bReader.GetReferenceCount() << endl;
  }
  
  //-------------------------------------------------------------------------
  //-------------------------------------------------------------------------
  // check for errors
  //-------------------------------------------------------------------------
  //-------------------------------------------------------------------------

  //-------------------------------------------------------------------------
  // ref does not exist
  //-------------------------------------------------------------------------
  if (RefId.count(ref) < 0) {
    cerr << "ERROR: reference sequence does not exist in BAM file: " << ref << ". Exiting..." << endl;
    exit(1);
  }

  //-------------------------------------------------------------------------
  // assign data ref data structure
  //-------------------------------------------------------------------------
  int refId = RefId[ref];
  RefData refData = refDatas[refId];
  int refLength = refData.RefLength;
      
  //-------------------------------------------------------------------------
  // ref does not have any alignments
  //-------------------------------------------------------------------------
  if (! refData.RefHasAlignments) {
    cerr << "WARNING: reference sequence does not have any alignments in BAM file: " << ref << ". Exiting..." << endl;
    exit(2);
  }

  //-------------------------------------------------------------------------
  //-------------------------------------------------------------------------
  // get ref pos base
  //-------------------------------------------------------------------------
  //-------------------------------------------------------------------------
  string ar = bases.substr(pos-1, 1);

  //-------------------------------------------------------------------------
  //-------------------------------------------------------------------------
  // data structures
  //-------------------------------------------------------------------------
  //-------------------------------------------------------------------------

  // contig-position specific coverage quantities
  map<int, int, less<int> > readCount, readCountPlus, readCountMinus;
  map<int, map_string_int, less<int> > alleleCount, alleleCountPlus, alleleCountMinus;
  map<int, map_string_longDouble, less<int> > alleleQual, alleleQualPlus, alleleQualMinus;
  
  // contig-position and individual specific coverage quantities
  map<int, map<string, vector<Basecall>, less<string> >, less<int> > individualBasecalls;
  
  // empty Basecall verctor
  vector<Basecall> emptyBasecalls;
  
  //-------------------------------------------------------------------------
  //-------------------------------------------------------------------------
  // jump to target
  //-------------------------------------------------------------------------
  //-------------------------------------------------------------------------

  //-------------------------------------------------------------------------
  // jumping to required position
  //-------------------------------------------------------------------------
  if (debug) {
    cerr << "Jumping to ref: " << ref << " position: " << pos  << endl;
  }
  if (! bReader.Jump(refId, pos)) {
    cerr << "WARNING: unable to jump to ref: " << ref << " position: " << pos << ". Exiting..." << endl;
    exit(2);
  }      
  //-------------------------------------------------------------------------
  // iterate thru every read overlapping pos
  //-------------------------------------------------------------------------

  if (debug) {
    cerr << "Iterating through reads aligned at this position." << endl;
  }

  BamAlignment ba;
  while (bReader.GetNextAlignment(ba) && (ba.Position + 1 <= pos)) {
    
    //-----------------------------------------------------------------------
    // checks
    //-----------------------------------------------------------------------

    // only process if mapped
    if (! ba.IsMapped()) {
      continue;
    }
	  
    // only process if not duplicate
    bool dupRead = false;
    if (ba.IsDuplicate()) {
      continue;
    }
	  
    //-----------------------------------------------------------------------
    // extract sample info
    //-----------------------------------------------------------------------
    string readName = ba.Name;
    SampleInfo sampleInfo = extractSampleInfo(readName, sample, sampleDel);
    string sampleName = sampleInfo.sampleId;
	  
    //-----------------------------------------------------------------------
    // parse cigar
    //-----------------------------------------------------------------------
    string rDna = ba.QueryBases;
    string rQual = ba.Qualities;

    /*  
    //-----------------------------------------------------------------------
    // filter on high quality mismatches between read and reference
    // make stats on duplicate reads
    //-----------------------------------------------------------------------
    int mru = 0;
    
    // initialize reference sequence position and read position
    // !!! +1 should be checked!!!
    unsigned int sp = ba.Position + 1;
    unsigned int rp = 1;
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
	  string b = rDna.substr(rp-1, 1);
	  char Q = rQual[rp-1];
	  
	  // convert base quality value into short int
	  short q = static_cast<short>(Q) - 33;
	  
	  // get reference allele
	  string sb = RefFastaData[refName].sequence.substr(sp-1, 1);	      
	  if (b != sb && q >= QML) {mru++;}
	  
	  // update positions
	  sp++;
	  rp++;
	}
      }
      else if (t == 'D' || t == 'N') { // deletion or skipped region
	for (int i=1; i<=l; i++) {
	  
	  // assign aligned base as gap
	  string b = "-";
	  char Q = '!';
	  
	  // convert base quality value into short int
	  short q = static_cast<short>(Q) - 33;
	  
	  // update ref position
	  sp++;
	}
      }
      else if (t == 'I') { // insertion
	rp += l;	    
      }
    }
    
    //----------------------------------------------------------------
    // skip read if too many mismatches
    //----------------------------------------------------------------
    if (mru > MRU) {
      continue;
    }
    */

    //----------------------------------------------------------------
    // register surviving alignment in coverage structures
    //----------------------------------------------------------------
	  
    // initialize reference sequence position and read position
    // !!! +1 should be checked!!!
    unsigned int sp = ba.Position + 1;
    unsigned int rp = 1;
    
    // iterate through Cigar operations and register read alignment
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
	  string b = rDna.substr(rp-1, 1);
	  char Q = rQual[rp-1];
	  
	  // convert base quality value into short int
	  short q = static_cast<short>(Q) - 33;
	  
	  // make base call data struct
	  Basecall bc;
	  bc.seqName = sampleName;
	  bc.base = b;
	  bc.qual = q;
	  bc.strand = "+";
	  if (ba.IsReverseStrand()) {
	    bc.strand = "-";
	  }
	  
	  // register base call if quality score sufficiently high
	  if (individualBasecalls[sp].count(sampleName) <= 0) {
	    individualBasecalls[sp][sampleName] = emptyBasecalls;
	  }
	  if (q >= QRL) {
	    individualBasecalls[sp][sampleName].push_back(bc);
	  }

	  // update positions
	  sp++;
	  rp++;
	}
      }
      else if (t == 'D' || t == 'N') { // deletion or skipped region
	for (int i=1; i<=l; i++) {
	  
	  // assign aligned base as gap
	  string b = "-";
	  char Q = '!';
	  
	  // convert base quality value into short int
	  short q = static_cast<short>(Q) - 33;
	  
	  // make base call data struct
	  Basecall bc;
	  bc.seqName = sampleName;
	  bc.base = b;
	  bc.qual = q;
	  bc.strand = "+";
	  if (ba.IsReverseStrand()) {
	    bc.strand = "-";
	  }
	  
	  // register base call
	  if (individualBasecalls[sp].count(sampleName) <= 0) {
	    individualBasecalls[sp][sampleName] = emptyBasecalls;
	  }
	  if (q >= QRL) {
	    individualBasecalls[sp][sampleName].push_back(bc);
	  }

	  // update ref position
	  sp++;
	}
      }
      else if (t == 'I') { // insertion
	rp += l;	    
      }
    }
  }

  //-------------------------------------------------------------------------
  //-------------------------------------------------------------------------
  // calculate genotype and snp probabilities
  //-------------------------------------------------------------------------
  //-------------------------------------------------------------------------

  //-------------------------------------------------------------------------
  // add ref allele to basecalls if needed
  //-------------------------------------------------------------------------
  if (useRefAllele) {

    if (debug) {
      cerr << "Adding reference allele." << endl;
    }

    // manufacture base call data structure for reference allele
    Basecall bc;
    bc.seqName = ref;
    bc.strand = "?";
    bc.base = ar;
    bc.qual = QR;
    
    // add to basecalls
    individualBasecalls[pos][ref].push_back(bc);
  }

  //-------------------------------------------------------------------------
  // calculate genotype data likelihoods
  //-------------------------------------------------------------------------

  if (debug) {
    cerr << "Calculating genotype likelihoods." << endl;
  }

  map<string, map<string, long double, less<string> >, less<string> > LnProbBiGivenGi;
  vector<string> sampleList;
  for (map<string, vector<Basecall>, less<string> >::const_iterator 
	 sampleIter = individualBasecalls[pos].begin();
       sampleIter != individualBasecalls[pos].end(); sampleIter++) {
    
    // retrieve
    string sample = sampleIter->first;
    vector<Basecall> basecalls = sampleIter->second;
    
    // record sample in list
    sampleList.push_back(sample);
    
    // calculate genotype log likelihoods
    map<string, long double, less<string> > logGl 
      = logGenotypeLikelihoods(
			       basecalls,
			       diploid,
			       0.9,
			       debug2
			       );
    
    // assign genotype log likelihoods in multi-individual data structure
    LnProbBiGivenGi[sample] = logGl;
  }	      
  
  //-------------------------------------------------------------------------
  // run posterior calculation for two best alleles
  //-------------------------------------------------------------------------
  if (debug) {
    cerr << "Calculating SNP and genotype probabilities." << endl;
  }

  Variation var = posteriorProb(
				LnProbBiGivenGi,
				sampleList,
				diploid,
				D,
				banded,
				WB,
				TB,
				false,
				TR,
				ar,
				an,
				debug2
				);

  //-------------------------------------------------------------------------
  //-------------------------------------------------------------------------
  // report
  //-------------------------------------------------------------------------
  //-------------------------------------------------------------------------

  //-------------------------------------------------------------------------
  // print reference base
  //-------------------------------------------------------------------------
  cout << "REF base at ref=" << ref << " pos=" << pos << " is: " << ar << endl;

  //-------------------------------------------------------------------------
  // print genotype probabilities and aligned read bases
  // calculate estimated allele frequency
  //-------------------------------------------------------------------------
  unsigned int numberAr = 0, numberAn = 0;
  for (map<string, vector<Basecall>, less<string> >::const_iterator 
	 sampleIter = individualBasecalls[pos].begin();
       sampleIter != individualBasecalls[pos].end(); sampleIter++) {

    string sample = sampleIter->first;
    vector<Basecall> basecalls = sampleIter->second;

    cout << "sample=" << sample << " #basecalls=" << basecalls.size() << " P(G):";
    string gBest = "";
    long double gBestP = 0.0;
    for (map<string, long double, less<string> >::const_iterator 
	   gIter = var.individualGenotypeProbability[sample].begin();
	 gIter != var.individualGenotypeProbability[sample].end();
	 gIter++) {
      string g = gIter->first;
      long double p = gIter->second;
      long double l = LnProbBiGivenGi[sample][IUPAC[g]];
      cout << " " << g << ": l=" << l << " p=" << p;

      // update best genotype
      if (p > gBestP) {
	gBestP = p;
	gBest = g;
      }
    }
    cout << endl;

    for (vector<Basecall>::const_iterator bcIter = basecalls.begin();
	 bcIter != basecalls.end(); bcIter++) {
      Basecall bc = *bcIter;
      cout << "  " << bc.base << " " << bc.qual << " " << bc.strand << endl;
    }

    // update allele frequency
    if (gBest == ar+ar) {
      numberAr += 2;
    }
    else if (gBest == an+an) {
      numberAn += 2;
    }
    else if (gBest == ar+an || gBest == an+ar) {
      numberAr++;
      numberAn++;
    }
  }
  
  //-------------------------------------------------------------------------
  // report P(SNP)
  //-------------------------------------------------------------------------
  cout << "P(SNP)=" << var.pSnp << endl;
  cout << "#ref alleles = " << numberAr << " #non-ref alleles = " << numberAn << endl;
  
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // finish up
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // report
  if (debug) {
    cerr << "Program completed." << endl;
  }

  // close BamReader
  bReader.Close();
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

