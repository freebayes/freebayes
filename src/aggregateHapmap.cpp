//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// aggregateHapmap.cpp
// aggregates HapMap genotypes across multiple HapMap genotype files
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

// private libraries
#include "Function-Generic.h"
#include "Class-BedReader.h"
#include "ReferenceSequenceReader.h"

// uses
using namespace std; 
using namespace TCLAP; 

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

static string ProgramName("rptStat");
static string ProgramDescription("Program to create data coverage and SNP stats.");
static string ProgramVersion("0.5.2");
static string ProgramDate("2009-04-29");
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

  // white space
  boost::regex re("\\s+");

  // HapMap file header line
  boost::regex patternHapmapHeaderLine("^rs#");

  // HapMap file data line
  boost::regex patternHapmapDataLine("^(rs\\d+)\\s+(\\S)/(\\S)\\s+chr(\\S+)\\s+(\\d+)");

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // variables
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // Regex match
  boost::smatch match;

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
  string mbr = cmd_mbr.getValue();
  int I = cmd_I.getValue();
  bool debug = cmd_debug.getValue();
  bool debug2 = cmd_debug2.getValue();

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // report command line and parameters used
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  map<bool, string, less<bool> > bool2String;
  bool2String[false] = "false";
  bool2String[true] = "true";

  if (debug) {
    cerr << "Command line:";
    for (int i=0; i<argc; i++) {
      cerr << " " << argv[i];
    }
    cerr << endl;
    cerr << endl;
    cerr << "Complete list of parameter values:" << endl;
    cerr << "  --mbr = " << mbr << endl;
    cerr << "  --I = " << I << endl;
    cerr << "  --debug = " <<  bool2String[debug] << endl;
    cerr << "  --debug2 = " <<  bool2String[debug2] << endl;
    cerr << endl;
  }

  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  // process MBR input file
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  
  //--------------------------------------------------------------------------
  // report
  //--------------------------------------------------------------------------
  if (debug) {cerr << "Opening MOSAIK binary format reference sequence archive input file: " << mbr << " ...";}
  
  //--------------------------------------------------------------------------
  // if this is not a valid MOSAIK reference file, complain and bomb
  //--------------------------------------------------------------------------
  if (! Mosaik::CReferenceSequenceReader::CheckFile(mbr, false)) {
    cerr << "ERROR: MOSAIK reference file not valid: " << mbr << ". Exiting..." << endl;
    exit(1);
  }

  //--------------------------------------------------------------------------
  // open the MOSAIK alignments file
  //--------------------------------------------------------------------------
  Mosaik::CReferenceSequenceReader rsr;
  rsr.Open(mbr);

  //--------------------------------------------------------------------------
  // retrieve all of the reference sequence metadata
  //--------------------------------------------------------------------------
  vector<Mosaik::ReferenceSequence> referenceSequences;
  rsr.GetReferenceSequences(referenceSequences);

  //--------------------------------------------------------------------------
  // load ref seq names into hash
  //--------------------------------------------------------------------------
  map<string, unsigned int, less<string> > refseqLength;
  map<string, string, less<string> > refseqDna;
  for (vector<Mosaik::ReferenceSequence>::const_iterator rsIter = referenceSequences.begin();
       rsIter != referenceSequences.end(); ++rsIter) {
    Mosaik::ReferenceSequence rs = *rsIter;
    string ref = rs.Name;
    refseqLength[ref] = rs.NumBases;
    
    // get the reference base sequence
    string bases;
    rsr.GetReferenceSequence(ref, bases);

    refseqDna[ref] = bases;
  }


  // close the reference sequence reader
  rsr.Close();
  
  // report success
  if (debug) {cerr << " done." << endl;}
  
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // read and process concatenated population-specific HapMap file
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  map<string, map<int, map<string, string, less<string> >, less<int> >, less<string> > G;
  map<string, map<int, string, less<int> >, less<string> > RS, A1, A2;
  map<string, bool, less<string> > sampleExists;
  vector<string> samples, allSamples;

  if (debug) {cerr << "Reading concatenated population-specific HapMap genotype file: " << endl;}
 
  //----------------------------------------------------------------------------
  // iterate through entries
  //----------------------------------------------------------------------------
  string line;
  unsigned long int lc = 1;
  unsigned int numberFields = 0;
  while (getline(cin, line)) {
      
    //--------------------------------------------------------------------------
    // increment line counter
    //--------------------------------------------------------------------------
    ++lc;

    //--------------------------------------------------------------------------
    // give life sign
    //--------------------------------------------------------------------------
    if (debug && lc % I == 0){
      cerr << "  processing line: " << lc << endl;
    }

    //--------------------------------------------------------------------------
    // parse line
    //--------------------------------------------------------------------------
    
    // if file header
    if (boost::regex_search(line, match, patternHapmapHeaderLine)) {	

      //------------------------------------------------------------------------
      // tokenize fields
      //------------------------------------------------------------------------
      boost::sregex_token_iterator rj;
      boost::sregex_token_iterator ri(line.begin(), line.end(), re, -1);
      vector<string> fields;
      while(ri != rj) {
	fields.push_back(*ri++);
      }
      
      //------------------------------------------------------------------------
      // set expected number of fields for data lines
      //------------------------------------------------------------------------
      numberFields = fields.size();
	
      //------------------------------------------------------------------------
      // get and register sample names
      //------------------------------------------------------------------------
      samples.clear();
      for (int i=11; i<numberFields; ++i) {
	string sample = fields[i];
	samples.push_back(sample);
	if (sampleExists.count(sample) <= 0) {
	  allSamples.push_back(sample);
	  sampleExists[sample] = true;
	}
      }
      
      //------------------------------------------------------------------------
      // report number samples if needed
      //------------------------------------------------------------------------
      if (debug) {
	cerr << "  number of samples: " << samples.size() << endl;
      }
    }

    // otherwise data line
    else if (boost::regex_search(line, match, patternHapmapDataLine)) {

      //------------------------------------------------------------------------
      // parse chromosome id, position, rs, and alleles
      //------------------------------------------------------------------------
      string rs = match[1];
      string allele1 = match[2];
      string allele2 = match[3];
      string ref = match[4];
      int pos = string2Int(match[5]);

      //------------------------------------------------------------------------
      // store rs, and alleles
      //------------------------------------------------------------------------
      RS[ref][pos] = rs;
      A1[ref][pos] = allele1;
      A2[ref][pos] = allele2;

      //------------------------------------------------------------------------
      // tokenize fields
      //------------------------------------------------------------------------
      boost::sregex_token_iterator rj;
      boost::sregex_token_iterator ri(line.begin(), line.end(), re, -1);
      vector<string> fields;
      while(ri != rj) {
	fields.push_back(*ri++);
      }
      
      //------------------------------------------------------------------------
      // check number of fields
      //------------------------------------------------------------------------
      if (fields.size() != numberFields) {
	cerr << "HapMap data line has incorrect number of fields: " << fields.size() 
	     << " instead of the expected number: " << numberFields<< endl;
	continue;
      }
 
      //------------------------------------------------------------------------
      // parse and register genotypes
      //------------------------------------------------------------------------
      int sc = 0;
      for (int i=11; i<numberFields; ++i) {
	string g = fields[i];
	string sample = samples[sc];
	G[ref][pos][sample] = g;
	++sc;
      }
    }
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // produce output
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  if (debug) {cerr << "Printing aggregate HapMap genotype file" << endl;}

  //----------------------------------------------------------------------------
  // print header line
  //----------------------------------------------------------------------------

  if (debug) {cerr << "  Printing header" << endl;}

  // print fixed fields
  cout << "refId\tpos\trsId\trefAllele\tallele1\tallele2\t#samples";

  // print sample names
  for (vector<string>::const_iterator sampleIter = allSamples.begin();
       sampleIter != allSamples.end(); ++sampleIter) {
    string sample = *sampleIter;
    cout << "\t" << sample;
  }
  cout << endl;

  //----------------------------------------------------------------------------
  // print metadata and genotypes for each ref position
  //----------------------------------------------------------------------------
  if (debug) {cerr << "  Printing genotypes" << endl;}
  unsigned int numSamples = allSamples.size();
  for (map<string, map<int, map<string, string, less<string> >, less<int> >, less<string> >::const_iterator refIter = G.begin();
       refIter != G.end(); ++refIter) {
    string ref = refIter->first;

    for (map<int, map<string, string, less<string> > >::const_iterator posIter = G[ref].begin();
       posIter != G[ref].end(); ++posIter) {
      int pos = posIter->first;
     
      // get reference allele
      string refAllele = refseqDna[ref].substr(pos-1, 1);	      

      cout << ref 
	   << "\t" << pos 
	   << "\t" << RS[ref][pos] 
	   << "\t" << refAllele 
	   << "\t" << A1[ref][pos] 
	   << "\t" << A2[ref][pos]
	   << "\t" << numSamples;

      // iterate through every sample in allSamples
      for (vector<string>::const_iterator sampleIter = allSamples.begin();
	   sampleIter != allSamples.end(); ++sampleIter) {
	string sample = *sampleIter;
	string g = "NE";
	if (G[ref][pos].count(sample) > 0) {
	  g = G[ref][pos][sample];
	}

	cout << "\t" << g;
      }
      cout << endl;
    } 
  }
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

