// standard includes
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <boost/algorithm/string.hpp>

using namespace std;

#ifndef _PARAMETERS_H
#define _PARAMETERS_H

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
// static variables
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// classes
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
class MyOutput;

// Encapsulates tracking and parsing of command line program parameters
class Parameters {

    friend ostream &operator<<(ostream &out, const Parameters &p);

public:

    // i/o parameters:
    string bam;
    vector<string> bams;
    string fasta;
    string targets;
    string samples;
    string log;

    // operation parameters
    bool useRefAllele;
    bool outputAlleles;
    bool forceRefAllele;
    int MQR;
    int BQR;
    string ploidy;
    string sampleNaming;
    string sampleDel;
    int MQL0;
    int BQL0;
    int MQL1;
    int BQL1;
    int BQL2;
    int RMU;
    int IDW;
    long double TH;
    long double PVL;
    string algorithm;
    long double RDF;
    int WB;
    int TB;
    bool includeMonoB;
    int TR;
    int I;
    bool debug;
    bool debug2;
    bool record;

    // functions
    Parameters(int argc, char** argv);

};

#endif
