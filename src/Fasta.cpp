// ***************************************************************************
// FastaIndex.cpp (c) 2010 Erik Garrison <erik.garrison@bc.edu>
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 9 February 2010 (EG)
// ---------------------------------------------------------------------------

#include "Fasta.h"

FastaIndexEntry::FastaIndexEntry(string name, int length, long long offset, int line_blen, int line_len)
    : name(name)
    , length(length)
    , offset(offset)
    , line_blen(line_blen)
    , line_len(line_len)
{}

FastaIndexEntry::FastaIndexEntry(void) // empty constructor
{ clear(); }

FastaIndexEntry::~FastaIndexEntry(void)
{}

void FastaIndexEntry::clear(void)
{
    name = "";
    length = NULL;
    offset = -1;  // no real offset will ever be below 0, so this allows us to
                  // check if we have already recorded a real offset
    line_blen = NULL;
    line_len = NULL;
}

ostream& operator<<(ostream& output, const FastaIndexEntry& e) {
    output << e.name << "\t" << e.length << "\t" << e.offset << "\t" <<
        e.line_blen << "\t" << e.line_len;
    return output;  // for multiple << operators.
}

FastaIndex::FastaIndex(void) 
{}

void FastaIndex::readIndexFile(string fname) {
    string line;
    long long linenum = 0;
    vector<string> fields;
    indexFile.open(fname.c_str(), ifstream::in);
    if (indexFile.is_open()) {
        while (getline (indexFile, line)) {
            ++linenum;
            fields.clear();
            // the fai format defined in samtools is tab-delimited, every line being:
            // fai->name[i], (int)x.len, (long long)x.offset, (int)x.line_blen, (int)x.line_len
            boost::split(fields, line, boost::is_any_of("\t")); // split on tab
            if (fields.size() == 5) {  // if we don't get enough fields then there is a problem with the file
                // note that fields[0] is the sequence name
                this->push_back(FastaIndexEntry(fields[0], lexical_cast<int>(fields[1]),
                                                    lexical_cast<long long>(fields[2]),
                                                    lexical_cast<int>(fields[3]),
                                                    lexical_cast<int>(fields[4])));
            } else {
                cerr << "Warning: malformed fasta index file " << fname << 
                    "does not have enough fields @ line " << linenum << endl;
                cerr << line << endl;
                exit(1);
            }
        }
    }
}

// for consistency this should be a class method
bool fastaIndexEntryCompare ( FastaIndexEntry a, FastaIndexEntry b) { return (a.offset<b.offset); }

ostream& operator<<(ostream& output, const FastaIndex& fastaIndex) {
    vector<FastaIndexEntry> sortedIndex;
    for(vector<FastaIndexEntry>::const_iterator it = fastaIndex.begin(); it != fastaIndex.end(); ++it)
    {
        sortedIndex.push_back(*it);
    }
    sort(sortedIndex.begin(), sortedIndex.end(), fastaIndexEntryCompare);
    for( vector<FastaIndexEntry>::iterator fit = sortedIndex.begin(); fit != sortedIndex.end(); ++fit) {
        output << *fit << endl;
    }
}

void FastaIndex::indexReference(string refname) {
    // overview:
    //  for line in the reference fasta file
    //  track byte offset from the start of the file
    //  if line is a fasta header, take the name and dump the last sequnece to the index
    //  if line is a sequence, add it to the current sequence
    cerr << "indexing fasta reference " << refname << endl;
    string line;
    FastaIndexEntry entry;  // an entry buffer used in processing
    entry.clear();
    int line_length = 0;
    long long offset = 0;  // byte offset from start of file
    long long line_number = 0; // current line number
    bool mismatchedLineLengths = false; // flag to indicate if our line length changes mid-file
                                        // this will be used to raise an error
                                        // if we have a line length change at
                                        // any line other than the last line in
                                        // the sequence
    bool emptyLine = false;  // flag to catch empty lines, which we allow for
                             // index generation only on the last line of the sequence
    ifstream refFile;
    refFile.open(refname.c_str());
    if (refFile.is_open()) {
        while (getline(refFile, line)) {
            ++line_number;
            line_length = line.length();
            if (line[0] == ';') {
                // fasta comment, skip
            } else if (line[0] == '+') {
                // fastq quality header
                getline(refFile, line);
                line_length = line.length();
                offset += line_length + 1;
                // get and don't handle the quality line
                getline(refFile, line);
                line_length = line.length();
            } else if (line[0] == '>' || line[0] == '@') { // fasta /fastq header
                // if we aren't on the first entry, push the last sequence into the index
                if (entry.name != "") {
                    mismatchedLineLengths = false; // reset line length error tracker for every new sequence
                    emptyLine = false;
                    flushEntryToIndex(entry);
                    entry.clear();
                }
                entry.name = line.substr(1, line_length - 1);
            } else { // we assume we have found a sequence line
                if (entry.offset == -1) // NB initially the offset is -1
                    entry.offset = offset;
                entry.length += line_length;
                if (entry.line_len) {
                    //entry.line_len = entry.line_len ? entry.line_len : line_length + 1;
                    if (mismatchedLineLengths || emptyLine) {
                        if (line_length == 0) {
                            emptyLine = true; // flag empty lines, raise error only if this is embedded in the sequence
                        } else {
                            if (emptyLine) {
                                cerr << "ERROR: embedded newline";
                            } else {
                                cerr << "ERROR: mismatched line lengths";
                            }
                            cerr << " at line " << line_number << " within sequence " << entry.name <<
                                endl << "File not suitable for fasta index generation." << endl;
                            exit(1);
                        }
                    }
                    // this flag is set here and checked on the next line
                    // because we may have reached the end of the sequence, in
                    // which case a mismatched line length is OK
                    if (entry.line_len != line_length + 1) {
                        mismatchedLineLengths = true;
                        if (line_length == 0) {
                            emptyLine = true; // flag empty lines, raise error only if this is embedded in the sequence
                        }
                    }
                } else {
                    entry.line_len = line_length + 1; // first line
                }
                entry.line_blen = entry.line_len - 1;
            }
            offset += line_length + 1;
        }
        // we've hit the end of the fasta file!
        // flush the last entry
        flushEntryToIndex(entry);
    } else {
        cerr << "could not open reference file " << refname << " for indexing!" << endl;
        exit(1);
    }
}

void FastaIndex::flushEntryToIndex(FastaIndexEntry& entry) {
    this->push_back(FastaIndexEntry(entry.name, entry.length,
                                entry.offset, entry.line_blen,
                                entry.line_len));

}

void FastaIndex::writeIndexFile(string fname) {
    cerr << "writing fasta index file " << fname << endl;
    ofstream file;
    file.open(fname.c_str()); 
    if (file.is_open()) {
        file << *this;
    } else { 
        cerr << "could not open index file " << fname << " for writing!" << endl;
        exit(1);
    }
}

FastaIndex::~FastaIndex(void) {
    indexFile.close();
}

FastaIndexEntry FastaIndex::entry(string name) {
    for (vector<FastaIndexEntry>::const_iterator it = this->begin(); it != this->end(); ++it) {
        if (it->name == name)
            return *it;
    }
}

string FastaIndex::indexFileExtension() { return ".fai"; }

FastaReference::FastaReference(string reffilename) {
    filename = reffilename;
    file = fopen(filename.c_str(), "r");
    index = new FastaIndex();
    struct stat stFileInfo; 
    string indexFileName = filename + index->indexFileExtension(); 
    // if we can find an index file, use it
    if(stat(indexFileName.c_str(), &stFileInfo) == 0) { 
        index->readIndexFile(indexFileName);
    } else { // otherwise, read the reference and generate the index file in the cwd
        cerr << "index file " << indexFileName << " not found, generating..." << endl;
        index->indexReference(filename);
        index->writeIndexFile(indexFileName);
    }
}

FastaReference::~FastaReference(void) {
    fclose(file);
    delete index;
}

string FastaReference::getSequence(string seqname) {
    FastaIndexEntry entry = index->entry(seqname);
    int newlines_in_sequence = entry.length / entry.line_blen;
    int total_length = newlines_in_sequence  + entry.length;
    char* seq = (char*) malloc (sizeof(char)*total_length);
    fseek64(file, entry.offset, SEEK_SET);
    fread(seq, sizeof(char), total_length, file);
    char* pbegin = seq;
    char* pend = seq + ((string)seq).size()/sizeof(char);
    pend = remove (pbegin, pend, '\n');
    // TODO check if this works for \r
    // pend = remove (pbegin, pend, '\r');
    string s = seq;
    s.resize((pend - pbegin)/sizeof(char));
    return s;
}

string FastaReference::sequenceNameStartingWith(string seqnameStart) {
    string teststr;
    string result = "";
    vector<string> fields;
    for(vector<FastaIndexEntry>::const_iterator it = index->begin(); it != index->end(); ++it)
    {
        boost::split(fields, it->name, boost::is_any_of("\t ")); // split on ws
        // check if the first field is the same as our startseq
        teststr = fields.at(0);
        if (teststr == seqnameStart) {
            if (result == "") {
                result = it->name;
            } else {
                cerr << " is not unique in fasta index" << endl;
                exit(1);
                //return "";
            }
        }
    }
    if (result != "") {
        return result;
    } else {
        cerr << "could not find sequence named " << seqnameStart << endl;
        exit(1);
        //return ""; // XXX returning an empty string is an unclear result; should raise an error!!!
    }
}

string FastaReference::getSubSequence(string seqname, int start, int length) {
    FastaIndexEntry entry = index->entry(seqname);
    if (start < 0 || length < 1) {
        cerr << "Error: cannot construct subsequence with negative offset or length < 1" << endl;
        exit(1);
    }
    // we have to handle newlines
    // approach: count newlines before start
    //           count newlines by end of read
    //             subtracting newlines before start find count of embedded newlines
    int newlines_before = start > 0 ? (start - 1) / entry.line_blen : 0;
    int newlines_by_end = (start + length - 1) / entry.line_blen;
    int newlines_inside = newlines_by_end - newlines_before;
    /*
    cout << "start: " << start << " length: " << length << endl;
    cout << "newlines before: " << newlines_before << endl;
    cout << "newlines inside: " << newlines_inside << endl;
    */
    int seqlen = length + newlines_inside;
    char* seq = (char*) malloc (sizeof(char)*seqlen + 1);
    fseek64(file, (off_t) (entry.offset + newlines_before + start), SEEK_SET);
    fread(seq, sizeof(char), (off_t) seqlen, file);
    seq[seqlen] = '\0';
    string s = seq;
    boost::erase_all(s, "\n");
    return s;
}

long unsigned int FastaReference::sequenceLength(string seqname) {
    FastaIndexEntry entry = index->entry(seqname);
    return entry.length;
}

