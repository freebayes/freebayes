/*
    vcflib C++ library for parsing and manipulating VCF files

    Copyright © 2010-2023 Erik Garrison
    Copyright © 2020-2023 Pjotr Prins

    This software is published under the MIT License. See the LICENSE file.
*/

#include "vcflib/Variant.h"
#include "vcflib/cigar.hpp"
// #include <utility>
#include "vcflib/multichoose.h"
#include <SmithWatermanGotoh.h>
#include "vcflib/ssw_cpp.hpp"
#include <regex>
#include "vcflib/join.h"

namespace vcflib {

static char rev_arr [26] = {84, 66, 71, 68, 69, 70, 67, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 65,
                           85, 86, 87, 88, 89, 90};

std::string reverse_complement(const std::string& seq) {
    // The old implementation of this function forgot to null-terminate its
    // returned string. This implementation uses heavier-weight C++ stuff that
    // may be slower but should ensure that that doesn't happen again.

    if (seq.size() == 0) {
        return seq;
    }

    string ret;
    ret.reserve(seq.size());

    std::transform(seq.rbegin(), seq.rend(), std::back_inserter(ret), [](char in) -> char {
        bool lower_case = (in >= 'a' && in <= 'z');
        if (lower_case) {
            // Convert to upper case
            in -= 32;
        }
        if (in < 'A' || in > 'Z') {
            throw std::runtime_error("Out of range character " + std::to_string((uint8_t)in) + " in inverted sequence");
        }
        // Compute RC in terms of letter identity, and then lower-case if necessary.
        return rev_arr[((int) in) - 'A'] + (lower_case ? 32 : 0);
    });

    return ret;
}

std::string toUpper(const std::string& seq) {
    if (seq.size() == 0) {
        return seq;
    }

    string ret;
    ret.reserve(seq.size());

    std::transform(seq.begin(), seq.end(), std::back_inserter(ret), [](char in) -> char {
        // If it's lower-case, bring it down in value to upper-case.
        return (in >= 'a' && in <= 'z') ? (in - 32) : in;
    });

    return ret;
}


bool allATGCN(const string& s, bool allowLowerCase){
    if (allowLowerCase){
       for (const auto c : s){
            if (c != 'A' && c != 'a' &&
                c != 'C' && c != 'c' &&
                c != 'T' && c != 't' &&
                c != 'G' && c != 'g' &&
                c != 'N' && c != 'n'){
                    return false;
            }
        }
    }
    else{
        for (const auto c : s){
            if (c != 'A' && c != 'C' && c != 'T' && c != 'G' && c != 'N'){
                return false;
            }
        }

    }
    return true;
}


/*
  Main VCF record parser
*/

void Variant::parse(string& line, bool parseSamples) {
    // clean up potentially variable data structures because the record may get reused(!)
    infoOrderedKeys.clear();
    info.clear();
    infoFlags.clear();
    format.clear();
    alt.clear();
    alleles.clear();

    // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT [SAMPLE1 .. SAMPLEN]
    vector<string> fields = split(line, '\t');
    if (fields.size() < 7) {
        cerr << "broken VCF record (less than 7 fields)" << endl
             << "Input line: " << line << endl;
        exit(1);
    }

    sequenceName = fields.at(0);
    char* end; // dummy variable for strtoll
    position = strtoll(fields.at(1).c_str(), &end, 10);
    id = fields.at(2);
    ref = fields.at(3);
    alt = split(fields.at(4), ","); // a comma-separated list of alternate alleles

    // make a list of all (ref + alts) alleles, allele[0] = ref, alleles[1:] = alts
    // add the ref allele ([0]), resize for the alt alleles, and then add the alt alleles
    alleles.push_back(ref);
    alleles.resize(alt.size()+1);
    std::copy(alt.begin(), alt.end(), alleles.begin()+1);

    // set up reverse lookup of allele index
    altAlleleIndexes.clear();
    int n = 0;
    for (vector<string>::iterator a = alt.begin();
            a != alt.end(); ++a, ++n) {
        altAlleleIndexes[*a] = n;
    }

    convert(fields.at(5), quality);
    filter = fields.at(6);
    // Process the INFO fields
    if (fields.size() > 7) {
        vector<string> infofields = split(fields.at(7), ';');
        for (const auto& field: infofields) {
            if (field == ".") {
                continue;
            }
            vector<string> kv = split(field, '='); // note that field gets split in place
            const auto& key = kv.at(0);
            if (kv.size() == 2) {
                split(kv.at(1), ',', info[key]); // value gets split in place
                infoOrderedKeys.push_back(key);
            } else if (kv.size() == 1) {
                infoFlags[key] = true;
                infoOrderedKeys.push_back(key);
            }
            // malformed fields with double '=' are silently skipped
        }
    }
    // check if we have samples specified
    // and that we are supposed to parse them
    if (parseSamples && fields.size() > 8) {
        format = split(fields.at(8), ':');
        // if the format changed, we have to rebuild the samples
        if (fields.at(8) != lastFormat) {
            samples.clear();
            lastFormat = fields.at(8);
        }
        vector<string>::iterator sampleName = sampleNames.begin();
        vector<string>::iterator sample = fields.begin() + 9;
        for (; sample != fields.end() && sampleName != sampleNames.end();
                ++sample, ++sampleName) {
			string& name = *sampleName;

			vector<string> samplefields = split(*sample, ':');
        	vector<string>::iterator i = samplefields.begin();

        	for (const auto& f : format) {
        		if(i != samplefields.end()){
        			samples[name][f] = split(*i, ',');
        			++i;
        		}
        		else{
        			std::vector<string> missing;
        			missing.push_back(".");
        			samples[name][f] = missing;
        		}
        	}
                }

        if (sampleName != sampleNames.end()) {
            cerr << "error: more sample names in header than sample fields" << endl;
            cerr << "samples: " << join(sampleNames, " ") << endl;
            cerr << "line: " << line << endl;
            exit(1);
        }
        if (sample != fields.end()) {
            cerr << "error: more sample fields than samples listed in header" << endl;
            cerr << "samples: " << join(sampleNames, " ") << endl;
            cerr << "line: " << line << endl;
            cerr << *sample << endl;
            exit(1);
        }
    } else if (!parseSamples) {
        originalLine = line;
    }

    //return true; // we should be catching exceptions...
}

bool Variant::hasSVTags() const{
    bool found_svtype = !getSVTYPE().empty();
    bool found_len = this->info.find("SVLEN") != this->info.end() || this->info.find("END") != this->info.end() || this->info.find("SPAN") != this->info.end();

   return found_svtype && found_len;
}

  /*
According to the VCF spec the ALT field can be use to indicate 'imprecise' structural
variants.
   */

bool Variant::isSymbolicSV() const{

    bool found_svtype = !getSVTYPE().empty();

    bool ref_valid = allATGCN(this->ref);
    bool alts_valid = true;
    for (const auto& a : this->alt){
        if (!allATGCN(a)){
            alts_valid = false;
        }
    }

    return (!ref_valid || !alts_valid) && (found_svtype);
}

string Variant::getSVTYPE(int altpos) const{

    if (altpos > 0){
        // TODO: Implement multi-alt SVs
        return "";
    }


    if (this->info.find("SVTYPE") != this->info.end()){
        if (altpos >= this->info.at("SVTYPE").size()) {
            return "";
        }
        return this->info.at("SVTYPE")[altpos];
    }

    return "";
};



int Variant::getMaxReferencePos(){
    if (this->canonical && this->info.find("END") != this->info.end()) {
        // We are cannonicalized and must have a correct END

        int end = 0;
        for (const auto& s : this->info.at("END")){
            // Get the latest one defined.
            end = max(abs(stoi(s)), end);
        }
        // Convert to 0-based.
        return end - 1;

    }

    if (!this->isSymbolicSV()){
        // We don't necessarily have an END, but we don't need one
        return this->zeroBasedPosition() + this->ref.length() - 1;
    }

    if (this->canonicalizable()){
        // We aren't canonical, but we could be.
        if (this->info.find("END") != this->info.end()){
            // We have an END; blindly trust it
            int end = 0;
            for (const auto& s : this->info.at("END")){
                // Get the latest one defined.
                end = max(abs(stoi(s)), end);
            }
            // Convert to 0-based.
            return end - 1;

        }
        else if (this->info.find("SVLEN") != this->info.end()){
            // There's no endpoint, but we know an SVLEN.
            // A negative SVLEN means a deletion, so if we find one we can say we delete that much.
            int deleted = 0;
            for (const auto& s : this->info.at("SVLEN")){
                int alt_len = stoi(s);
                if (alt_len > 0){
                    // Not a deletion, so doesn't affect any ref bases
                    continue;
                }
                deleted = max(-alt_len, deleted);
            }

            // The anchoring base at POS gets added in (because it isn't
            // deleted) but then subtracted out (because we have to do that to
            // match non-SV deletions). For insertions, deleted is 0 and we
            // return 0-based POS. Inversions must have an END.
            return this->zeroBasedPosition() + deleted;
        }
        else{
            cerr << "Warning: insufficient length information for " << *this << endl;
            return -1;
        }
    }
    else {
        cerr << "Warning: can't get end of non-canonicalizeable variant " << *this << endl;
    }
    return -1;
}




// To canonicalize a variant, we need either both REF and ALT seqs filled in
// or SVTYPE and SVLEN or END or SPAN or SEQ sufficient to define the variant.
bool Variant::canonicalizable(){
    bool pre_canon = allATGCN(this->ref);

    for (auto& a : this->alt){
        if (!allATGCN(a)){
            pre_canon = false;
        }
    }

    if (pre_canon){
        // It came in in a fully specified way.
        // TODO: ideally, we'd check to make sure ref/alt lengths, svtypes, and ends line up right here.
        return true;
    }

    string svtype = getSVTYPE();

    if (svtype.empty()){
        // We have no SV type, so we can't interpret things.
        return false;
    }

    // Check the tags
    bool has_len = this->info.count("SVLEN") && !this->info.at("SVLEN").empty();
    bool has_seq = this->info.count("SEQ") && !this->info.at("SEQ").empty();
    bool has_span = this->info.count("SPAN") && !this->info.at("SPAN").empty();
    bool has_end = this->info.count("END") && !this->info.at("END").empty();


    if (svtype == "INS"){
        // Insertions need a SEQ, SVLEN, or SPAN
        return has_seq || has_len || has_span;
    }
    else if (svtype == "DEL"){
        // Deletions need an SVLEN, SPAN, or END
        return has_len || has_span || has_end;
    }
    else if (svtype == "INV"){
        // Inversions need a SPAN or END
        return has_span || has_end;
    }
    else{
        // Other SV types are unsupported
        // TODO: DUP
        return false;
    }
}

bool Variant::canonicalize(FastaReference& fasta_reference, vector<FastaReference*> insertions, bool place_seq, int min_size){

    // Nobody should call this without checking
    assert(canonicalizable());

    // Nobody should call this twice
    assert(!this->canonical);

    // Find where the inserted sequence can come from for insertions
    bool do_external_insertions = !insertions.empty();
    FastaReference* insertion_fasta;
    if (do_external_insertions){
        insertion_fasta = insertions[0];
    }

    bool ref_valid = allATGCN(ref);

    if (!ref_valid && !place_seq){
        // If the reference is invalid, and we aren't allowed to change the ref sequence,
        // we can't canonicalize the variant.
        return false;
    }

    // Check the alts to see if they are not symbolic
    vector<bool> alt_i_atgcn (alt.size());
    for (int i = 0; i < alt.size(); ++i){
        alt_i_atgcn[i] = allATGCN(alt[i]);
    }

    // Only allow single-alt variants
    bool single_alt = alt.size() == 1;
    if (!single_alt){
        // TODO: this will need to be remove before supporting multiple alleles
        cerr << "Warning: multiple ALT alleles not yet allowed for SVs" << endl;
        return false;
    }

    // Fill in the SV tags
    string svtype = getSVTYPE();
    bool has_len = this->info.count("SVLEN") && !this->info.at("SVLEN").empty();
    bool has_seq = this->info.count("SEQ") && !this->info.at("SEQ").empty();
    bool has_span = this->info.count("SPAN") && !this->info.at("SPAN").empty();
    bool has_end = this->info.count("END") && !this->info.at("END").empty();

    // Where is the end, or where should it be?
    long info_end = 0;
    if (has_end) {
        // Get the END from the tag
        info_end = stol(this->info.at("END")[0]);
    }
    else if(ref_valid && !place_seq) {
        // Get the END from the reference sequence, which is ready.
        info_end = this->position + this->ref.length() - 1;
    }
    else if ((svtype == "DEL" || svtype == "INV") && has_span) {
        // For deletions and inversions, we can get the END from the SPAN
        info_end = this->position + abs(stol(this->info.at("SPAN")[0]));
    }
    else if (svtype == "DEL" && has_len) {
        // For deletions, we can get the END from the SVLEN
        info_end = this->position + abs(stol(this->info.at("SVLEN")[0]));
    }
    else if (svtype == "INS"){
        // For insertions, END is just POS if not specified
        info_end = this->position;
    }
    else{
        cerr << "Warning: could not set END info " << *this << endl;
        return false;
    }

    // Commit back the END
    this->info["END"].resize(1);
    this->info["END"][0] = to_string(info_end);
    has_end = true;

    // What is the variant length change?
    // We store it as absolute value
    long info_len = 0;
    if (has_len){
        // Get the SVLEN from the tag
        info_len = abs(stol(this->info.at("SVLEN")[0]));
    }
    else if ((svtype == "INS" || svtype == "DEL") && has_span){
        info_len = abs(stol(this->info.at("SPAN")[0]));
    }
    else if (svtype == "DEL"){
        // We always have the end by now
        // Deletion ends give you length change
        info_len = info_end - this->position;
    }
    else if (svtype == "INV"){
        // Inversions have 0 length change unless otherwise specified.
        info_len = 0;
    }
    else if (svtype == "INS" && has_seq) {
        // Insertions can let us pick it up from the SEQ tag
        info_len = this->info.at("SEQ").at(0).size();
    }
    else{
        cerr << "Warning: could not set SVLEN info " << *this << endl;
        return false;
    }

    // Commit the SVLEN back
    if (svtype == "DEL"){
        // Should be saved as negative
        this->info["SVLEN"].resize(1);
        this->info["SVLEN"][0] = to_string(-info_len);
    }
    else{
        // Should be saved as positive
        this->info["SVLEN"].resize(1);
        this->info["SVLEN"][0] = to_string(info_len);
    }
    // Now the length change is known
    has_len = true;

    // We also compute a span
    long info_span = 0;
    if (has_span){
        // Use the specified span
        info_span = abs(stol(this->info.at("SVLEN")[0]));
    }
    else if (svtype == "INS" || svtype == "DEL"){
        // has_len is always true here
        // Insertions and deletions let us determine the span from the length change, unless they are complex.
        info_span = info_len;
    }
    else if (svtype == "INV"){
        // has_end is always true here
        // Inversion span is start to end
        info_span = info_end - this->position;
    }
    else{
        cerr << "Warning: could not set SPAN info " << *this << endl;
        return false;
    }

    // Commit the SPAN back
    this->info["SPAN"].resize(1);
    this->info["SPAN"][0] = to_string(info_span);
    // Now the span change is known
    has_span = true;

    if (info_end < this->position) {
        cerr << "Warning: SV END is before POS [canonicalize] " <<
        *this << endl << "END: " << info_end << "  " << "POS: " << this->position << endl;
        return false;
    }

    if (has_seq) {
        // Force the SEQ to upper case, if already present
        this->info["SEQ"].resize(1);
        this->info["SEQ"][0] = toUpper(this->info["SEQ"][0]);
    }

    // Set the other necessary SV Tags (SVTYPE, SEQ (if insertion))
    // Also check for agreement in the position tags
    if (svtype == "INS"){
        if (info_end != this->position){
            cerr << "Warning: insertion END and POS do not agree (complex insertions not canonicalizeable) [canonicalize] " <<
            *this << endl << "END: " << info_end << "  " << "POS: " << this->position << endl;

            if (info_end == this->position + info_len) {
                // We can probably guess what they meant here.
                cerr << "Warning: VCF writer incorrecty produced END = POS + SVLEN for an insertion. Fixing END to POS." << endl;
                info_end = this->position;
                this->info["END"][0] = to_string(info_end);
            } else {
                return false;
            }
        }

        if (info_len != info_span){
            cerr << "Warning: insertion SVLEN and SPAN do not agree (complex insertions not canonicalizeable) [canonicalize] " <<
            *this << endl << "SVLEN: " << info_len << "  " << "SPAN: " << info_span << endl;
            return false;
        }

        if (has_seq && allATGCN(this->info.at("SEQ")[0]) && this->info.at("SEQ")[0].size() != info_len){
            cerr << "Warning: insertion SVLEN and SEQ do not agree (complex insertions not canonicalizeable) [canonicalize] " <<
            *this << endl << "SVLEN: " << info_len << "  " << "SEQ length: " << this->info.at("SEQ")[0].size() << endl;
            return false;
        }

        // Set REF
        string ref_base = toUpper(fasta_reference.getSubSequence(this->sequenceName, this->zeroBasedPosition(), 1));
        if (place_seq){
            this->ref.assign(ref_base);
        }

        if (has_seq &&
                 alt[0] != this->info.at("SEQ")[0] &&
                 allATGCN(this->info.at("SEQ")[0])){
            // Try to remove prepended ref sequence, assuming it's left-aligned
            string s = this->alt[0];
            s = toUpper(s.substr(this->ref.length()));
            if (s != this->info.at("SEQ")[0] && !place_seq){
                cerr << "Warning: INS sequence in alt field does not match SEQ tag" << endl <<
                this->alt[0] << " " << this->info.at("SEQ")[0] << endl;
                return false;
            }
            if (place_seq){
                this->alt[0].assign( ref_base + this->info.at("SEQ")[0] );
            }

        }
        else if (alt_i_atgcn[0] && !has_seq){
            string s = this->alt[0];
            s = toUpper(s.substr(this->ref.length()));
            this->info["SEQ"].resize(1);
            this->info.at("SEQ")[0].assign(s);

            if (s.size() != info_len){
                cerr << "Warning: insertion SVLEN and added bases do not agree (complex insertions not canonicalizeable) [canonicalize] " <<
                *this << endl << "SVLEN: " << info_len << "  " << "added bases: " << s.size() << endl;
                return false;
            }

        }
        else if (alt[0][0] == '<' && do_external_insertions){

            string ins_seq;
            string seq_id = alt[0].substr(1, alt[0].size() - 2);

            if (insertion_fasta->index->find(seq_id) != insertion_fasta->index->end()){
                ins_seq = toUpper(insertion_fasta->getSequence(seq_id));
                if (allATGCN(ins_seq)){
                    this->info["SEQ"].resize(1);
                    this->info["SEQ"][0].assign(ins_seq);
                    if (place_seq){
                        this->alt[0].assign(ref_base + ins_seq);
                    }
                }
                else {
                    cerr << "Warning: Loaded invalid alt sequence for: " << *this << endl;
                    return false;
                }

                if (ins_seq.size() != info_len){
                    cerr << "Warning: insertion SVLEN and FASTA do not agree (complex insertions not canonicalizeable) [canonicalize] " <<
                    *this << endl << "SVLEN: " << info_len << "  " << "FASTA bases: " << ins_seq.size() << endl;
                    return false;
                }
            }
            else{
                cerr << "Warning: could not locate alt sequence for: " << *this << endl;
                return false;
            }

        }
        else{
            cerr << "Warning: could not set SEQ [canonicalize]. " << *this << endl;
            return false;
        }
    }
    else if (svtype == "DEL"){
        // Note that info_len has been abs'd and is always positive
        if (this->position + info_len != info_end){
            cerr << "Warning: deletion END and SVLEN do not agree [canonicalize] " << *this << endl <<
            "END: " << info_end << "  " << "SVLEN: " << info_len << endl;
            return false;
        }

        if (this->position + info_span != info_end){
            cerr << "Warning: deletion END and SPAN do not agree [canonicalize] " << *this << endl <<
            "END: " << info_end << "  " << "SPAN: " << info_span << endl;
            return false;
        }

        if (info_end > fasta_reference.sequenceLength(this->sequenceName)) {
            cerr << "Warning: deletion END is past end of sequence [canonicalize] " << *this << endl <<
            "END: " << info_end << "  " << "length: " << fasta_reference.sequenceLength(this->sequenceName) << endl;
            return false;
        }

        // Set REF
        if (place_seq){
            string del_seq = toUpper(fasta_reference.getSubSequence(this->sequenceName, this->zeroBasedPosition(), info_len + 1));
            string ref_base = toUpper(fasta_reference.getSubSequence(this->sequenceName, this->zeroBasedPosition(), 1));
            this->ref.assign( del_seq );
            this->alt[0].assign( ref_base );
        }
    }
    else if (svtype == "INV"){
        if (this->position + info_span != info_end){
            cerr << "Warning: inversion END and SPAN do not agree [canonicalize] " << *this << endl <<
            "END: " << info_end << "  " << "SPAN: " << info_span << endl;
            return false;
        }

        if (info_len != 0){
            cerr << "Warning: inversion SVLEN specifies nonzero length change (complex inversions not canonicalizeable) [canonicalize] " <<
            *this << endl << "SVLEN: " << info_len << endl;

            if (info_end == this->position + info_len) {
                // We can probably guess what they meant here.
                cerr << "Warning: VCF writer incorrecty produced END = POS + SVLEN for an inversion. Fixing SVLEN to 0." << endl;
                info_len = 0;
                this->info["SVLEN"][0] = to_string(info_len);
            } else {
                return false;
            }
        }

        if (info_end > fasta_reference.sequenceLength(this->sequenceName)) {
            cerr << "Warning: inversion END is past end of sequence [canonicalize] " << *this << endl <<
            "END: " << info_end << "  " << "length: " << fasta_reference.sequenceLength(this->sequenceName) << endl;
            return false;
        }

        if (place_seq){
            string ref_seq = toUpper(fasta_reference.getSubSequence(this->sequenceName, this->zeroBasedPosition(), info_span + 1));
            // Note that inversions still need an anchoring left base at POS
            string inv_seq = ref_seq.substr(0, 1) + reverse_complement(ref_seq.substr(1));
            this->ref.assign(ref_seq);
            this->alt[0].assign(inv_seq);
        }

    }
    else{
        cerr << "Warning: invalid SV type [canonicalize]:" << *this << endl;
        return false;
    }


    this->updateAlleleIndexes();

    // Check for harmony between ref / alt / tags
    if (this->position > stol(this->info.at("END").at(0))){
        cerr << "Warning: position > END. Possible reference genome mismatch." << endl;
        return false;
    }

    if (svtype == "INS"){
        assert(!this->info.at("SEQ")[0].empty());
    }

    this->canonical = true;
    return true;
}

void Variant::setVariantCallFile(VariantCallFile& v) {
    sampleNames = v.sampleNames;
    outputSampleNames = v.sampleNames;
    vcf = &v;
}

void Variant::setVariantCallFile(VariantCallFile* v) {
    sampleNames = v->sampleNames;
    outputSampleNames = v->sampleNames;
    vcf = v;
}

ostream& operator<<(ostream& out, VariantFieldType type) {
    switch (type) {
        case FIELD_INTEGER:
            out << "integer";
            break;
        case FIELD_FLOAT:
            out << "float";
            break;
        case FIELD_BOOL:
            out << "bool";
            break;
        case FIELD_STRING:
            out << "string";
            break;
        default:
            out << "unknown";
            break;
    }
    return out;
}

VariantFieldType typeStrToVariantFieldType(string& typeStr) {
    if (typeStr == "Integer") {
        return FIELD_INTEGER;
    } else if (typeStr == "Float") {
        return FIELD_FLOAT;
    } else if (typeStr == "Flag") {
        return FIELD_BOOL;
    } else if (typeStr == "String") {
        return FIELD_STRING;
    } else {
        return FIELD_UNKNOWN;
    }
}

VariantFieldType Variant::infoType(const string& key) {
    map<string, VariantFieldType>::iterator s = vcf->infoTypes.find(key);
    if (s == vcf->infoTypes.end()) {
        if (key == "FILTER") { // hack to use FILTER as an "info" field (why the hack?)
            return FIELD_STRING;
        }
        if (key == "QUAL") { // hack to use QUAL as an "info" field
            return FIELD_INTEGER;
        }
        cerr << "no info field " << key << endl;
        exit(1);
    } else {
        return s->second;
    }
}

    VariantFieldType Variant::formatType(const string& key) {
        map<string, VariantFieldType>::iterator s = vcf->formatTypes.find(key);
        if (s == vcf->formatTypes.end()) {
            cerr << "no format field " << key << endl;
            exit(1);
        } else {
            return s->second;
        }
    }

    bool Variant::getInfoValueBool(const string& key, int index) {
        map<string, VariantFieldType>::iterator s = vcf->infoTypes.find(key);
        if (s == vcf->infoTypes.end()) {
            cerr << "no info field " << key << endl;
            exit(1);
        } else {
            int count = vcf->infoCounts[key];
            // XXX TODO, fix for Genotype variants...
            if (count != ALLELE_NUMBER) {
                index = 0;
            }
            if (index == INDEX_NONE) {
                if (count != 1) {
                    cerr << "no field index supplied and field count != 1" << endl;
                    exit(1);
                } else {
                    index = 0;
                }
            }
            VariantFieldType type = s->second;
            if (type == FIELD_BOOL) {
                map<string, bool>::iterator b = infoFlags.find(key);
                if (b == infoFlags.end())
                    return false;
                else
                    return true;
            } else {
                cerr << "not flag type " << key << endl;
                exit(1);
            }
        }
    }

    string Variant::getInfoValueString(const string& key, int index) {
        map<string, VariantFieldType>::iterator s = vcf->infoTypes.find(key);
        if (s == vcf->infoTypes.end()) {
            if (key == "FILTER") {
              return filter;
            }
            cerr << "no info field " << key << endl;
            exit(1);
        } else {
            int count = vcf->infoCounts[key];
            // XXX TODO, fix for Genotype variants...
            if (count != ALLELE_NUMBER) {
                index = 0;
            }
            if (index == INDEX_NONE) {
                if (count != 1) {
                    cerr << "no field index supplied and field count != 1" << endl;
                    exit(1);
                } else {
                    index = 0;
                }
            }
            VariantFieldType type = s->second;
            if (type == FIELD_STRING) {
                map<string, vector<string> >::iterator b = info.find(key);
                if (b == info.end())
                    return "";
                return b->second.at(index);
            } else {
                cerr << "not string type " << key << endl;
                return "";
            }
        }
    }

    double Variant::getInfoValueFloat(const string& key, int index) {
        map<string, VariantFieldType>::iterator s = vcf->infoTypes.find(key);
        if (s == vcf->infoTypes.end()) {
            if (key == "QUAL") {
                return quality;
            }
            cerr << "no info field " << key << endl;
            exit(1);
        } else {
            int count = vcf->infoCounts[key];
            // XXX TODO, fix for Genotype variants...
            if (count != ALLELE_NUMBER) {
                index = 0;
            }
            if (index == INDEX_NONE) {
                if (count != 1) {
                    cerr << "no field index supplied and field count != 1" << endl;
                    exit(1);
                } else {
                    index = 0;
                }
            }
            VariantFieldType type = s->second;
            if (type == FIELD_FLOAT || type == FIELD_INTEGER) {
                map<string, vector<string> >::iterator b = info.find(key);
                if (b == info.end())
                    return false;
                double r;
                if (!convert(b->second.at(index), r)) {
                    cerr << "could not convert field " << key << "=" << b->second.at(index) << " to " << type << endl;
                    exit(1);
                }
                return r;
            } else {
                cerr << "unsupported type for variant record " << type << endl;
                exit(1);
            }
        }
    }

    int Variant::getNumSamples(void) {
        return sampleNames.size();
    }

    int Variant::getNumValidGenotypes(void) {
        int valid_genotypes = 0;
        map<string, map<string, vector<string> > >::const_iterator s     = samples.begin();
        map<string, map<string, vector<string> > >::const_iterator sEnd  = samples.end();
        for (; s != sEnd; ++s) {
            map<string, vector<string> > sample_info = s->second;
            if (sample_info["GT"].front() != "./.") {
                valid_genotypes++;
            }
        }
        return valid_genotypes;
    }

    bool Variant::getSampleValueBool(const string& key, string& sample, int index) {
        map<string, VariantFieldType>::iterator s = vcf->formatTypes.find(key);
        if (s == vcf->infoTypes.end()) {
            cerr << "no info field " << key << endl;
            exit(1);
        } else {
            int count = vcf->formatCounts[key];
            // XXX TODO, fix for Genotype variants...
            if (count != ALLELE_NUMBER) {
                index = 0;
            }
            if (index == INDEX_NONE) {
                if (count != 1) {
                    cerr << "no field index supplied and field count != 1" << endl;
                    exit(1);
                } else {
                    index = 0;
                }
            }
            VariantFieldType type = s->second;
            map<string, vector<string> >& sampleData = samples[sample];
            if (type == FIELD_BOOL) {
                map<string, vector<string> >::iterator b = sampleData.find(key);
                if (b == sampleData.end())
                    return false;
                else
                    return true;
            } else {
                cerr << "not bool type " << key << endl;
                exit(1);
            }
        }
    }

    string Variant::getSampleValueString(const string& key, string& sample, int index) {
        map<string, VariantFieldType>::iterator s = vcf->formatTypes.find(key);
        if (s == vcf->infoTypes.end()) {
            cerr << "no info field " << key << endl;
            exit(1);
        } else {
            int count = vcf->formatCounts[key];
            // XXX TODO, fix for Genotype variants...
            if (count != ALLELE_NUMBER) {
                index = 0;
            }
            if (index == INDEX_NONE) {
                if (count != 1) {
                    cerr << "no field index supplied and field count != 1" << endl;
                    exit(1);
                } else {
                    index = 0;
                }
            }
            VariantFieldType type = s->second;
            map<string, vector<string> >& sampleData = samples[sample];
            if (type == FIELD_STRING) {
                map<string, vector<string> >::iterator b = sampleData.find(key);
                if (b == sampleData.end()) {
                    return "";
                } else {
                    return b->second.at(index);
                }
            } else {
                cerr << "not string type " << key << endl;
                exit(1);
            }
        }
    }

    double Variant::getSampleValueFloat(const string& key, string& sample, int index) {
        map<string, VariantFieldType>::iterator s = vcf->formatTypes.find(key);
        if (s == vcf->infoTypes.end()) {
            cerr << "no info field " << key << endl;
            exit(1);
        } else {
            // XXX TODO wrap this with a function call
            int count = vcf->formatCounts[key];
            // XXX TODO, fix for Genotype variants...
            if (count != ALLELE_NUMBER) {
                index = 0;
            }
            if (index == INDEX_NONE) {
                if (count != 1) {
                    cerr << "no field index supplied and field count != 1" << endl;
                    exit(1);
                } else {
                    index = 0;
                }
            }
            VariantFieldType type = s->second;
            map<string, vector<string> >& sampleData = samples[sample];
            if (type == FIELD_FLOAT || type == FIELD_INTEGER) {
                map<string, vector<string> >::iterator b = sampleData.find(key);
                if (b == sampleData.end())
                    return false;
                double r;
                if (!convert(b->second.at(index), r)) {
                    cerr << "could not convert field " << key << "=" << b->second.at(index) << " to " << type << endl;
                    exit(1);
                }
                return r;
            } else {
                cerr << "unsupported type for sample " << type << endl;
                exit(1);
            }
        }
    }

    bool Variant::getValueBool(const string& key, string& sample, int index) {
        if (sample.empty()) { // an empty sample name means
            return getInfoValueBool(key, index);
        } else {
            return getSampleValueBool(key, sample, index);
        }
    }

    double Variant::getValueFloat(const string& key, string& sample, int index) {
        if (sample.empty()) { // an empty sample name means
            return getInfoValueFloat(key, index);
        } else {
            return getSampleValueFloat(key, sample, index);
        }
    }

    string Variant::getValueString(const string& key, string& sample, int index) {
        if (sample.empty()) { // an empty sample name means
            return getInfoValueString(key, index);
        } else {
            return getSampleValueString(key, sample, index);
        }
    }

    int Variant::getAltAlleleIndex(const string& allele) {
        map<string, int>::iterator f = altAlleleIndexes.find(allele);
        if (f == altAlleleIndexes.end()) {
            cerr << "no such allele \'" << allele << "\' in record " << sequenceName << ":" << position << endl;
            exit(1);
        } else {
            return f->second;
        }
    }

    void Variant::addFilter(const string& tag) {
        if (filter == "" || filter == ".")
            filter = tag;
        else
            filter += "," + tag;
    }

    void Variant::addFormatField(const string& key) {
        bool hasTag = false;
        for (const auto& t : format) {
            if (t == key) {
                hasTag = true;
                break;
            }
        }
        if (!hasTag) {
            format.push_back(key);
        }
    }

    void Variant::printAlt(ostream& out) const {
        for (vector<string>::const_iterator i = alt.begin(); i != alt.end(); ++i) {
            out << *i;
            // add a comma for all but the last alternate allele
            if (i != (alt.end() - 1)) out << ",";
        }
    }

    void Variant::printAlleles(ostream& out) const {
        for (vector<string>::const_iterator i = alleles.begin(); i != alleles.end(); ++i) {
            out << *i;
            // add a comma for all but the last alternate allele
            if (i != (alleles.end() - 1)) out << ",";
        }
    }

    /*
      This is the main outputter of VCF records/lines
    */
    ostream& operator<<(ostream& out, Variant& var) {
        // ensure there are no empty fields
        if (var.sequenceName.empty()) var.sequenceName = ".";
        if (var.id.empty()) var.id = ".";
        if (var.ref.empty()) var.ref = ".";
        if (var.alt.empty()) var.alt.push_back(".");
        if (var.filter.empty()) var.filter = ".";

        out << var.sequenceName << "\t"
            << var.position << "\t"
            << var.id << "\t"
            << var.ref << "\t";
        // report the list of alternate alleles.
        var.printAlt(out);

        out << "\t"
            << var.quality << "\t"
            << var.filter << "\t";
        if (var.info.empty() && var.infoFlags.empty()) {
            out << ".";
        } else {
            // We want to display the info fields in the original
            // order.  Because the actual info list may have been
            // modified since the record was read, we need to recreate
            // a valid ordered key list.
            map<string,bool> lookup_keys; // for quick lookup in 2nd step
            vector<string> ordered_keys, missing_keys;  // the output list
            // first lookup the keys that appear both in infoOrdered keys
            // and the info field:
            for (const auto& name: var.infoOrderedKeys)
            {
                lookup_keys[name] = true;
                if (!var.info[name].empty()) ordered_keys.push_back(name);
                if (var.infoFlags[name]) ordered_keys.push_back(name);
            };
            // next add the keys that are not in the original list:
            for (const auto& [name1, value]: var.info)
                if (!lookup_keys[name1]) missing_keys.push_back(name1);
            for (const auto& [name2, value]: var.infoFlags)
                if (lookup_keys[name2] == false) missing_keys.push_back(name2);

            // append sorted missing keys
            std::sort(missing_keys.begin(), missing_keys.end());

            ordered_keys.insert(ordered_keys.end(), missing_keys.begin(), missing_keys.end());
            // output the ordered info fields
            string s = "";
            for (const auto& name: ordered_keys) {
                const auto& value = var.info[name];
                if (!value.empty()) {
                    s += name + "=" + join(value, ",") + ";" ;
                } else {
                    const auto infoflag = var.infoFlags[name];
                    if (infoflag == true)
                        s += name + ";";
                }
            }
            auto len = s.length();
            if (len)
                out << s.substr(0, len-1); // chop s1.substr(0, i-1);
        }
        if (!var.format.empty()) {
            out << "\t";
            string format = "";
            for (const auto& f: var.format) {
                format += f + ":";
            }
            auto len = format.length();
            if (len)
                out << format.substr(0, len-1); // chop s1.substr(0, i-1);
            for (const auto& s: var.outputSampleNames) {
                out << "\t";
                const auto sampleItr = var.samples.find(s);
                if (sampleItr == var.samples.end()) {
                    out << ".";
                } else {
                    const map<string, vector<string> >& sample = sampleItr->second;
                    if (sample.empty()) {
                        out << ".";
                    } else {
                        for (vector<string>::iterator f = var.format.begin(); f != var.format.end(); ++f) {
                            const auto g = sample.find(*f);
                            out << ((f == var.format.begin()) ? "" : ":");
                            if (g != sample.end() && !g->second.empty()) {
                                out << join(g->second, ",");
                            } else {
                                out << ".";
                            }
                        }
                    }
                }
            }
        }
        return out;
    }

    void Variant::setOutputSampleNames(const vector<string>& samplesToOutput) {
        outputSampleNames = samplesToOutput;
    }


// shunting yard algorithm
    void infixToPrefix(queue<RuleToken> tokens, queue<RuleToken>& prefixtokens) {
        stack<RuleToken> ops;
        while (!tokens.empty()) {
            RuleToken& token = tokens.front();
            if (isOperator(token)) {
                //cerr << "found operator " << token.value << endl;
                while (ops.size() > 0 && isOperator(ops.top())
                       && (   (isLeftAssociative(token)  && priority(token) <= priority(ops.top()))
                              || (isRightAssociative(token) && priority(token) <  priority(ops.top())))) {
                    prefixtokens.push(ops.top());
                    ops.pop();
                }
                ops.push(token);
            } else if (isLeftParenthesis(token)) {
                //cerr << "found paran " << token.value << endl;
                ops.push(token);
            } else if (isRightParenthesis(token)) {
                //cerr << "found paran " << token.value << endl;
                while (ops.size() > 0 && !isLeftParenthesis(ops.top())) {
                    prefixtokens.push(ops.top());
                    ops.pop();
                }
                if (ops.size() == 0) {
                    cerr << "error: mismatched parentheses" << endl;
                    exit(1);
                }
                if (isLeftParenthesis(ops.top())) {
                    ops.pop();
                }
            } else {
                //cerr << "found operand " << token.value << endl;
                prefixtokens.push(token);
            }
            tokens.pop();
        }
        while (ops.size() > 0) {
            if (isRightParenthesis(ops.top()) || isLeftParenthesis(ops.top())) {
                cerr << "error: mismatched parentheses" << endl;
                exit(1);
            }
            prefixtokens.push(ops.top());
            ops.pop();
        }
    }

    RuleToken::RuleToken(const string& tokenstr, map<string, VariantFieldType>& variables) {
        isVariable = false;
        if (tokenstr == "!") {
            type = RuleToken::NOT_OPERATOR;
        } else if (tokenstr == "&") {
            type = RuleToken::AND_OPERATOR;
        } else if (tokenstr == "|") {
            type = RuleToken::OR_OPERATOR;
        } else if (tokenstr == "+") {
            type = RuleToken::ADD_OPERATOR;
        } else if (tokenstr == "-") {
            type = RuleToken::SUBTRACT_OPERATOR;
        } else if (tokenstr == "*") {
            type = RuleToken::MULTIPLY_OPERATOR;
        } else if (tokenstr == "/") {
            type = RuleToken::DIVIDE_OPERATOR;
        } else if (tokenstr == "=") {
            type = RuleToken::EQUAL_OPERATOR;
        } else if (tokenstr == ">") {
            type = RuleToken::GREATER_THAN_OPERATOR;
        } else if (tokenstr == "<") {
            type = RuleToken::LESS_THAN_OPERATOR;
        } else if (tokenstr == "(") {
            type = RuleToken::LEFT_PARENTHESIS;
        } else if (tokenstr == ")") {
            type = RuleToken::RIGHT_PARENTHESIS;
        } else { // operand
            type = RuleToken::OPERAND;
            if (variables.find(tokenstr) == variables.end()) {
                if (convert(tokenstr, number)) {
                    type = RuleToken::NUMBER;
                } else if (tokenstr == "QUAL") {
                    isVariable = true;
                } else if (tokenstr == "FILTER") {
                    isVariable = true;
                } else {
                    type = RuleToken::STRING_VARIABLE;
                }
            } else {
                isVariable = true;
            }
        }
        value = tokenstr;
    }


    void tokenizeFilterSpec(string& filterspec, queue<RuleToken>& tokens, map<string, VariantFieldType>& variables) {
        string lastToken = "";
        bool inToken = false;
        for (unsigned int i = 0; i <  filterspec.size(); ++i) {
            char c = filterspec.at(i);
            if (c == ' ' || c == '\n') {
                inToken = false;
                if (!inToken && lastToken.size() > 0) {
                    tokens.push(RuleToken(lastToken, variables));
                    lastToken = "";
                }
            } else if (!inToken && (isOperatorChar(c) || isParanChar(c))) {
                inToken = false;
                if (lastToken.size() > 0) {
                    tokens.push(RuleToken(lastToken, variables));
                    lastToken = "";
                }
                tokens.push(RuleToken(filterspec.substr(i,1), variables));
            } else {
                inToken = true;
                lastToken += c;
            }
        }
        // get the last token
        if (inToken) {
            tokens.push(RuleToken(lastToken, variables));
        }
    }

// class which evaluates filter expressions
// allow filters to be defined using boolean infix expressions e.g.:
//
// "GQ > 10 & (DP < 3 | DP > 5) & SAMPLE = NA12878"
// or
// "GT = 1/1 | GT = 0/0"
//
// on initialization, tokenizes the input sequence, and converts it from infix to postfix
// on call to
//


    VariantFilter::VariantFilter(string filterspec, VariantFilterType filtertype, map<string, VariantFieldType>& variables) {
        type = filtertype;
        spec = filterspec;
        tokenizeFilterSpec(filterspec, tokens, variables);
        infixToPrefix(tokens, rules);
        /*while (!rules.empty()) {
          cerr << " " << rules.front().value << ((isNumeric(rules.front())) ? "f" : "");
          rules.pop();
          }
        */
        //cerr << endl;
        //cerr << join(" ", tokens) << endl;
    }

// all alts pass
    bool VariantFilter::passes(Variant& var, string& sample) {
        for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
            string& allele = *a;
            if (!passes(var, sample, allele)) {
                return false;
            }
        }
        return true;
    }

    bool VariantFilter::passes(Variant& var, string& sample, string& allele) {
        // to evaluate a rpn boolean queue with embedded numbers and variables
        // make a result stack, use float to allow comparison of floating point
        // numbers, booleans, and integers
        stack<RuleToken> results;
        queue<RuleToken> rulesCopy = rules; // copy

        int index;
        if (allele.empty()) {
            index = 0; // apply to the whole record
        } else {
            // apply to a specific allele
            index = var.getAltAlleleIndex(allele);
        }

        while (!rulesCopy.empty()) {
            RuleToken token = rulesCopy.front();
            rulesCopy.pop();
        // pop operands from the front of the queue and push them onto the stack
        if (isOperand(token)) {
            //cout << "is operand: " << token.value << endl;
            // if the token is variable, i.e. not evaluated in this context, we
            // must evaluate it before pushing it onto the stack
            if (token.isVariable) {
                //cout << "is variable" << endl;
                // look up the variable using the Variant, depending on our filter type
                //cout << "token.value " << token.value << endl;
                VariantFieldType vtype;
                if (sample.empty()) { // means we are record-specific
                    vtype = var.infoType(token.value);
                } else {
                    vtype = var.formatType(token.value);
                    //cout << "type = " << type << endl;
                }
                //cout << "type: " << type << endl;

                if (vtype == FIELD_INTEGER || vtype == FIELD_FLOAT) {
                    token.type = RuleToken::NUMERIC_VARIABLE;
                    token.number = var.getValueFloat(token.value, sample, index);
                    //cerr << "number: " << token.number << endl;
                } else if (vtype == FIELD_BOOL) {
                    token.type = RuleToken::BOOLEAN_VARIABLE;
                    token.state = var.getValueBool(token.value, sample, index);
                    //cerr << "state: " << token.state << endl;
                } else if (vtype == FIELD_STRING) {
                    //cout << "token.value = " << token.value << endl;
                    token.type = RuleToken::STRING_VARIABLE;
                    token.str = var.getValueString(token.value, sample, index);
                } else if (isString(token)) {
                    token.type = RuleToken::STRING_VARIABLE;
                    token.str = var.getValueString(token.value, sample, index);
                    //cerr << "string: " << token.str << endl;
                }
            } else {
                double f;
                string s;
                //cerr << "parsing operand" << endl;
                if (convert(token.value, f)) {
                    token.type = RuleToken::NUMERIC_VARIABLE;
                    token.number = f;
                    //cerr << "number: " << token.number << endl;
                } else if (convert(token.value, s)) {
                    token.type = RuleToken::STRING_VARIABLE;
                    token.str = s;
                    //cerr << "string: " << token.str << endl;
                } else {
                    cerr << "could not parse non-variable operand " << token.value << endl;
                    exit(1);
                }
            }
            results.push(token);
        }
        // apply operators to the first n elements on the stack and push the result back onto the stack
        else if (isOperator(token)) {
            //cerr << "is operator: " << token.value << endl;
            RuleToken a, b, r;
            // is it a not-operator?
            switch (token.type) {
                case ( RuleToken::NOT_OPERATOR ):
                    a = results.top();
                    results.pop();
                    if (!isBoolean(a)) {
                        cerr << "cannot negate a non-boolean" << endl;
                    } else {
                        a.state = !a.state;
                        results.push(a);
                    }
                    break;

                case ( RuleToken::EQUAL_OPERATOR ):
                    a = results.top(); results.pop();
                    b = results.top(); results.pop();
                    if (a.type == b.type) {
                        switch (a.type) {
                            case (RuleToken::STRING_VARIABLE):
                                r.state = (a.str == b.str);
                                break;
                            case (RuleToken::NUMERIC_VARIABLE):
                                r.state = (a.number == b.number);
                                break;
                            case (RuleToken::BOOLEAN_VARIABLE):
                                r.state = (a.state == b.state);
                                break;
                            default:
                                cerr << "should not get here" << endl; exit(1);
                                break;
                        }
                    } else if (a.type == RuleToken::STRING_VARIABLE && b.type == RuleToken::NUMERIC_VARIABLE) {
                        r.state = (convert(b.number) == a.str);
                    } else if (b.type == RuleToken::STRING_VARIABLE && a.type == RuleToken::NUMERIC_VARIABLE) {
                        r.state = (convert(a.number) == b.str);
                    }
                    results.push(r);
                    break;

                case ( RuleToken::GREATER_THAN_OPERATOR ):
                    a = results.top(); results.pop();
                    b = results.top(); results.pop();
                    if (a.type == b.type && a.type == RuleToken::NUMERIC_VARIABLE) {
                        r.state = (b.number > a.number);
                    } else {
                        cerr << "cannot compare (>) objects of dissimilar types" << endl;
                        cerr << a.type << " " << b.type << endl;
                        exit(1);
                    }
                    results.push(r);
                    break;

                case ( RuleToken::LESS_THAN_OPERATOR ):
                    a = results.top(); results.pop();
                    b = results.top(); results.pop();
                    if (a.type == b.type && a.type == RuleToken::NUMERIC_VARIABLE) {
                        r.state = (b.number < a.number);
                    } else {
                        cerr << "cannot compare (<) objects of dissimilar types" << endl;
                        cerr << a.type << " " << b.type << endl;
                        exit(1);
                    }
                    results.push(r);
                    break;

                case ( RuleToken::ADD_OPERATOR ):
                    a = results.top(); results.pop();
                    b = results.top(); results.pop();
                    if (a.type == b.type && a.type == RuleToken::NUMERIC_VARIABLE) {
                        r.number = (b.number + a.number);
                        r.type = RuleToken::NUMERIC_VARIABLE;
                    } else {
                        cerr << "cannot add objects of dissimilar types" << endl;
                        cerr << a.type << " " << b.type << endl;
                        exit(1);
                    }
                    results.push(r);
                    break;

                case ( RuleToken::SUBTRACT_OPERATOR ):
                    a = results.top(); results.pop();
                    b = results.top(); results.pop();
                    if (a.type == b.type && a.type == RuleToken::NUMERIC_VARIABLE) {
                        r.number = (b.number - a.number);
                        r.type = RuleToken::NUMERIC_VARIABLE;
                    } else {
                        cerr << "cannot subtract objects of dissimilar types" << endl;
                        cerr << a.type << " " << b.type << endl;
                        exit(1);
                    }
                    results.push(r);
                    break;

                case ( RuleToken::MULTIPLY_OPERATOR ):
                    a = results.top(); results.pop();
                    b = results.top(); results.pop();
                    if (a.type == b.type && a.type == RuleToken::NUMERIC_VARIABLE) {
                        r.number = (b.number * a.number);
                        r.type = RuleToken::NUMERIC_VARIABLE;
                    } else {
                        cerr << "cannot multiply objects of dissimilar types" << endl;
                        cerr << a.type << " " << b.type << endl;
                        exit(1);
                    }
                    results.push(r);
                    break;

                case ( RuleToken::DIVIDE_OPERATOR):
                    a = results.top(); results.pop();
                    b = results.top(); results.pop();
                    if (a.type == b.type && a.type == RuleToken::NUMERIC_VARIABLE) {
                        r.number = (b.number / a.number);
                        r.type = RuleToken::NUMERIC_VARIABLE;
                    } else {
                        cerr << "cannot divide objects of dissimilar types" << endl;
                        cerr << a.type << " " << b.type << endl;
                        exit(1);
                    }
                    results.push(r);
                    break;

                case ( RuleToken::AND_OPERATOR ):
                case ( RuleToken::OR_OPERATOR ):
                    a = results.top(); results.pop();
                    b = results.top(); results.pop();
                    if (a.type == b.type && a.type == RuleToken::BOOLEAN_VARIABLE) {
                        if (token.type == RuleToken::AND_OPERATOR) {
                            r.state = (a.state && b.state);
                        } else {
                            r.state = (a.state || b.state);
                        }
                    } else {
                        cerr << "cannot compare (& or |) objects of dissimilar types" << endl;
                        exit(1);
                    }
                    results.push(r);
                    break;
                default:
                    cerr << "should not get here!" << endl; exit(1);
                    break;
            }
        }
    }
    // at the end you should have only one value on the stack, return it as a boolean
    if (results.size() == 1) {
        if (isBoolean(results.top())) {
            return results.top().state;
        } else {
            cerr << "error, non-boolean value left on stack" << endl;
            //cerr << results.top().value << endl;
            exit(1);
        }
    } else if (results.size() > 1) {
        cerr << "more than one value left on results stack!" << endl;
        while (!results.empty()) {
            cerr << results.top().value << endl;
            results.pop();
        }
        exit(1);
    } else {
        cerr << "results stack empty" << endl;
        exit(1);
    }
}

void VariantFilter::removeFilteredGenotypes(Variant& var, bool keepInfo) {

    for (vector<string>::iterator s = var.sampleNames.begin(); s != var.sampleNames.end(); ++s) {
        string& name = *s;
        if (!passes(var, name)) {
        	if (keepInfo) {
				var.samples[name]["GT"].clear();
				var.samples[name]["GT"].push_back("./.");
        	}
        	else {
			    var.samples.erase(name);
        	}
        }
    }
}

/*
bool VariantCallFile::openVCF(string& filename) {
    file.open(filename.c_str(), ifstream::in);
    if (!file.is_open()) {
        cerr << "could not open " << filename << endl;
        return false;
    } else {
        return parseHeader();
    }
}

bool VariantCallFile::openVCF(ifstream& stream) {
    file = stream;
    if (!file.is_open()) {
        cerr << "provided file is not open" << endl;
        return false;
    } else {
        return parseHeader();
    }
}
*/

void VariantCallFile::updateSamples(vector<string>& newSamples) {
    sampleNames = newSamples;
    // regenerate the last line of the header
    vector<string> headerLines = split(header, '\n');
    vector<string> colnames = split(headerLines.at(headerLines.size() - 1), '\t'); // get the last, update the samples
    vector<string> newcolnames;
    newcolnames.resize(9 + sampleNames.size());
    copy(colnames.begin(), colnames.begin() + 9, newcolnames.begin());
    copy(sampleNames.begin(), sampleNames.end(), newcolnames.begin() + 9);
    headerLines.at(headerLines.size() - 1) = join(newcolnames, "\t");
    header = join(headerLines, "\n");
}

// non-destructive version of above
string VariantCallFile::headerWithSampleNames(vector<string>& newSamples) {
    // regenerate the last line of the header
    if (newSamples.empty()) return header;
    vector<string> headerLines = split(header, '\n');
    vector<string> colnames = split(headerLines.at(headerLines.size() - 1), '\t'); // get the last, update the samples
    vector<string> newcolnames;
    unsigned int colCount = colnames.size(); // used to be hard-coded 9, hopefully the dynamic colCount isn't an issue
    if (colCount < 8)
    {
        cout << "VCF file is not suitable for use because it does not have a format field." << endl;
        exit(0);
    }
    newcolnames.resize(colCount + newSamples.size());
    copy(colnames.begin(), colnames.begin() + colCount, newcolnames.begin());
    copy(newSamples.begin(), newSamples.end(), newcolnames.begin() + colCount);
    headerLines.at(headerLines.size() - 1) = join(newcolnames, "\t");
    return join(headerLines, "\n");
}

// TODO cleanup, store header lines instead of bulk header
void VariantCallFile::addHeaderLine(string line) {
    vector<string> headerLines = split(header, '\n');
    headerLines.insert(headerLines.end() - 1, line);
    header = join(unique(headerLines), "\n");
}

// helper to addHeaderLine
vector<string>& unique(vector<string>& strings) {
    set<string> uniq;
    vector<string> res;
    for (const auto& s : strings) {
        if (uniq.find(s) == uniq.end()) {
            res.push_back(s);
            uniq.insert(s);
        }
    }
    strings = res;
    return strings;
}

vector<string> VariantCallFile::infoIds(void) {
    vector<string> tags;
    vector<string> headerLines = split(header, '\n');
    for (const auto& line : headerLines) {
        if (line.find("##INFO") == 0) {
            size_t pos = line.find("ID=");
            if (pos != string::npos) {
                pos += 3;
                size_t tagend = line.find(",", pos);
                if (tagend != string::npos) {
                    tags.push_back(line.substr(pos, tagend - pos));
                }
            }
        }
    }
    return tags;
}

vector<string> VariantCallFile::formatIds(void) {
    vector<string> tags;
    vector<string> headerLines = split(header, '\n');
    for (const auto& line : headerLines) {
        if (line.find("##FORMAT") == 0) {
            size_t pos = line.find("ID=");
            if (pos != string::npos) {
                pos += 3;
                size_t tagend = line.find(",", pos);
                if (tagend != string::npos) {
                    tags.push_back(line.substr(pos, tagend - pos));
                }
            }
        }
    }
    return tags;
}

void VariantCallFile::removeInfoHeaderLine(string const & tag) {
    vector<string> headerLines = split(header, '\n');
    vector<string> newHeader;
    string id = "ID=" + tag + ",";
    for (const auto& line : headerLines) {
        if (line.find("##INFO") == 0) {
            if (line.find(id) == string::npos) {
                newHeader.push_back(line);
            }
        } else {
            newHeader.push_back(line);
        }
    }
    header = join(newHeader, "\n");
}

void VariantCallFile::removeGenoHeaderLine(string const & tag) {
    vector<string> headerLines = split(header, '\n');
    vector<string> newHeader;
    string id = "ID=" + tag + ",";
    for (const auto& headerLine : headerLines) {
        if (headerLine.find("##FORMAT") == 0) {
            if (headerLine.find(id) == string::npos) {
                newHeader.push_back(headerLine);
            }
        } else {
            newHeader.push_back(headerLine);
        }
    }
    header = join(newHeader, "\n");
}

vector<string> VariantCallFile::getHeaderLinesFromFile()
{
    string headerStr = "";

    if (usingTabix) {
        tabixFile->getHeader(headerStr);
        if (headerStr.empty()) {
            cerr << "error: no VCF header" << endl;
            exit(1);
        }
        tabixFile->getNextLine(line);
        firstRecord = true;
    } else {
        while (std::getline(*file, line)) {
            if (line.substr(0,1) == "#") {
                headerStr += line + '\n';
            } else {
                // done with header
                if (headerStr.empty()) {
                    cerr << "error: no VCF header" << endl;
                    return vector<string>();
                }
                firstRecord = true;
                break;
            }
        }
    }
    return split(headerStr, "\n");
}

bool VariantCallFile::parseHeader(void) {

    string headerStr = "";

    if (usingTabix) {
        tabixFile->getHeader(headerStr);
        if (headerStr.empty()) {
            cerr << "error: no VCF header" << endl;
            exit(1);
        }
        tabixFile->getNextLine(line);
        firstRecord = true;
    } else {
        while (std::getline(*file, line)) {
            if (line.substr(0,1) == "#") {
                headerStr += line + '\n';
            } else {
                // done with header
                if (headerStr.empty()) {
                    cerr << "error: no VCF header" << endl;
                    return false;
                }
                firstRecord = true;
                break;
            }
        }
    }
    this->vcf_header = headerStr;

    return parseHeader(headerStr);

}

bool VariantCallFile::parseHeader(string& hs) {

    if (hs.empty()) return false;
    if (hs.substr(hs.size() - 1, 1) == "\n") {
	hs.erase(hs.size() - 1, 1); // remove trailing newline
    }
    header = hs; // stores the header in the object instance

    vector<string> headerLines = split(header, "\n");
    for (const auto& headerLine : headerLines) {
        if (headerLine.substr(0,2) == "##") {
            // meta-information headerLines
            // TODO parse into map from info/format key to type
            // ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
            // ##FORMAT=<ID=CB,Number=1,Type=String,Description="Called by S(Sanger), M(UMich), B(BI)">
            size_t found = headerLine.find_first_of("=");
            string entryType = headerLine.substr(2, found - 2);
            // handle reference here, no "<" and ">" given
                //} else if (entryType == "reference") {
            size_t dataStart = headerLine.find_first_of("<");
            size_t dataEnd = headerLine.find_first_of(">");
            if (dataStart != string::npos && dataEnd != string::npos) {
                string entryData = headerLine.substr(dataStart + 1, dataEnd - dataStart - 1);
                // this will break if it includes a "long form" string
                // including either a = or , in the first or second
                // field
                if (entryType == "INFO" || entryType == "FORMAT") {
                    vector<string> fields = split(entryData, "=,");
                    if (fields.size() < 8) {
                        cerr << "header line does not have all of the required fields: ID, Number, Type, and Description" << endl
                             << headerLine << endl;
                        exit(1);
                    }
                    // get the required fields from the header line
                    auto id_field = find(fields.begin(), fields.begin() + 8, "ID");
                    auto num_field = find(fields.begin(), fields.begin() + 8, "Number");
                    auto type_field = find(fields.begin(), fields.begin() + 8, "Type");
                    auto desc_field = find(fields.begin(), fields.begin() + 8, "Description");
                    for (auto it : {id_field, num_field, type_field, desc_field}) {
                        // make sure we found the field and that all of the keys have a value associated
                        if (it == fields.begin() + 8 || ((it - fields.begin()) % 2 == 1)) {
                            if (it == desc_field) {
                                // we don't actually record / use the description, so we'll just give a warning
                                cerr << "warning: ";
                            }
                            cerr << "header line does not have all of the required fields (ID, Number, Type, and Description) in the first 4 fields" << endl
                                 << headerLine << endl;
                            if (it != desc_field) {
                                exit(1);
                            }
                        }
                    }
                    string id = *(id_field + 1);
                    int number;
                    string numberstr = *(num_field + 1);
                    // string numberstr = mapper["Number"].c_str();

                    // XXX TODO VCF has variable numbers of fields...
                    if (numberstr == "A") {
                        number = ALLELE_NUMBER;
                    } else if (numberstr == "G") {
                        number = GENOTYPE_NUMBER;
                    } else if (numberstr == ".") {
                        number = 1;
                    } else {
                        convert(numberstr, number);
                    }
                    VariantFieldType type = typeStrToVariantFieldType(*(type_field + 1));
                    // VariantFieldType type = typeStrToVariantFieldType(mapper["TYPE"]);
                    if (entryType == "INFO") {
                        infoCounts[id] = number;
                        infoTypes[id] = type;
                    } else if (entryType == "FORMAT") {
                        formatCounts[id] = number;
                        formatTypes[id] = type;
                    }
                }
            }
        } else if (headerLine.substr(0,1) == "#") {
            // field name headerLine
            vector<string> fields = split(headerLine, '\t');
            if (fields.size() > 8) {
                sampleNames.resize(fields.size() - 9);
                copy(fields.begin() + 9, fields.end(), sampleNames.begin());
            }
        }
    }

    return true;
}

bool VariantCallFile::getNextVariant(Variant& var) {
        if (firstRecord && !justSetRegion) {
            if (!line.empty() && line.substr(0,1) != "#") {
                var.parse(line, parseSamples);
                firstRecord = false;
                _done = false;
                return true;
            } else {
                return false;
            }
        }
        if (usingTabix) {
            if (justSetRegion && !line.empty() && line.substr(0,1) != "#") {
                if (firstRecord) {
                    firstRecord = false;
                }
                var.parse(line, parseSamples);
                line.clear();
                justSetRegion = false;
                _done = false;
                return true;
            } else if (tabixFile->getNextLine(line)) {
                var.parse(line, parseSamples);
                _done = false;
                return true;
            } else {
                _done = true;
                return false;
            }
        } else {
            if (std::getline(*file, line)) {
                var.parse(line, parseSamples);
                _done = false;
                return true;
            } else {
                _done = true;
                return false;
            }
        }
}

bool VariantCallFile::setRegion(const string& seq, long int start, long int end) {
    stringstream regionstr;
    if (end) {
        regionstr << seq << ":" << start << "-" << end;
    } else {
        regionstr << seq << ":" << start;
    }
    return setRegion(regionstr.str());
}

bool VariantCallFile::setRegion(const string& region) {
    if (!usingTabix) {
        cerr << "cannot setRegion on a non-tabix indexed file" << endl;
        exit(1);
    }
    // convert between bamtools/freebayes style region string and tabix/samtools style
    regex txt_regex("(\\d+)\\.\\.(\\d+)$");
    string tabix_region = regex_replace(region, txt_regex, "$1-$2");

    if (tabixFile->setRegion(tabix_region)) {
        if (tabixFile->getNextLine(line)) {
	    justSetRegion = true;
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}


// genotype manipulation
/*
map<string, int> decomposeGenotype(string& genotype) {
    string splitter = "/";
    if (genotype.find("|") != string::npos) {
        splitter = "|";
    }
    vector<string> haps = split(genotype, splitter);
    map<string, int> decomposed;
    for (vector<string>::iterator h = haps.begin(); h != haps.end(); ++h) {
        ++decomposed[*h];
    }
    return decomposed;
}
*/

map<int, int> decomposeGenotype(const string& genotype) {
    string splitter = "/";
    if (genotype.find('|') != string::npos) {
        splitter = "|";
    }
    vector<string> haps = split(genotype, splitter);
    map<int, int> decomposed;
    for (const auto& h : haps) {
        int alt;
        if (h == ".") {
            ++decomposed[NULL_ALLELE];
        } else {
            convert(h, alt);
            ++decomposed[alt];
        }
    }
    return decomposed;
}

vector<int> decomposePhasedGenotype(const string& genotype) {
    string splitter = "/";
    if (genotype.find('|') != string::npos) {
        splitter = "|";
    }
    vector<string> haps = split(genotype, splitter);
    if (haps.size() > 1 && splitter == "/") {
        cerr << "could not find '|' in genotype, cannot decomposePhasedGenotype on unphased genotypes" << endl;
        exit(1);
    }
    vector<int> decomposed;
    for (const auto& h : haps) {
        int alt;
        if (h == ".") {
            decomposed.push_back(NULL_ALLELE);
        } else {
            convert(h, alt);
            decomposed.push_back(alt);
        }
    }
    return decomposed;
}

string genotypeToString(const map<int, int>& genotype) {
    vector<int> s;
    for (const auto& g : genotype) {
        int a = g.first;
        int c = g.second;
        for (int i = 0; i < c; ++i) s.push_back(a);
    }
    sort(s.begin(), s.end());
    vector<string> r;
    r.reserve(s.size());

    for (const auto& i : s) {
        if (i == NULL_ALLELE) r.push_back(".");
        else r.push_back(convert(i));
    }
    return join(r, "/"); // TODO adjust for phased/unphased
}

string phasedGenotypeToString(const vector<int>& genotype) {
    vector<string> r;
    for (const auto& i : genotype) {
        if (i == NULL_ALLELE) r.push_back(".");
        else r.push_back(convert(i));
    }
    return join(r, "|");
}

bool isHet(const map<int, int>& genotype) {
    return genotype.size() > 1;
}

bool isHom(const map<int, int>& genotype) {
    return genotype.size() == 1;
}

bool hasNonRef(const map<int, int>& genotype) {
    for (const auto& g : genotype) {
        if (g.first != 0) {
            return true;
        }
    }
    return false;
}

bool isHomRef(const map<int, int>& genotype) {
    return isHom(genotype) && !hasNonRef(genotype);
}

bool isHomNonRef(const map<int, int>& genotype) {
    return isHom(genotype) && hasNonRef(genotype);
}

bool isNull(const map<int, int>& genotype) {
    return genotype.find(NULL_ALLELE) != genotype.end();
}

int ploidy(const map<int, int>& genotype) {
    int i = 0;
    for (const auto& g : genotype) {
        i += g.second;
    }
    return i;
}




map<string, vector<VariantAllele> > Variant::flatAlternates(void) {
    map<string, vector<VariantAllele> > variantAlleles;
    for (const auto& alternate: alt) {
        vector<VariantAllele>& variants = variantAlleles[alternate];
        variants.push_back(VariantAllele(ref, alternate, position));
    }
    return variantAlleles;
}

set<string> Variant::altSet(void) {
    set<string> altset(alt.begin(), alt.end());
    return altset;
}

map<pair<int, int>, int> Variant::getGenotypeIndexesDiploid(void) {

    map<pair<int, int>, int> genotypeIndexes;
    //map<int, map<Genotype*, int> > vcfGenotypeOrder;
    vector<int> indexes;
    indexes.reserve(alleles.size());
    for (int i = 0; i < alleles.size(); ++i) {
        indexes.push_back(i);
    }
    int ploidy = 2; // ONLY diploid
    vector<vector<int> > genotypes = multichoose(ploidy, indexes);
    for (auto& g : genotypes) {
        sort(g.begin(), g.end());  // enforce e.g. 0/1, 0/2, 1/2 ordering over reverse
        // XXX this does not handle non-diploid!!!!
        int j = g.front();
        int k = g.back();
        genotypeIndexes[make_pair(j, k)] = (k * (k + 1) / 2) + j;
    }
    return genotypeIndexes;
}

void Variant::updateAlleleIndexes(void) {
    // adjust the allele index
    altAlleleIndexes.clear();
    int m = 0;
    for (vector<string>::iterator a = alt.begin();
            a != alt.end(); ++a, ++m) {
        altAlleleIndexes[*a] = m;
    }
}

// TODO only works on "A"llele variant fields
  void Variant::removeAlt(const string& altAllele) {

    int altIndex = getAltAlleleIndex(altAllele);  // this is the alt-relative index, 0-based

    for (const auto& c: vcf->infoCounts) {
      int count = c.second;
      if (count == ALLELE_NUMBER) {
	const string& key = c.first;
	map<string, vector<string> >::iterator v = info.find(key);
	if (v != info.end()) {
	  vector<string>& vals = v->second;
	  vector<string> tokeep;
	  int i = 0;
	  for (vector<string>::iterator a = vals.begin();
	       a != vals.end(); ++a, ++i) {
	    if (i != altIndex) {
	      tokeep.push_back(*a);
	    }
	  }
	  vals = tokeep;
	}
      }
    }

    for (const auto& c : vcf->formatCounts) {
      int count = c.second;
      if (count == ALLELE_NUMBER) {
      	const string& key = c.first;
      	for (auto& [_, sample] : samples) {
      		map<string, vector<string> >::iterator v = sample.find(key);
      		if (v != sample.end()) {
      			vector<string>& vals = v->second;
      			vector<string> tokeep;
      			int i = 0;
      			for (vector<string>::iterator a = vals.begin();
                    a != vals.end(); ++a, ++i) {
      				if (i != altIndex) {
      					tokeep.push_back(*a);
      				}
                    }
      			vals = tokeep;
      		}
      	}
      }
    }

    int altSpecIndex = altIndex + 1; // this is the genotype-spec index, ref=0, 1-based for alts

    vector<string> newalt;
    map<int, int> alleleIndexMapping;
    // setup the new alt string
    alleleIndexMapping[0] = 0; // reference allele remains the same
    alleleIndexMapping[NULL_ALLELE] = NULL_ALLELE; // null allele remains the same
    int i = 1; // current index
    int j = 1; // new index
    for (vector<string>::iterator a = alt.begin(); a != alt.end(); ++a, ++i) {
        if (i != altSpecIndex) {
            newalt.push_back(*a);
            // get the mapping between new and old allele indexes
            alleleIndexMapping[i] = j;
            ++j;
        } else {
            alleleIndexMapping[i] = NULL_ALLELE;
        }
    }

    // fix the sample genotypes, removing reference to the old allele
    map<string, int> samplePloidy;
    for (map<string, map<string, vector<string> > >::iterator s = samples.begin(); s != samples.end(); ++s) {
        map<string, vector<string> >& sample = s->second;
        if (sample.find("GT") != sample.end()) {
            string& gt = sample["GT"].front();
            string splitter = "/";
            if (gt.find('|') != string::npos) {
                splitter = "|";
            }

            if (splitter == "/") {
                samplePloidy[s->first] = split(gt, splitter).size();
                map<int, int> genotype = decomposeGenotype(sample["GT"].front());
                map<int, int> newGenotype;
                for (map<int, int>::iterator g = genotype.begin(); g != genotype.end(); ++g) {
                    newGenotype[alleleIndexMapping[g->first]] += g->second;
                }
                sample["GT"].clear();
                sample["GT"].push_back(genotypeToString(newGenotype));
            } else {
                samplePloidy[s->first] = split(gt, splitter).size();
                vector<int> genotype = decomposePhasedGenotype(sample["GT"].front());
                vector<int> newGenotype;
                for (vector<int>::iterator g = genotype.begin(); g != genotype.end(); ++g) {
                    newGenotype.push_back(alleleIndexMapping[*g]);
                }
                sample["GT"].clear();
                sample["GT"].push_back(phasedGenotypeToString(newGenotype));
            }
        }
    }

    set<int> ploidies;
    for (map<string, int>::iterator p = samplePloidy.begin(); p != samplePloidy.end(); ++p) {
        ploidies.insert(p->second);
    }

    // fix the sample genotype likelihoods, removing reference to the old allele
    // which GL fields should we remove?
    vector<int> toRemove;
    toRemove.push_back(altSpecIndex);
    map<int, map<int, int> > glMappingByPloidy;
    for (const auto p : ploidies) {
        glMappingByPloidy[p] = glReorder(p, alt.size() + 1, alleleIndexMapping, toRemove);
    }

    for (auto& s : samples) {
        auto& sample = s.second;
        auto glsit = sample.find("GL");
        if (glsit != sample.end()) {
            vector<string>& gls = glsit->second; // should be split already
            map<int, string> newgls;
            map<int, int>& newOrder = glMappingByPloidy[samplePloidy[s.first]];
            int i = 0;
            for (vector<string>::iterator g = gls.begin(); g != gls.end(); ++g, ++i) {
                int j = newOrder[i];
                if (j != -1) {
                    newgls[i] = *g;
                }
            }
            // update the gls
            gls.clear();
            for (const auto& g : newgls) {
                gls.push_back(g.second);
            }
        }
    }

    // reset the alt
    alt = newalt;

    // and the alleles
    alleles.clear();
    alleles.push_back(ref);
    alleles.insert(alleles.end(), alt.begin(), alt.end());

    updateAlleleIndexes();

}

// union of lines in headers of input files
string unionInfoHeaderLines(string& s1, string& s2) {
    vector<string> lines1 = split(s1, "\n");
    vector<string> lines2 = split(s2, "\n");
    vector<string> result;
    set<string> l2;
    string lastHeaderLine; // this one needs to be at the end
    for (const auto& s : lines2) {
        if (s.substr(0,6) == "##INFO") {
            l2.insert(s);
        }
    }
    for (const auto& s : lines1) {
        if (l2.count(s)) {
            l2.erase(s);
        }
        if (s.substr(0,6) == "#CHROM") {
            lastHeaderLine = s;
        } else {
            result.push_back(s);
        }
    }
    for (const auto& s : l2) {
        result.push_back(s);
    }
    if (lastHeaderLine.empty()) {
        cerr << "could not find CHROM POS ... header line" << endl;
        exit(1);
    }
    result.push_back(lastHeaderLine);
    return join(result, "\n");
}

list<list<int> > _glorder(int ploidy, int alts) {
    if (ploidy == 1) {
        list<list<int> > results;
        for (int n = 0; n < alts; ++n) {
            list<int> v;
            v.push_back(n);
            results.push_back(v);
        }
        return results;
    } else {
        list<list<int> > results;
        for (int n = 0; n < alts; ++n) {
            list<list<int> > x = _glorder(ploidy - 1, alts);
            for (auto& v : x) {
                if (v.front() <= n) {
                    v.push_front(n);
                    results.push_back(v);
                }
            }
        }
        return results;
    }
}

// genotype likelihood-ordering of genotypes, where each genotype is a
// list of integers (as written in the GT field)
list<list<int> > glorder(int ploidy, int alts) {
    list<list<int> > results = _glorder(ploidy, alts);
    for (auto& v : results) {
        v.reverse();
    }
    return results;
}

// which genotype likelihoods would include this alternate allele
list<int> glsWithAlt(int alt, int ploidy, int numalts) {
    list<int> gls;
    list<list<int> > orderedGenotypes = glorder(ploidy, numalts);
    int i = 0;
    for (list<list<int> >::iterator v = orderedGenotypes.begin(); v != orderedGenotypes.end(); ++v, ++i) {
        for (const auto& q : *v) {
            if (q == alt) {
                gls.push_back(i);
                break;
            }
        }
    }
    return gls;
}

// describes the mapping between the old gl ordering and and a new
// one in which the GLs including the old alt have been removed
// a map to -1 means "remove"
map<int, int> glReorder(int ploidy, int numalts, map<int, int>& alleleIndexMapping, vector<int>& altsToRemove) {
    map<int, int> mapping;
    list<list<int> > orderedGenotypes = glorder(ploidy, numalts);
    for (auto& v : orderedGenotypes) {
        for (auto& n : v) {
            n = alleleIndexMapping[n];
        }
    }
    list<list<int> > newOrderedGenotypes = glorder(ploidy, numalts - altsToRemove.size());
    map<list<int>, int> newOrderedGenotypesMapping;
    int i = 0;
    // mapping is wrong...
    for (list<list<int> >::iterator v = newOrderedGenotypes.begin(); v != newOrderedGenotypes.end(); ++v, ++i) {
        newOrderedGenotypesMapping[*v] = i;
    }
    i = 0;
    for (list<list<int> >::iterator v = orderedGenotypes.begin(); v != orderedGenotypes.end(); ++v, ++i) {
        map<list<int>, int>::iterator m = newOrderedGenotypesMapping.find(*v);
        if (m != newOrderedGenotypesMapping.end()) {
            //cout << "new gl order of " << i << " is " << m->second << endl;
            mapping[i] = m->second;
        } else {
            //cout << i << " will be removed" << endl;
            mapping[i] = -1;
        }
    }
    return mapping;
}

string Variant::getGenotype(const string& sample) {
    map<string, map<string, vector<string> > >::iterator s = samples.find(sample);
    if (s != samples.end()) {
        map<string, vector<string> >::iterator f = s->second.find("GT");
        if (f != s->second.end()) {
            return f->second.front();
        }
    }
    return "";
}

bool Variant::isPhased(void) {
    for (map<string, map<string, vector<string> > >::iterator s = samples.begin(); s != samples.end(); ++s) {
        map<string, vector<string> >& sample = s->second;
        map<string, vector<string> >::iterator g = sample.find("GT");
        if (g != sample.end()) {
            string gt = g->second.front();
            if (gt.size() > 1 && gt.find('|') == string::npos) {
                return false;
            }
        }
    }
    return true;
}

long Variant::zeroBasedPosition(void) const {
    return position - 1;
}

string Variant::vrepr(void) {
    return sequenceName + "\t" + convert(position) + "\t" + join(alleles, ",");
}

// TODO
/*
vector<Variant*> Variant::matchingHaplotypes() {

    int haplotypeStart = var.position;
    int haplotypeEnd = var.position + var.ref.size();

    for (vector<Variant*>::iterator v = overlapping.begin(); v != overlapping.end(); ++v) {
        haplotypeStart = min((*v)->position, (long int) haplotypeStart);
        haplotypeEnd = max((*v)->position + (*v)->ref.size(), (long unsigned int) haplotypeEnd);
    }

    // for everything overlapping and the current variant, construct the local haplotype within the bounds
    // if there is an exact match, the allele in the current VCF does intersect

    string referenceHaplotype = reference.getSubSequence(var.sequenceName, haplotypeStart - 1, haplotypeEnd - haplotypeStart);
    map<string, vector<pair<Variant*, int> > > haplotypes; // map to variant and alt index

    for (vector<Variant*>::iterator v = overlapping.begin(); v != overlapping.end(); ++v) {
        Variant& variant = **v;
        int altindex = 0;
        for (vector<string>::iterator a = variant.alt.begin(); a != variant.alt.end(); ++a, ++altindex) {
            string haplotype = referenceHaplotype;
            // get the relative start and end coordinates for the variant alternate allele
            int relativeStart = variant.position - haplotypeStart;
            haplotype.replace(relativeStart, variant.ref.size(), *a);
            haplotypes[haplotype].push_back(make_pair(*v, altindex));
        }
    }

    Variant originalVar = var;

    // determine the non-intersecting alts
    vector<string> altsToRemove;
    vector<int> altIndexesToRemove;
    for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
        string haplotype = referenceHaplotype;
        int relativeStart = var.position - haplotypeStart;
        haplotype.replace(relativeStart, var.ref.size(), *a);
        map<string, vector<pair<Variant*, int> > >::iterator h = haplotypes.find(haplotype);
        if ((intersecting && !invert && h == haplotypes.end())
            || (intersecting && invert && h != haplotypes.end())
            || (unioning && h != haplotypes.end())) {
            if (tag.empty() && mergeToTag.empty()) {
                altsToRemove.push_back(*a);
            } else {
                if (!tag.empty()) {
                    var.info[tag].push_back(".");
                }
                if (!mergeToTag.empty()) {
                    var.info[mergeToTag].push_back(".");
                }
            }
        } else {
            if (!tag.empty()) {
                var.info[tag].push_back(tagValue);
            }
            // NB: just take the first value for the mergeFromTag
            if (!mergeToTag.empty()) {
                Variant* v = h->second.front().first;
                int index = h->second.front().second;
                if (v->info.find(mergeFromTag) != v->info.end()) {
                    // now you have to find the exact allele...
                    string& otherValue = v->info[mergeFromTag].at(index);
                    var.info[mergeToTag].push_back(otherValue);
                } else if (mergeFromTag == "QUAL") {
                    var.info[mergeToTag].push_back(convert(v->quality));
                } else {
                    var.info[mergeToTag].push_back(".");
                }
            }
        }
    }

    // remove the non-overlapping (intersecting) or overlapping (unioning) alts
    if (intersecting && loci && altsToRemove.size() != var.alt.size()) {
        // we have a match in loci mode, so we should output the whole loci, not just the matching sequence
    } else {
        for (vector<string>::iterator a = altsToRemove.begin(); a != altsToRemove.end(); ++a) {
            var.removeAlt(*a);
        }
    }

    if (unioning) {

        // somehow sort the records and combine them?
        map<long int, vector<Variant*> > variants;
        for (vector<Variant*>::iterator o = overlapping.begin(); o != overlapping.end(); ++o) {
            if ((*o)->position <= var.position && // check ensures proper ordering of variants on output
                outputVariants.find(*o) == outputVariants.end()) {
                outputVariants.insert(*o);
                variants[(*o)->position].push_back(*o);
            }
        }
        // add in the current variant, if it has alts left
        if (!var.alt.empty()) {
            vector<Variant*>& vars = variants[var.position];
            int numalts = 0;
            for (vector<Variant*>::iterator v = vars.begin(); v != vars.end(); ++v) {
                numalts += (*v)->alt.size();
            }
            if (numalts + var.alt.size() == originalVar.alt.size()) {
                variants[var.position].clear();
                variants[var.position].push_back(&originalVar);
            } else {
                variants[var.position].push_back(&var);
            }
        }

        for (map<long int, vector<Variant*> >::iterator v = variants.begin(); v != variants.end(); ++v) {
            for (vector<Variant*>::iterator o = v->second.begin(); o != v->second.end(); ++o) {
                cout << **o << endl;
                lastOutputPosition = max(lastOutputPosition, (*o)->position);
            }
        }
    } else {
        // if any alts remain, output the variant record
        if (!var.alt.empty()) {
            cout << var << endl;
            lastOutputPosition = max(lastOutputPosition, var.position);
        }
    }

}
*/


    VCFHeader::VCFHeader()
    {

        // add manditory fields
        this->header_columns.push_back("#CHROM");
        this->header_columns.push_back("POS");
        this->header_columns.push_back("ID");
        this->header_columns.push_back("REF");
        this->header_columns.push_back("ALT");
        this->header_columns.push_back("QUAL");
        this->header_columns.push_back("FILTER");
        this->header_columns.push_back("INFO");

        // add the line names in order
        // the order is used when outputting as a string
        this->header_line_names_ordered.push_back("##fileFormat");
        this->header_line_names_ordered.push_back("##fileDate");
        this->header_line_names_ordered.push_back("##source");
        this->header_line_names_ordered.push_back("##reference");
        this->header_line_names_ordered.push_back( "##contig");
        this->header_line_names_ordered.push_back("##phasing");
        this->header_line_names_ordered.push_back( "##assembly");

        // add the list names in order
        // the order is used when outputting as a string (getHeaderString)
        this->header_list_names_ordered.push_back("##info");
        this->header_list_names_ordered.push_back("##filter");
        this->header_list_names_ordered.push_back("##format");
        this->header_list_names_ordered.push_back("##alt");
        this->header_list_names_ordered.push_back("##sample");
        this->header_list_names_ordered.push_back("##pedigree");
        this->header_list_names_ordered.push_back("##pedigreedb");

        // initialize the header_lines with the above vector.
        // Set the key as the ##_type_ and the value as an empty string
        // Empty strings are ignored when outputting as string (getHeaderString)
        for (const auto& header_lines_iter : this->header_line_names_ordered)
        {
            this->header_lines[header_lines_iter] = "";
        }

        // initialize the header_lines with the above vector.
        // Set the key as the ##_type_ and the value as an empty vector<string>
        // Empty vectors are ignored when outputting as string (getHeaderString)
        for (const auto& header_lists_iter : header_list_names_ordered)
        {
            this->header_lists[header_lists_iter] = vector<string>(0);
        }

    }

    void VCFHeader::addMetaInformationLine(const string& meta_line)
    {
        // get the meta_line unique key (first chars before the =)
        unsigned int meta_line_index = meta_line.find('=', 0);
        string meta_line_prefix = meta_line.substr(0, meta_line_index);

        // check if the meta_line_prefix is in the header_lines, if so add it to the appropirate list
        if (this->header_lines.find(meta_line_prefix) != header_lines.end()) // the meta_line is a header line so replace what was there
        {
            this->header_lines[meta_line_prefix] = meta_line;
        }
        else if (header_lists.find(meta_line_prefix) != header_lists.end() &&
            !metaInfoIdExistsInVector(meta_line, this->header_lists[meta_line_prefix])) // check if the metalineprefix is in the headerLists, if so add it to the appropirate list
        {
            this->header_lists[meta_line_prefix].push_back(meta_line);
        }
    }

    string VCFHeader::getHeaderString()
    {
        // getHeaderString generates the string each time it is called
        string header_string;

        // start by adding the header_lines
        for (const auto& header_lines_iter : header_line_names_ordered)
        {
            if (this->header_lines[header_lines_iter] != "")
            {
                header_string += this->header_lines[header_lines_iter] + "\n";
            }
        }

        // next add header_lists
        for (const auto& header_lists_iter : header_list_names_ordered)
        {
            vector<string> tmp_header_lists = this->header_lists[header_lists_iter];
            for (const auto& header_list : tmp_header_lists)
            {
                header_string += header_list + "\n";
            }
        }

        // last add header columns
        const auto last_element = this->header_columns.end() - 1;
        for (auto header_column_iter = this->header_columns.begin(); header_column_iter != this->header_columns.end(); ++header_column_iter)
        {
            string delimiter = (header_column_iter == last_element) ? "\n" : "\t";
            header_string += (*header_column_iter) + delimiter;
        }
        return header_string;
    }

    bool VCFHeader::metaInfoIdExistsInVector(const string& meta_line, vector<string>& meta_lines)
    {
        // extract the id from meta_line
        size_t meta_line_id_start_idx = meta_line.find("ID=", 0); // used for the start of the substring index
        size_t meta_line_id_end_idx = meta_line.find(',', meta_line_id_start_idx); // used for end of the substring index
        string meta_line_id = (meta_line_id_start_idx < meta_line_id_end_idx) ? meta_line.substr(meta_line_id_start_idx, meta_line_id_end_idx - meta_line_id_start_idx) : "";

        for (const auto& meta_line : meta_lines)
        {
            // extract the id from meta_line string
            size_t meta_line_id_start_idx = meta_line.find("ID=", 0);
            size_t meta_line_id_end_idx = meta_line.find(",", meta_line_id_start_idx);
            string meta_line_id = (meta_line_id_start_idx < meta_line_id_end_idx) ? meta_line.substr(meta_line_id_start_idx, meta_line_id_end_idx - meta_line_id_start_idx) : "";
            // compare the meta_line_id with the meta_line_id
            if (strcasecmp(meta_line_id.c_str(), meta_line_id.c_str()) == 0)
            {
                return true;
            }
        }
        return false;
    }

    void VCFHeader::addHeaderColumn(const string& header_column)
    {
        // don't add duplicates
        //  vector<string>::iterator test = find(this->header_columns.begin(), this->header_columns.end(), header_column);
        if (find(this->header_columns.begin(), this->header_columns.end(), header_column) == this->header_columns.end())
        {
            this->header_columns.push_back(header_column);
        }
    }


// Reintroduce old version for freebayes

map<string, vector<VariantAllele> > Variant::parsedAlternates(bool includePreviousBaseForIndels,
                                                              bool useMNPs,
                                                              bool useEntropy,
                                                              float matchScore,
                                                              float mismatchScore,
                                                              float gapOpenPenalty,
                                                              float gapExtendPenalty,
                                                              float repeatGapExtendPenalty,
                                                              const string& flankingRefLeft,
                                                              const string& flankingRefRight) {

    map<string, vector<VariantAllele> > variantAlleles;

    if (isSymbolicSV()){
        // Don't ever align SVs. It just wrecks things.
        return this->flatAlternates();
    }
    // add the reference allele
    variantAlleles[ref].push_back(VariantAllele(ref, ref, position));

    // single SNP case, no ambiguity possible, no need to spend a lot of
    // compute aligning ref and alt fields
    if (alt.size() == 1 && ref.size() == 1 && alt.front().size() == 1) {
        variantAlleles[alt.front()].push_back(VariantAllele(ref, alt.front(), position));
        return variantAlleles;
    }

    // padding is used to ensure a stable alignment of the alternates to the reference
    // without having to go back and look at the full reference sequence
    int paddingLen = max(10, (int) (ref.size()));  // dynamically determine optimum padding length
    for (const auto& alternate : alt) {
        paddingLen = max(paddingLen, (int) (alternate.size()));
    }
    char padChar = 'Z';
    char anchorChar = 'Q';
    string padding(paddingLen, padChar);

    // this 'anchored' string is done for stability
    // the assumption is that there should be a positional match in the first base
    // this is true for VCF 4.1, and standard best practices
    // using the anchor char ensures this without other kinds of realignment
    string reference_M;
    if (flankingRefLeft.empty() && flankingRefRight.empty()) {
        reference_M = padding + ref + padding;
        reference_M[paddingLen] = anchorChar;
    } else {
        reference_M = flankingRefLeft + ref + flankingRefRight;
        paddingLen = flankingRefLeft.size();
    }

    // passed to sw.Align
    unsigned int referencePos;

    string cigar;

    for (const auto& alternate : alt) {

      vector<VariantAllele>& variants = variantAlleles[alternate];
      string alternateQuery_M;
      if (flankingRefLeft.empty() && flankingRefRight.empty()) {
	alternateQuery_M = padding + alternate + padding;
	alternateQuery_M[paddingLen] = anchorChar;
      } else {
	alternateQuery_M = flankingRefLeft + alternate + flankingRefRight;
      }
      //const unsigned int alternateLen = alternate.size();

      if (true) {
	CSmithWatermanGotoh sw(matchScore,
			       mismatchScore,
			       gapOpenPenalty,
			       gapExtendPenalty);
	if (useEntropy) sw.EnableEntropyGapPenalty(1);
	if (repeatGapExtendPenalty != 0){
	  sw.EnableRepeatGapExtensionPenalty(repeatGapExtendPenalty);
	}
	sw.Align(referencePos, cigar, reference_M, alternateQuery_M);
      } else {  // disabled for now
	StripedSmithWaterman::Aligner aligner;
	StripedSmithWaterman::Filter sswFilter;
	StripedSmithWaterman::Alignment alignment;
	aligner.Align(alternateQuery_M.c_str(),
		      reference_M.c_str(),
		      reference_M.size(), sswFilter, &alignment);
	cigar = alignment.cigar_string;
      }

      // left-realign the alignment...

      vector<pair<int, string> > cigarData = old_splitCigar(cigar);

      if (cigarData.front().second != "M"
	  || cigarData.back().second != "M"
	  || cigarData.front().first < paddingLen
	  || cigarData.back().first < paddingLen) {
	cerr << "parsedAlternates: alignment does not start with match over padded sequence" << endl;
	cerr << cigar << endl;
	cerr << reference_M << endl;
	cerr << alternateQuery_M << endl;
	exit(1);
      } else {
	cigarData.front().first -= paddingLen;
	cigarData.back().first -= paddingLen;;
      }
      //cigarData = cleanCigar(cigarData);
      cigar = old_joinCigar(cigarData);

      int altpos = 0;
      int refpos = 0;

      for (const auto& e : cigarData) {

	int len = e.first;
	const string& type = e.second;

	switch (type.at(0)) {
	case 'I':
	  if (includePreviousBaseForIndels) {
	    if (!variants.empty() &&
		variants.back().ref != variants.back().alt) {
	      VariantAllele a =
		VariantAllele("",
			      alternate.substr(altpos, len),
			      refpos + position);
	      variants.back() = variants.back() + a;
	    } else {
	      VariantAllele a =
		VariantAllele(ref.substr(refpos - 1, 1),
			      alternate.substr(altpos - 1, len + 1),
			      refpos + position - 1);
	      variants.push_back(a);
	    }
	  } else {
	    variants.push_back(VariantAllele("",
					     alternate.substr(altpos, len),
					     refpos + position));
	  }
	  altpos += len;
	  break;
	case 'D':
	  if (includePreviousBaseForIndels) {
	    if (!variants.empty() &&
		variants.back().ref != variants.back().alt) {
	      VariantAllele a
		= VariantAllele(ref.substr(refpos, len)
				, "", refpos + position);
	      variants.back() = variants.back() + a;
	      } else {
	      VariantAllele a
		= VariantAllele(ref.substr(refpos - 1, len + 1),
				alternate.substr(altpos - 1, 1),
				refpos + position - 1);
	      variants.push_back(a);
	    }
	  } else {
	    variants.push_back(VariantAllele(ref.substr(refpos, len),
					     "", refpos + position));
	  }
	  refpos += len;
	  break;

	  // zk has added (!variants.empty()) solves the seg fault in
          // vcfstats, but need to test
	case 'M':
	  {
	    for (int i = 0; i < len; ++i) {
	      VariantAllele a
		= VariantAllele(ref.substr(refpos + i, 1),
				alternate.substr(altpos + i, 1),
				(refpos + i + position));
	      if (useMNPs && (!variants.empty()) &&
		  variants.back().ref.size() == variants.back().alt.size()
		  && variants.back().ref != variants.back().alt) {
		  variants.back() = variants.back() + a;
	      } else {
		variants.push_back(a);
	      }
	    }
	  }
	  refpos += len;
	  altpos += len;
	  break;
	case 'S':
	  {
	    refpos += len;
	    altpos += len;
	    break;
	  }
	default:
	  {
	    break;
	  }
	}
      }
    }
    return variantAlleles;
}


} // end namespace vcf
