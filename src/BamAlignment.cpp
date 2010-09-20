// ***************************************************************************
// BamAlignment.cpp (c) 2009 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 19 September 2010 (DB)
// ---------------------------------------------------------------------------
// Provides the BamAlignment data structure
// ***************************************************************************

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <map>
#include <utility>
#include "BamAlignment.h"
using namespace BamTools;
using namespace std;

// default ctor
BamAlignment::BamAlignment(void) 
    : RefID(-1)
    , Position(-1)
    , MateRefID(-1)
    , MatePosition(-1)
    , InsertSize(0)
{ }

// copy ctor
BamAlignment::BamAlignment(const BamAlignment& other)
    : Name(other.Name)
    , Length(other.Length)
    , QueryBases(other.QueryBases)
    , AlignedBases(other.AlignedBases)
    , Qualities(other.Qualities)
    , TagData(other.TagData)
    , RefID(other.RefID)
    , Position(other.Position)
    , Bin(other.Bin)
    , MapQuality(other.MapQuality)
    , AlignmentFlag(other.AlignmentFlag)
    , CigarData(other.CigarData)
    , MateRefID(other.MateRefID)
    , MatePosition(other.MatePosition)
    , InsertSize(other.InsertSize)
    , SupportData(other.SupportData)
{ }

// dtor
BamAlignment::~BamAlignment(void) { }

// Queries against alignment flags
bool BamAlignment::IsDuplicate(void) const         { return ( (AlignmentFlag & DUPLICATE)     != 0 ); }
bool BamAlignment::IsFailedQC(void) const          { return ( (AlignmentFlag & QC_FAILED)     != 0 ); }
bool BamAlignment::IsFirstMate(void) const         { return ( (AlignmentFlag & READ_1)        != 0 ); }
bool BamAlignment::IsMapped(void) const            { return ( (AlignmentFlag & UNMAPPED)      == 0 ); }
bool BamAlignment::IsMateMapped(void) const        { return ( (AlignmentFlag & MATE_UNMAPPED) == 0 ); }
bool BamAlignment::IsMateReverseStrand(void) const { return ( (AlignmentFlag & MATE_REVERSE)  != 0 ); }
bool BamAlignment::IsPaired(void) const            { return ( (AlignmentFlag & PAIRED)        != 0 ); }
bool BamAlignment::IsPrimaryAlignment(void) const  { return ( (AlignmentFlag & SECONDARY)     == 0 ); }
bool BamAlignment::IsProperPair(void) const        { return ( (AlignmentFlag & PROPER_PAIR)   != 0 ); }
bool BamAlignment::IsReverseStrand(void) const     { return ( (AlignmentFlag & REVERSE)       != 0 ); }
bool BamAlignment::IsSecondMate(void) const        { return ( (AlignmentFlag & READ_2)        != 0 ); }

// Manipulate alignment flags 
void BamAlignment::SetIsDuplicate(bool ok)          { if (ok) AlignmentFlag |= DUPLICATE;     else AlignmentFlag &= ~DUPLICATE; }
void BamAlignment::SetIsFailedQC(bool ok)           { if (ok) AlignmentFlag |= QC_FAILED;     else AlignmentFlag &= ~QC_FAILED; }
void BamAlignment::SetIsFirstMate(bool ok)          { if (ok) AlignmentFlag |= READ_1;        else AlignmentFlag &= ~READ_1; }
void BamAlignment::SetIsMateUnmapped(bool ok)       { if (ok) AlignmentFlag |= MATE_UNMAPPED; else AlignmentFlag &= ~MATE_UNMAPPED; }
void BamAlignment::SetIsMateReverseStrand(bool ok)  { if (ok) AlignmentFlag |= MATE_REVERSE;  else AlignmentFlag &= ~MATE_REVERSE; }
void BamAlignment::SetIsPaired(bool ok)             { if (ok) AlignmentFlag |= PAIRED;        else AlignmentFlag &= ~PAIRED; }
void BamAlignment::SetIsProperPair(bool ok)         { if (ok) AlignmentFlag |= PROPER_PAIR;   else AlignmentFlag &= ~PROPER_PAIR; }
void BamAlignment::SetIsReverseStrand(bool ok)      { if (ok) AlignmentFlag |= REVERSE;       else AlignmentFlag &= ~REVERSE; }
void BamAlignment::SetIsSecondaryAlignment(bool ok) { if (ok) AlignmentFlag |= SECONDARY;     else AlignmentFlag &= ~SECONDARY; }
void BamAlignment::SetIsSecondMate(bool ok)         { if (ok) AlignmentFlag |= READ_2;        else AlignmentFlag &= ~READ_2; }
void BamAlignment::SetIsUnmapped(bool ok)           { if (ok) AlignmentFlag |= UNMAPPED;      else AlignmentFlag &= ~UNMAPPED; }

// calculates alignment end position, based on starting position and CIGAR operations
int BamAlignment::GetEndPosition(bool usePadded, bool zeroBased) const {

    // initialize alignment end to starting position
    int alignEnd = Position;

    // iterate over cigar operations
    vector<CigarOp>::const_iterator cigarIter = CigarData.begin();
    vector<CigarOp>::const_iterator cigarEnd  = CigarData.end();
    for ( ; cigarIter != cigarEnd; ++cigarIter) {
        const char cigarType = (*cigarIter).Type;
        if ( cigarType == 'M' || cigarType == 'D' || cigarType == 'N' )
            alignEnd += (*cigarIter).Length;
        else if ( usePadded && cigarType == 'I' )
            alignEnd += (*cigarIter).Length;
    }
    
    // adjust for zeroBased, if necessary
    if (zeroBased) 
        return alignEnd - 1;
    else 
        return alignEnd;
}

bool BamAlignment::AddTag(const string& tag, const string& type, const string& value) {
  
    if ( SupportData.HasCoreOnly ) return false;
    if ( tag.size() != 2 || type.size() != 1 ) return false;
    if ( type != "Z" && type != "H" ) return false;
  
    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    
    // if tag already exists, return false
    // use EditTag explicitly instead
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) return false;
  
    // otherwise, copy tag data to temp buffer
    string newTag = tag + type + value;
    const int newTagDataLength = tagDataLength + newTag.size() + 1; // leave room for null-term
    char originalTagData[newTagDataLength];
    memcpy(originalTagData, TagData.c_str(), tagDataLength + 1);    // '+1' for TagData null-term
    
    // append newTag
    strcat(originalTagData + tagDataLength, newTag.data());  // removes original null-term, appends newTag + null-term
    
    // store temp buffer back in TagData
    const char* newTagData = (const char*)originalTagData;
    TagData.assign(newTagData, newTagDataLength);
    
    // return success
    return true;
}

bool BamAlignment::AddTag(const string& tag, const string& type, const uint32_t& value) {
  
    if ( SupportData.HasCoreOnly ) return false;
    if ( tag.size() != 2 || type.size() != 1 ) return false;
    if ( type == "f" || type == "Z" || type == "H" ) return false;
  
    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    
    // if tag already exists, return false
    // use EditTag explicitly instead
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) return false;
  
    // otherwise, convert value to string
    union { unsigned int value; char valueBuffer[sizeof(unsigned int)]; } un;
    un.value = value;

    // copy original tag data to temp buffer
    string newTag = tag + type;
    const int newTagDataLength = tagDataLength + newTag.size() + 4; // leave room for new integer
    char originalTagData[newTagDataLength];
    memcpy(originalTagData, TagData.c_str(), tagDataLength + 1);    // '+1' for TagData null-term
    
    // append newTag
    strcat(originalTagData + tagDataLength, newTag.data());
    memcpy(originalTagData + tagDataLength + newTag.size(), un.valueBuffer, sizeof(unsigned int));
    
    // store temp buffer back in TagData
    const char* newTagData = (const char*)originalTagData;
    TagData.assign(newTagData, newTagDataLength);
    
    // return success
    return true;
}

bool BamAlignment::AddTag(const string& tag, const string& type, const int32_t& value) {
    return AddTag(tag, type, (const uint32_t&)value);
}

bool BamAlignment::AddTag(const string& tag, const string& type, const float& value) {
  
    if ( SupportData.HasCoreOnly ) return false;
    if ( tag.size() != 2 || type.size() != 1 ) return false;
    if ( type == "Z" || type == "H" ) return false;
  
    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    
    // if tag already exists, return false
    // use EditTag explicitly instead
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) return false;
  
    // otherwise, convert value to string
    union { float value; char valueBuffer[sizeof(float)]; } un;
    un.value = value;

    // copy original tag data to temp buffer
    string newTag = tag + type;
    const int newTagDataLength = tagDataLength + newTag.size() + 4; // leave room for new float
    char originalTagData[newTagDataLength];
    memcpy(originalTagData, TagData.c_str(), tagDataLength + 1);    // '+1' for TagData null-term
    
    // append newTag
    strcat(originalTagData + tagDataLength, newTag.data());
    memcpy(originalTagData + tagDataLength + newTag.size(), un.valueBuffer, sizeof(float));
    
    // store temp buffer back in TagData
    const char* newTagData = (const char*)originalTagData;
    TagData.assign(newTagData, newTagDataLength);
    
    // return success
    return true;
}

bool BamAlignment::EditTag(const string& tag, const string& type, const string& value) {
  
    if ( SupportData.HasCoreOnly ) return false;
    if ( tag.size() != 2 || type.size() != 1 ) return false;
    if ( type != "Z" && type != "H" ) return false;
  
    // localize the tag data
    char* pOriginalTagData = (char*)TagData.data();
    char* pTagData = pOriginalTagData;
    const unsigned int originalTagDataLength = TagData.size();
    
    unsigned int newTagDataLength = 0;
    unsigned int numBytesParsed = 0;
    
    // if tag found, store data in readGroup, return success
    if ( FindTag(tag, pTagData, originalTagDataLength, numBytesParsed) ) {
        
        // make sure array is more than big enough
        char newTagData[originalTagDataLength + value.size()];  

        // copy original tag data up til desired tag
        const unsigned int beginningTagDataLength = numBytesParsed;
        newTagDataLength += beginningTagDataLength;
        memcpy(newTagData, pOriginalTagData, numBytesParsed);
      
        // copy new VALUE in place of current tag data
        const unsigned int dataLength = strlen(value.c_str());
        memcpy(newTagData + beginningTagDataLength, (char*)value.c_str(), dataLength+1 );
        
        // skip to next tag (if tag for removal is last, return true) 
        const char* pTagStorageType = pTagData - 1;
        if ( !SkipToNextTag(*pTagStorageType, pTagData, numBytesParsed) ) return true;
         
        // copy everything from current tag (the next one after tag for removal) to end
        const unsigned int skippedDataLength = (numBytesParsed - beginningTagDataLength);
        const unsigned int endTagOffset      = beginningTagDataLength + dataLength + 1;
        const unsigned int endTagDataLength  = originalTagDataLength - beginningTagDataLength - skippedDataLength;
        memcpy(newTagData + endTagOffset, pTagData, endTagDataLength);
        
        // ensure null-terminator
        newTagData[ endTagOffset + endTagDataLength + 1 ] = 0;
        
        // save new tag data
        TagData.assign(newTagData, endTagOffset + endTagDataLength);
        return true;
    }
    
    // tag not found, attempt AddTag
    else return AddTag(tag, type, value);
}

bool BamAlignment::EditTag(const string& tag, const string& type, const uint32_t& value) {
  
    if ( SupportData.HasCoreOnly ) return false;
    if ( tag.size() != 2 || type.size() != 1 ) return false;
    if ( type == "f" || type == "Z" || type == "H" ) return false;
    
     // localize the tag data
    char* pOriginalTagData = (char*)TagData.data();
    char* pTagData = pOriginalTagData;
    const unsigned int originalTagDataLength = TagData.size();
    
    unsigned int newTagDataLength = 0;
    unsigned int numBytesParsed = 0;
    
    // if tag found, store data in readGroup, return success
    if ( FindTag(tag, pTagData, originalTagDataLength, numBytesParsed) ) {
        
        // make sure array is more than big enough
        char newTagData[originalTagDataLength + sizeof(value)];  

        // copy original tag data up til desired tag
        const unsigned int beginningTagDataLength = numBytesParsed;
        newTagDataLength += beginningTagDataLength;
        memcpy(newTagData, pOriginalTagData, numBytesParsed);
      
        // copy new VALUE in place of current tag data
        union { unsigned int value; char valueBuffer[sizeof(unsigned int)]; } un;
        un.value = value;
        memcpy(newTagData + beginningTagDataLength, un.valueBuffer, sizeof(unsigned int));
        
        // skip to next tag (if tag for removal is last, return true) 
        const char* pTagStorageType = pTagData - 1;
        if ( !SkipToNextTag(*pTagStorageType, pTagData, numBytesParsed) ) return true;
         
        // copy everything from current tag (the next one after tag for removal) to end
        const unsigned int skippedDataLength = (numBytesParsed - beginningTagDataLength);
        const unsigned int endTagOffset      = beginningTagDataLength + sizeof(unsigned int);
        const unsigned int endTagDataLength  = originalTagDataLength - beginningTagDataLength - skippedDataLength;
        memcpy(newTagData + endTagOffset, pTagData, endTagDataLength);
        
        // ensure null-terminator
        newTagData[ endTagOffset + endTagDataLength + 1 ] = 0;
        
        // save new tag data
        TagData.assign(newTagData, endTagOffset + endTagDataLength);
        return true;
    }
    
    // tag not found, attempt AddTag
    else return AddTag(tag, type, value);
}

bool BamAlignment::EditTag(const string& tag, const string& type, const int32_t& value) {
    return EditTag(tag, type, (const uint32_t&)value);
}

bool BamAlignment::EditTag(const string& tag, const string& type, const float& value) {
  
    if ( SupportData.HasCoreOnly ) return false;
    if ( tag.size() != 2 || type.size() != 1 ) return false;
    if ( type == "Z" || type == "H" ) return false;
    
     // localize the tag data
    char* pOriginalTagData = (char*)TagData.data();
    char* pTagData = pOriginalTagData;
    const unsigned int originalTagDataLength = TagData.size();
    
    unsigned int newTagDataLength = 0;
    unsigned int numBytesParsed = 0;
    
    // if tag found, store data in readGroup, return success
    if ( FindTag(tag, pTagData, originalTagDataLength, numBytesParsed) ) {
        
        // make sure array is more than big enough
        char newTagData[originalTagDataLength + sizeof(value)];  

        // copy original tag data up til desired tag
        const unsigned int beginningTagDataLength = numBytesParsed;
        newTagDataLength += beginningTagDataLength;
        memcpy(newTagData, pOriginalTagData, numBytesParsed);
      
        // copy new VALUE in place of current tag data
        union { float value; char valueBuffer[sizeof(float)]; } un;
        un.value = value;
        memcpy(newTagData + beginningTagDataLength, un.valueBuffer, sizeof(float));
        
        // skip to next tag (if tag for removal is last, return true) 
        const char* pTagStorageType = pTagData - 1;
        if ( !SkipToNextTag(*pTagStorageType, pTagData, numBytesParsed) ) return true;
         
        // copy everything from current tag (the next one after tag for removal) to end
        const unsigned int skippedDataLength = (numBytesParsed - beginningTagDataLength);
        const unsigned int endTagOffset      = beginningTagDataLength + sizeof(float);
        const unsigned int endTagDataLength  = originalTagDataLength - beginningTagDataLength - skippedDataLength;
        memcpy(newTagData + endTagOffset, pTagData, endTagDataLength);
        
        // ensure null-terminator
        newTagData[ endTagOffset + endTagDataLength + 1 ] = 0;
        
        // save new tag data
        TagData.assign(newTagData, endTagOffset + endTagDataLength);
        return true;
    }
    
    // tag not found, attempt AddTag
    else return AddTag(tag, type, value);
}

// get "NM" tag data - originally contributed by Aaron Quinlan
// stores data in 'editDistance', returns success/fail
bool BamAlignment::GetEditDistance(uint32_t& editDistance) const { 
    return GetTag("NM", (uint32_t&)editDistance);
}

// get "RG" tag data
// stores data in 'readGroup', returns success/fail
bool BamAlignment::GetReadGroup(string& readGroup) const {
    return GetTag("RG", readGroup);
}

bool BamAlignment::GetTag(const string& tag, string& destination) const {

    // make sure tag data exists
    if ( SupportData.HasCoreOnly || TagData.empty() ) 
        return false;

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    
    // if tag found, store data in readGroup, return success
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
        const unsigned int dataLength = strlen(pTagData);
        destination.clear();
        destination.resize(dataLength);
        memcpy( (char*)destination.data(), pTagData, dataLength );
        return true;
    }
    
    // tag not found, return failure
    return false;
}

bool BamAlignment::GetTag(const string& tag, uint32_t& destination) const {
  
    // make sure tag data exists
    if ( SupportData.HasCoreOnly || TagData.empty() ) 
        return false;

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    
    // if tag found, determine data byte-length, store data in readGroup, return success
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
        
        // determine data byte-length
        const char type = *(pTagData - 1);
        int destinationLength = 0;
        switch (type) {
            // 1 byte data
            case 'A':
            case 'c':
            case 'C':
                destinationLength = 1;
                break;

            // 2 byte data
            case 's':
            case 'S':
                destinationLength = 2;
                break;

            // 4 byte data
            case 'i':
            case 'I':
                destinationLength = 4;
                break;

            // unsupported type for integer destination (float or var-length strings)
            case 'f':
            case 'Z':
            case 'H':
                fprintf(stderr, "ERROR: Cannot store tag of type %c in integer destination\n", type);
                return false;

            // unknown tag type
            default:
                fprintf(stderr, "ERROR: Unknown tag storage class encountered: [%c]\n", type);
                return false;
        }
          
        // store in destination
        destination = 0;
        memcpy(&destination, pTagData, destinationLength);
        return true;
    }
    
    // tag not found, return failure
    return false;
}

bool BamAlignment::GetTag(const string& tag, int32_t& destination) const {
    return GetTag(tag, (uint32_t&)destination);
}

bool BamAlignment::GetTag(const string& tag, float& destination) const {
  
    // make sure tag data exists
    if ( SupportData.HasCoreOnly || TagData.empty() ) 
        return false;

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    
    // if tag found, determine data byte-length, store data in readGroup, return success
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
        
        // determine data byte-length
        const char type = *(pTagData - 1);
        int destinationLength = 0;
        switch(type) {

            // 1 byte data
            case 'A':
            case 'c':
            case 'C':
                destinationLength = 1;
                break;

            // 2 byte data
            case 's':
            case 'S':
                destinationLength = 2;
                break;

            // 4 byte data
            case 'f':
            case 'i':
            case 'I':
                destinationLength = 4;
                break;
            
            // unsupported type (var-length strings)
            case 'Z':
            case 'H':
                fprintf(stderr, "ERROR: Cannot store tag of type %c in integer destination\n", type);
                return false;

            // unknown tag type
            default:
                fprintf(stderr, "ERROR: Unknown tag storage class encountered: [%c]\n", type);
                return false;
        }
          
        // store in destination
        destination = 0.0;
        memcpy(&destination, pTagData, destinationLength);
        return true;
    }
    
    // tag not found, return failure
    return false;
}

bool BamAlignment::GetTagType(const string& tag, char& type) const {
  
    // make sure tag data exists
    if ( SupportData.HasCoreOnly || TagData.empty() ) 
        return false;

    // localize the tag data
    char* pTagData = (char*)TagData.data();
    const unsigned int tagDataLength = TagData.size();
    unsigned int numBytesParsed = 0;
    
    // lookup tag
    if ( FindTag(tag, pTagData, tagDataLength, numBytesParsed) ) {
        
        // retrieve tag type code
        type = *(pTagData - 1);
        
        // validate that type is a proper BAM tag type
        switch(type) {
            case 'A':
            case 'c':
            case 'C':
            case 's':
            case 'S':
            case 'f':
            case 'i':
            case 'I':
            case 'Z':
            case 'H':
                return true;

            // unknown tag type
            default:
                fprintf(stderr, "ERROR: Unknown tag storage class encountered: [%c]\n", type);
                return false;
        }
    }
    
    // tag not found, return failure
    return false;
}

bool BamAlignment::RemoveTag(const string& tag) {
  
    // BamAlignments fetched using BamReader::GetNextAlignmentCore() are not allowed
    // also, return false if no data present to remove
    if ( SupportData.HasCoreOnly || TagData.empty() ) return false;
  
    // localize the tag data
    char* pOriginalTagData = (char*)TagData.data();
    char* pTagData = pOriginalTagData;
    const unsigned int originalTagDataLength = TagData.size();
    unsigned int newTagDataLength = 0;
    unsigned int numBytesParsed = 0;
    
    // if tag found, store data in readGroup, return success
    if ( FindTag(tag, pTagData, originalTagDataLength, numBytesParsed) ) {
        
        char newTagData[originalTagDataLength];

        // copy original tag data up til desired tag
        pTagData -= 3;
        numBytesParsed -= 3;
        const unsigned int beginningTagDataLength = numBytesParsed;
        newTagDataLength += beginningTagDataLength;
        memcpy(newTagData, pOriginalTagData, numBytesParsed);
        
        // skip to next tag (if tag for removal is last, return true) 
        const char* pTagStorageType = pTagData + 2;
        pTagData       += 3;
        numBytesParsed += 3;
        if ( !SkipToNextTag(*pTagStorageType, pTagData, numBytesParsed) ) return true;
         
        // copy everything from current tag (the next one after tag for removal) to end
        const unsigned int skippedDataLength = (numBytesParsed - beginningTagDataLength);
        const unsigned int endTagDataLength = originalTagDataLength - beginningTagDataLength - skippedDataLength;
        memcpy(newTagData + beginningTagDataLength, pTagData, endTagDataLength );
        
        // save new tag data
        TagData.assign(newTagData, beginningTagDataLength + endTagDataLength);
        return true;
    }
    
    // tag not found, no removal - return failure
    return false;
}

bool BamAlignment::FindTag(const string& tag, char* &pTagData, const unsigned int& tagDataLength, unsigned int& numBytesParsed) {

    while ( numBytesParsed < tagDataLength ) {

        const char* pTagType        = pTagData;
        const char* pTagStorageType = pTagData + 2;
        pTagData       += 3;
        numBytesParsed += 3;

        // check the current tag, return true on match
        if ( strncmp(pTagType, tag.c_str(), 2) == 0 ) 
            return true;

        // get the storage class and find the next tag
        if ( *pTagStorageType == '\0' ) return false; 
        if ( !SkipToNextTag(*pTagStorageType, pTagData, numBytesParsed) ) return false;
        if ( *pTagData == '\0' ) return false;
    }
  
    // checked all tags, none match
    return false;
}

bool BamAlignment::SkipToNextTag(const char storageType, char* &pTagData, unsigned int& numBytesParsed) {
    
    switch(storageType) {

        case 'A':
        case 'c':
        case 'C':
            ++numBytesParsed;
            ++pTagData;
            break;

        case 's':
        case 'S':
            numBytesParsed += 2;
            pTagData       += 2;
            break;

        case 'f':
        case 'i':
        case 'I':
            numBytesParsed += 4;
            pTagData       += 4;
            break;

        case 'Z':
        case 'H':
            while(*pTagData) {
                ++numBytesParsed;
                ++pTagData;
            }
            // increment for null-terminator
            ++numBytesParsed;
            ++pTagData;
            break;

        default: 
            // error case
            fprintf(stderr, "ERROR: Unknown tag storage class encountered: [%c]\n", storageType);
            return false;
    }
    
    // return success
    return true;
}