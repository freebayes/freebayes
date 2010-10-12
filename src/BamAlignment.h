// ***************************************************************************
// BamAlignment.h (c) 2009 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 9 October 2010 (DB)
// ---------------------------------------------------------------------------
// Provides the BamAlignment data structure
// ***************************************************************************

#ifndef BAMALIGNMENT_H
#define BAMALIGNMENT_H

#include <string>
#include <vector>
#include "BamAux.h"

namespace BamTools {

// BamAlignment data structure
// explicitly labeled as 'struct' to indicate that (most of) its fields are public
struct BamAlignment {

    // constructors & destructor
    public:
        BamAlignment(void);
        BamAlignment(const BamAlignment& other);
        ~BamAlignment(void);

    // Queries against alignment flags
    public:        
        bool IsDuplicate(void) const;           // Returns true if this read is a PCR duplicate       
        bool IsFailedQC(void) const;            // Returns true if this read failed quality control      
        bool IsFirstMate(void) const;           // Returns true if alignment is first mate on read        
        bool IsMapped(void) const;              // Returns true if alignment is mapped        
        bool IsMateMapped(void) const;          // Returns true if alignment's mate is mapped        
        bool IsMateReverseStrand(void) const;   // Returns true if alignment's mate mapped to reverse strand        
        bool IsPaired(void) const;              // Returns true if alignment part of paired-end read        
        bool IsPrimaryAlignment(void) const;    // Returns true if reported position is primary alignment       
        bool IsProperPair(void) const;          // Returns true if alignment is part of read that satisfied paired-end resolution     
        bool IsReverseStrand(void) const;       // Returns true if alignment mapped to reverse strand
        bool IsSecondMate(void) const;          // Returns true if alignment is second mate on read

    // Manipulate alignment flags
    public:        
        void SetIsDuplicate(bool ok);           // Sets "PCR duplicate" flag        
        void SetIsFailedQC(bool ok);            // Sets "failed quality control" flag        
        void SetIsFirstMate(bool ok);           // Sets "alignment is first mate" flag        
        void SetIsMateUnmapped(bool ok);        // Sets "alignment's mate is mapped" flag        
        void SetIsMateReverseStrand(bool ok);   // Sets "alignment's mate mapped to reverse strand" flag        
        void SetIsPaired(bool ok);              // Sets "alignment part of paired-end read" flag        
        void SetIsProperPair(bool ok);          // Sets "alignment is part of read that satisfied paired-end resolution" flag        
        void SetIsReverseStrand(bool ok);       // Sets "alignment mapped to reverse strand" flag        
        void SetIsSecondaryAlignment(bool ok);  // Sets "position is primary alignment" flag        
        void SetIsSecondMate(bool ok);          // Sets "alignment is second mate on read" flag        
        void SetIsUnmapped(bool ok);            // Sets "alignment is mapped" flag

    // Tag data access methods
    public:
        // -------------------------------------------------------------------------------------
        // N.B. - The following tag access methods may not be used on BamAlignments fetched
        // using BamReader::GetNextAlignmentCore().  Attempting to use them will not result in 
        // error message (to keep output clean) but will ALWAYS return false.  Only user-created
        // BamAlignments or those retrieved using BamReader::GetNextAlignment() are valid here.

        // add tag data (create new TAG entry with TYPE and VALUE)
        // TYPE is one of {A, i, f, Z, H} depending on VALUE - see SAM/BAM spec for details
        // returns true if new data added, false if error or TAG already exists
        // N.B. - will NOT modify existing tag. Use EditTag() instead
        // @tag   - two character tag name
        // @type  - single character tag type (see SAM/BAM spec for details)
        // @value - value to associate with tag
        bool AddTag(const std::string& tag, const std::string& type, const std::string& value); // type must be Z or H
        bool AddTag(const std::string& tag, const std::string& type, const uint32_t& value);    // type must be A or i
        bool AddTag(const std::string& tag, const std::string& type, const int32_t& value);     // type must be A or i
        bool AddTag(const std::string& tag, const std::string& type, const float& value);       // type must be A, i, or f
        
        // edit tag data (sets existing TAG with TYPE to VALUE or adds new TAG if not already present)
        // TYPE is one of {A, i, f, Z, H} depending on VALUE - see SAM/BAM spec for details
        // returns true if edit was successfaul, false if error
        // @tag   - two character tag name
        // @type  - single character tag type (see SAM/BAM spec for details)
        // @value - new value for tag
        bool EditTag(const std::string& tag, const std::string& type, const std::string& value); // type must be Z or H
        bool EditTag(const std::string& tag, const std::string& type, const uint32_t& value);    // type must be A or i
        bool EditTag(const std::string& tag, const std::string& type, const int32_t& value);     // type must be A or i
        bool EditTag(const std::string& tag, const std::string& type, const float& value);       // type must be A, i, or f

        // specific tag data access methods - these only remain for legacy support
        // returns whether specific tag could be retrieved
        bool GetEditDistance(uint32_t& editDistance) const; // get "NM" tag data (equivalent to GetTag("NM", editDistance))
        bool GetReadGroup(std::string& readGroup) const;    // get "RG" tag data (equivalent to GetTag("RG", readGroup)) 
        
        // generic tag data access methods 
        // returns whether tag is found & tag type is compatible with DESTINATION
        // @tag - two character tag name
        // @destination - if found, tag value is stored here
        bool GetTag(const std::string& tag, std::string& destination) const;    // access variable-length char or hex strings 
        bool GetTag(const std::string& tag, uint32_t& destination) const;       // access unsigned integer data
        bool GetTag(const std::string& tag, int32_t& destination) const;        // access signed integer data
        bool GetTag(const std::string& tag, float& destination) const;          // access floating point data
        
        // retrieve the tag type code for TAG
        // returns true if tag could be found and type determined
        bool GetTagType(const std::string& tag, char& type) const;
        
        // remove tag data
        // returns true if removal was successful, false if error
        // N.B. - returns false if TAG does not exist (no removal can occur)
        // @tag - two character tag name
        bool RemoveTag(const std::string& tag);

    // Additional data access methods
    public:
        // calculates & returns alignment end position, based on starting position and CIGAR operations
        // @usePadded - if true, counts inserted bases. Default is false, so that alignment end position matches the last base's position in reference
        // @zeroBased - if true, returns 0-based coordinate; else returns 1-based. Setting this to false is useful when using BAM data along with other, half-open formats.
        int GetEndPosition(bool usePadded = false, bool zeroBased = true) const;  

    // 'internal' utility methods 
    private:
        static bool FindTag(const std::string& tag, char* &pTagData, const unsigned int& tagDataLength, unsigned int& numBytesParsed);
        static bool SkipToNextTag(const char storageType, char* &pTagData, unsigned int& numBytesParsed);

    // Data members
    public:
        std::string Name;              // Read name
        int32_t     Length;            // Query length
        std::string QueryBases;        // 'Original' sequence (as reported from sequencing machine)
        std::string AlignedBases;      // 'Aligned' sequence (includes any indels, padding, clipping)
        std::string Qualities;         // FASTQ qualities (ASCII characters, not numeric values)
        std::string TagData;           // Tag data (accessor methods will pull the requested information out)
        int32_t     RefID;             // ID number for reference sequence
        int32_t     Position;          // Position (0-based) where alignment starts
        uint16_t    Bin;               // Bin in BAM file where this alignment resides
        uint16_t    MapQuality;        // Mapping quality score
        uint32_t    AlignmentFlag;     // Alignment bit-flag - see Is<something>() methods to query this value, SetIs<something>() methods to manipulate 
        std::vector<CigarOp> CigarData; // CIGAR operations for this alignment
        int32_t     MateRefID;         // ID number for reference sequence where alignment's mate was aligned
        int32_t     MatePosition;      // Position (0-based) where alignment's mate starts
        int32_t     InsertSize;        // Mate-pair insert size
          
    public:
        struct BamAlignmentSupportData {
      
            // data members
            std::string AllCharData;
            uint32_t    BlockLength;
            uint32_t    NumCigarOperations;
            uint32_t    QueryNameLength;
            uint32_t    QuerySequenceLength;
            bool        HasCoreOnly;
            
            // constructor
            BamAlignmentSupportData(void)
                : BlockLength(0)
                , NumCigarOperations(0)
                , QueryNameLength(0)
                , QuerySequenceLength(0)
                , HasCoreOnly(false)
            { }
        };
        
        // ** THIS IS INTERNAL DATA! DO NOT ACCESS OR EDIT FROM CLIENT CODE **
	//
	// Intended for use by BamReader & BamWriter ONLY. No, really, I mean it.
	//
	// BamTools makes some assumptions about this data being pristine, so please don't tinker with it.
	// The regular data fields above should be sufficient for client code.
	//
	// Technical/design note - Ideally, this would be a private data member with BamReader & BamWriter 
	// allowed direct 'friend' access. However older compilers (especially gcc before v4.1 ) do not 
	// propagate the friend access to BamReader/Writer's implementation inner classes. 
        BamAlignmentSupportData SupportData;   
        
    // Alignment flag query constants
    // Use the get/set methods above instead
    private:
        enum { PAIRED        = 1
             , PROPER_PAIR   = 2
             , UNMAPPED      = 4
             , MATE_UNMAPPED = 8
             , REVERSE       = 16
             , MATE_REVERSE  = 32
             , READ_1        = 64
             , READ_2        = 128
             , SECONDARY     = 256
             , QC_FAILED     = 512
             , DUPLICATE     = 1024 
             };
};

// convenience typedef(s)
typedef std::vector<BamAlignment> BamAlignmentVector;

} // namespace BamTools

#endif // BAMALIGNMENT_H
