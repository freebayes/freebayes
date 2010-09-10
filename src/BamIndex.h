// ***************************************************************************
// BamIndex.h (c) 2009 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 10 September 2010 (DB)
// ---------------------------------------------------------------------------
// Provides index functionality - both for the standardized BAM index format
// (".bai") as well as a BamTools-specific (nonstandard) index format (".bti").
// ***************************************************************************

#ifndef BAM_INDEX_H
#define BAM_INDEX_H

#include <iostream>
#include <string>
#include <vector>
#include "BamAux.h"

namespace BamTools {

class BamReader;
class BgzfData;
  
// --------------------------------------------------  
// BamIndex base class
class BamIndex {

    // ctor & dtor
    public:
        BamIndex(BamTools::BgzfData* bgzf, BamTools::BamReader* reader, bool isBigEndian);
        virtual ~BamIndex(void) { }
        
    // index interface
    public:
        // creates index data (in-memory) from current reader data
        virtual bool Build(void) =0;
        // returns supported file extension
        virtual const std::string Extension(void) const =0;
        // attempts to use index to jump to region; returns success/fail
        virtual bool Jump(const BamTools::BamRegion& region) =0;
        // returns whether reference has alignments or no
        virtual bool HasAlignments(const int& referenceID); 
        // loads existing data from file into memory
        virtual bool Load(const std::string& filename)  =0;
        // writes in-memory index data out to file 
        // N.B. - (this is the original BAM filename, method will modify it to use applicable extension)
        virtual bool Write(const std::string& bamFilename) =0;
        
    // factory methods for returning proper BamIndex-derived type based on available index files
    public:
      
        // returns index based on BAM filename 'stub'
        // checks first for preferred type, returns that type if found
        // (if not found, attmempts to load other type(s), returns 0 if NONE found)
        //
        // ** default preferred type is BamToolsIndex ** use this anytime it exists
         enum PreferredIndexType { BAMTOOLS = 0, STANDARD };
        static BamIndex* FromBamFilename(const std::string& bamFilename, 
                                         BamTools::BgzfData* bgzf, 
                                         BamTools::BamReader* reader, 
                                         bool isBigEndian, 
                                         const BamIndex::PreferredIndexType& type = BamIndex::BAMTOOLS);
        
        // returns index based on explicitly named index file (or 0 if not found)
        static BamIndex* FromIndexFilename(const std::string& indexFilename, 
                                           BamTools::BgzfData* bgzf, 
                                           BamTools::BamReader* reader, 
                                           bool isBigEndian);

    // data members
    protected:
        BamTools::BgzfData*  m_BGZF;
        BamTools::BamReader* m_reader;
        BamTools::RefVector  m_references;
        bool m_isBigEndian;
};

// --------------------------------------------------
// BamStandardIndex class
// 
// implements standardized (per SAM/BAM spec) index file ops
class BamStandardIndex : public BamIndex {

  
    // ctor & dtor
    public:
        BamStandardIndex(BamTools::BgzfData*  bgzf, 
                        BamTools::BamReader* reader,
                        bool isBigEndian);
        ~BamStandardIndex(void);
        
    // interface (implements BamIndex virtual methods)
    public:
        // creates index data (in-memory) from current reader data
        bool Build(void);
        // returns supported file extension
        const std::string Extension(void) const { return std::string(".bai"); }
        // attempts to use index to jump to region; returns success/fail
        bool Jump(const BamTools::BamRegion& region);
         // loads existing data from file into memory
        bool Load(const std::string& filename);
        // writes in-memory index data out to file 
        // N.B. - (this is the original BAM filename, method will modify it to use applicable extension)
        bool Write(const std::string& bamFilename);
      
    // internal implementation
    private:
        struct BamStandardIndexPrivate;
        BamStandardIndexPrivate* d;
};

// --------------------------------------------------
// BamToolsIndex class
//
// implements BamTools-specific index file ops
class BamToolsIndex : public BamIndex {

    // ctor & dtor
    public:
        BamToolsIndex(BamTools::BgzfData*  bgzf, 
                      BamTools::BamReader* reader,
                      bool isBigEndian);
        ~BamToolsIndex(void);
        
    // interface (implements BamIndex virtual methods)
    public:
        // creates index data (in-memory) from current reader data
        bool Build(void);
        // returns supported file extension
        const std::string Extension(void) const { return std::string(".bti"); }
        // attempts to use index to jump to region; returns success/fail
        bool Jump(const BamTools::BamRegion& region);
         // loads existing data from file into memory
        bool Load(const std::string& filename);
        // writes in-memory index data out to file 
        // N.B. - (this is the original BAM filename, method will modify it to use applicable extension)
        bool Write(const std::string& bamFilename);
    
    // internal implementation
    private:
        struct BamToolsIndexPrivate;
        BamToolsIndexPrivate* d;
};

// --------------------------------------------------
// BamIndex factory methods
//
// return proper BamIndex-derived type based on available index files

inline
BamIndex* BamIndex::FromBamFilename(const std::string& bamFilename, 
                                    BamTools::BgzfData* bgzf, 
                                    BamTools::BamReader* reader, 
                                    bool isBigEndian, 
                                    const BamIndex::PreferredIndexType& type)
{
    // ---------------------------------------------------
    // attempt to load preferred type first
    
    const std::string bamtoolsIndexFilename = bamFilename + ".bti";
    const bool bamtoolsIndexExists = BamTools::FileExists(bamtoolsIndexFilename);
    if ( (type == BamIndex::BAMTOOLS) && bamtoolsIndexExists )
        return new BamToolsIndex(bgzf, reader, isBigEndian);

    const std::string standardIndexFilename = bamFilename + ".bai";
    const bool standardIndexExists = BamTools::FileExists(standardIndexFilename);
    if ( (type == BamIndex::STANDARD) && standardIndexExists ) 
        return new BamStandardIndex(bgzf, reader, isBigEndian);
    
    // ----------------------------------------------------
    // preferred type could not be found, try other (non-preferred) types
    // if none found, return 0
    
    if ( bamtoolsIndexExists ) return new BamToolsIndex(bgzf, reader, isBigEndian);
    if ( standardIndexExists ) return new BamStandardIndex(bgzf, reader, isBigEndian);
    return 0;
}

inline
BamIndex* BamIndex::FromIndexFilename(const std::string& indexFilename,
                                      BamTools::BgzfData* bgzf, 
                                      BamTools::BamReader* reader, 
                                      bool isBigEndian) 
{
    // see if specified file exists
    const bool indexExists = BamTools::FileExists(indexFilename);
    if ( !indexExists ) return 0;
  
    const std::string bamtoolsIndexExtension(".bti");
    const std::string standardIndexExtension(".bai");
  
    // if has bamtoolsIndexExtension
    if ( indexFilename.find(bamtoolsIndexExtension) == (indexFilename.length() - bamtoolsIndexExtension.length()) )
        return new BamToolsIndex(bgzf, reader, isBigEndian);
    
     // if has standardIndexExtension
    if ( indexFilename.find(standardIndexExtension) == (indexFilename.length() - standardIndexExtension.length()) )
        return new BamStandardIndex(bgzf, reader, isBigEndian);

    // otherwise, unsupported file type
    return 0;
}

} // namespace BamTools

#endif // BAM_INDEX_H
