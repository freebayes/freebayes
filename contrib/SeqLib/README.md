[![Build Status](https://travis-ci.org/walaj/SeqLib.svg?branch=master)](https://travis-ci.org/walaj/SeqLib)
[![Coverage Status](https://coveralls.io/repos/github/walaj/SeqLib/badge.svg?branch=master)](https://coveralls.io/github/walaj/SeqLib?branch=master)

C++ interface to HTSlib, BWA-MEM and Fermi

**License:** [Apache2][license]

API Documentation
-----------------
[API Documentation][htmldoc]

Citation
--------
If you use SeqLib in your applications, please cite: http://bioinformatics.oxfordjournals.org/content/early/2016/12/21/bioinformatics.btw741.full.pdf+html

Note that the values for the SeqAn benchmarking in Table 2 should be corrected to 7.7 Gb memory and 33.92 seconds in CPU time, when compiling SeqAn with ``-O3 -DNDEBUG``. SeqAn also does full string decompression.
Wall times for SeqAn may be shorter than CPU time because it uses embedded multi-threading during BAM IO.

Table of contents
=================

  * [Installation](#installation)
  * [Integrating into build system](#integrating-into-build-system)
  * [Description](#description)
    * [Memory management](#memory-management)
    * [Other C++ APIs](#other-c++-apis)
  * [Command line usage](#command-line-usage)
  * [Examples](#examples)
  * [Attributions](#attributions)

Installation
------------

#######
```bash
git clone --recursive https://github.com/walaj/SeqLib.git
cd SeqLib
## cd htslib && ./configure --enable-libcurl && cd .. # support for remote (FTP/HTTPS/Google etc) BAM access
./configure ## or: ./configure LDFLAGS='-lcurl -lcrypto' # for remote support
make ## for c++11 (req. for AhoCorasick), run as: make CXXFLAGS='-std=c++11'
make install
make seqtools ## for the command line version
```
 
I have successfully compiled with GCC-4.5+ and Clang on Linux and OSX.

SeqLib is compatible with c++98 and later.

Integrating into build system
-----------------------------

After building, you will need to add the relevant header directories:
```bash
SEQ=<path_to_seqlib_git_repos>
C_INCLUDE_PATH=$C_INCLUDE_PATH:$SEQ:$SEQ/htslib
```

And need to link the SeqLib static library and Fermi, BWA and HTSlib libraries
```bash
SEQ=<path_to_seqlib>
LDFLAGS="$LDFLAGS $SEQ/bin/libseqlib.a $SEQ/bin/libbwa.a $SEQ/bin/libfml.a $SEQ/bin/libhts.a"
```

To add support for reading BAMs, etc with HTTPS, FTP, S3, Google cloud, etc, you must compile and link with libcurl.
```bash
## set hts to build with libcurl links and hfile_libcurl.c
cd SeqLib/htslib
./configure --enable-libcurl 
## compile seqlib with libcurl support
cd ../ # back to SeqLib main directory
./configure LDFLAGS="-lcurl -lcrypto"
make 
make install
```
Remember then to then link any projects made with SeqLib with the additional ``-lcurl -lcrypto`` flags.

Description
-----------

SeqLib is a C++ library for querying BAM/SAM/CRAM files, performing 
BWA-MEM operations in memory, and performing sequence assembly. Core operations
in SeqLib are peformed by:
* [HTSlib][htslib]
* [BWA-MEM][BWA] (Apache2 branch)
* [FermiKit][fermi]

The primary developer for these three projects is Heng Li.

SeqLib also has support for storing and manipulating genomic intervals via ``GenomicRegion`` and ``GenomicRegionCollection``. 
It uses an [interval tree][int] (provided by Erik Garrison @ekg) to provide for rapid interval queries.

SeqLib is built to be extendable. See [VariantBam][var] for examples of how to take advantage of C++
class extensions to build off of the SeqLib base functionality. 
 
Memory management
-----------------
SeqLib is built to automatically handle memory management of C code from BWA-MEM and HTSlib by using C++ smart
pointers that handle freeing memory automatically. One of the 
main motivations behind SeqLib is that all access to sequencing reads, BWA, etc should
completely avoid ``malloc`` and ``free``. In SeqLib all the mallocs/frees are handled automatically in the constructors and
destructors.

Other C++ APIs
------------------------------
There are overlaps between this project and [BamTools][BT] from Derek Barnett, [Gamgee][gam] 
from the Broad Institute, and [SeqAn][seqan] from Freie Universitat Berlin. These projects 
provide excellent and high quality APIs. SeqLib provides further performance enhancement and new capabilites for certain classes of 
bioinformatics problems.

Some differences:
* SeqLib has ~2-4x faster read/write speed over BamTools and lower memory footprint.
* SeqLib has support for CRAM file
* SeqLib provides in memory access to BWA-MEM, BLAT, chromosome aware interval tree, read correction, and sequence assembly with Fermi.
* SeqAn provide a substantial amount of additional capabilites not in SeqLib, including graph operations and an expanded suite of multi-sequence alignments.
* SeqAn embeds multi-threading into some functionality like BAM IO to improve wall times.

For your particular application, our hope is that SeqLib will provide a comprehensive and powerful envrionment to develop 
bioinformatics tools, or to be used in conjuction with the capablities in SeqAn and BamTools. Feature requests and comments are welcomed.

Command Line Usage
------------------
```bash
## BFC correction (input mode -m is b (BAM/SAM/CRAM), output mode -w is SAM stream
samtools view in.bam -h 1:1,000,000-1,002,000 | seqtools bfc - -G $REF | samtools sort - -m 4g -o corrected.bam

## Without a pipe, write to BAM
seqtools bfc in.bam -G $REF -b > corrected.bam

## Skip realignment, send to fasta
seqtools bfc in.bam -f > corrected.fasta

## Input as gzipped or plain fasta (or fastq), send to aligned BAM
seqtools bfc --infasta in.fasta -G $REG -b > corrected.bam


##### ASSEMBLY (same patterns as above)
samtools view in.bam -h 1:1,000,000-1,002,000 | seqtools fml - -G $REF | samtools sort - -m 4g -o assembled.bam

```

Examples
--------
##### Targeted re-alignment of reads to a given region with BWA-MEM
```
#include "SeqLib/RefGenome.h"
#include "SeqLib/BWAWrapper.h"
using namespace SeqLib;
RefGenome ref;
ref.LoadIndex("hg19.fasta");

// get sequence at given locus
std::string seq = ref.QueryRegion("1", 1000000,1001000);

// Make an in-memory BWA-MEM index of region
BWAWrapper bwa;
UnalignedSequenceVector usv = {{"chr_reg1", seq}};
bwa.ConstructIndex(usv);

// align an example string with BWA-MEM
std::string querySeq = "CAGCCTCACCCAGGAAAGCAGCTGGGGGTCCACTGGGCTCAGGGAAG";
BamRecordVector results;
// hardclip=false, secondary score cutoff=0.9, max secondary alignments=10
bwa.AlignSequence(querySeq, "my_seq", results, false, 0.9, 10); 

// print results to stdout
for (auto& i : results)
    std::cout << i << std::endl;
```

##### Read a BAM line by line, realign reads with BWA-MEM, write to new BAM
```
#include "SeqLib/BamReader.h"
#include "SeqLib/BWAWrapper.h"
using namespace SeqLib;

// open the reader BAM/SAM/CRAM
BamReader bw;
bw.Open("test.bam");

// open a new interface to BWA-MEM
BWAWrapper bwa;
bwa.LoadIndex("hg19.fasta");

// open the output BAM
BamWriter writer; // or writer(SeqLib::SAM) or writer(SeqLib::CRAM) 
writer.SetWriteHeader(bwa.HeaderFromIndex());
writer.Open("out.bam");
writer.WriteHeader();

BamRecord r;
bool hardclip = false;
float secondary_cutoff = 0.90; // secondary alignments must have score >= 0.9*top_score
int secondary_cap = 10; // max number of secondary alignments to return
while (GetNextRecord(r)) {
      BamRecordVector results; // alignment results (can have multiple alignments)
      bwa.AlignSequence(r.Sequence(), r.Qname(), results, hardclip, secondary_cutoff, secondary_cap);

      for (auto& i : results)
        writer.WriteRecord(i);
}
```


##### Perform sequence assembly with Fermi directly from a BAM
```

#include "SeqLib/FermiAssembler.h"
using namespace SeqLib;

FermiAssembler f;

// read in data from a BAM
BamReader br;
br.Open("test_data/small.bam");

// retreive sequencing reads (up to 20,000)
BamRecord r;
BamRecordVector brv;
size_t count = 0;
while(br.GetNextRead(r) && count++ < 20000) 
  brv.push_back(r);

// add the reads and error correct them  
f.AddReads(brv);
f.CorrectReads();

// peform the assembly
f.PerformAssembly();

// retrieve the contigs
std::vector<std::string> contigs = f.GetContigs();

// write as a fasta to stdout
for (size_t i = 0; i < contigs.size(); ++i)
    std::cout << ">contig" << i << std::endl << contigs[i] << std::endl;
```

##### Plot a collection of gapped alignments
```
using namespace SeqLib;
BamReader r;
r.Open("test_data/small.bam");

GenomicRegion gr("X:1,002,942-1,003,294", r.Header());
r.SetRegion(gr);

SeqPlot s;
s.SetView(gr);

BamRecord rec;
BamRecordVector brv;
while(r.GetNextRecord(rec))
  if (!rec.CountNBases() && rec.MappedFlag())
    brv.push_back(rec);
s.SetPadding(20);

std::cout << s.PlotAlignmentRecords(brv);
```

Trimmed output from above (NA12878):
```
CTATCTATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTGTCCATCCATCCATCCATCCA
CTATCTATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTAT                    CATCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATCTATCAATCCATCCATCCATCCA
 TATCTATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCA              
 TATCTATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTGTCCATCCATCCATCCATC                        
 TATCTATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCA              
  ATCTATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCAT             
   TCTATCTATCTCTTCTTCTGTCCGCTCATGTGTCTGTCCATCTATCTATC                    GTCCATCCATCCATCCATCCATCCATCCACCCATTCATCCATCCACCTATCCATCTATCAATCCATCCATCCATCCATCCGTCTATCTTATGCATCACAGC
   TCTATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATT                                                                  
    CTATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTC                                                                 
    CTATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTC                                                                 
      ATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCAT                                                               
      ATCTATCTCTTCTTCTGTCCGTTCATGTGTCTGTCCATCTATCTATCCATCTATCTATCATCTAACTATCTG----TCCATCCATCCATCCATCCACCCATTCAT                                                               
```

##### Read simultaneously from a BAM, CRAM and SAM file and send to stdout
```
using namespace SeqLib;
#include "SeqLib/BamReader.h"
BamReader r;
  
// read from multiple streams coordinate-sorted order
r.Open("test_data/small.bam");
r.Open("test_data/small.cram");
r.Open("test_data/small.sam");

BamWriter w(SeqLib::SAM); // set uncompressed output
w.Open("-");              // write to stdout
w.SetHeader(r.Header());  // specify the header
w.WriteHeader();          // write out the header

BamRecord rec;
while(r.GetNextRecord(rec)) 
  w.WriteRecord(rec);
w.Close();               // Optional. Will close on destruction
```

##### Perform error correction on reads, using [BFC][bfc]
```
#include "SeqLib/BFC.h"
using namespace SeqLib;

// brv is some set of reads to train the error corrector
for (BamRecordVector::const_iterator r = brv.begin(); r != brv.end(); ++r)
    b.AddSequence(r->Sequence().c_str(), r->Qualities().c_str(), r->Qname().c_str());
b.Train();
b.clear(); // clear the training sequences. Training parameters saved in BFC object

// brv2 is some set to correct
for (BamRecordVector::const_iterator r = brv2.begin(); r != brv2.end(); ++r)
    b.AddSequence(r->Sequence().c_str(), r->Qualities().c_str(), r->Qname().c_str());
b.ErrorCorrect();

// retrieve the sequences
UnalignedSequenceVector v;
std::string name, seq;
while (b.GetSequences(seq, name))
  v.push_back({name, seq});      			   

```

Support
-------
This project is being actively developed and maintained by Jeremiah Wala (jwala@broadinstitute.org). 

Attributions
------------
We would like to thank Heng Li (htslib/bwa/fermi), Erik Garrison (interval tree), Christopher Gilbert (aho corasick), 
and Mengyao Zhao (sw alignment), for providing open-source and robust bioinformatics solutions.

Development, support, guidance, testing:
* Steve Huang - Research Scientist, Broad Institute
* Steve Schumacher - Computational Biologist, Dana Farber Cancer Institute
* Cheng-Zhong Zhang - Research Scientist, Broad Institute
* Marcin Imielinski - Assistant Professor, Cornell University
* Rameen Beroukhim - Assistant Professor, Harvard Medical School

[htslib]: https://github.com/samtools/htslib.git

[SGA]: https://github.com/jts/sga

[BLAT]: https://genome.ucsc.edu/cgi-bin/hgBlat?command=start

[BWA]: https://github.com/lh3/bwa

[license]: https://github.com/walaj/SeqLib/blob/master/LICENSE

[BamTools]: https://raw.githubusercontent.com/wiki/pezmaster31/bamtools/Tutorial_Toolkit_BamTools-1.0.pdf

[API]: http://pezmaster31.github.io/bamtools/annotated.html

[htmldoc]: http://walaj.github.io/seqlibdocs/doxygen

[var]: https://github.com/walaj/VariantBam

[BT]: https://github.com/pezmaster31/bamtools

[seqan]: https://www.seqan.de

[gam]: https://github.com/broadinstitute/gamgee

[int]: https://github.com/ekg/intervaltree.git

[fermi]: https://github.com/lh3/fermi-lite

[bfc]: https://github.com/lh3/bfc
