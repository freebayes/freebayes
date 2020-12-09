#include "SeqLib/BamHeader.h"
#include "SeqLib/BamWalker.h"

#include <sstream>
#include <stdexcept>
#include <iostream>

#include "htslib/khash.h"

namespace SeqLib {

  BamHeader::BamHeader(const HeaderSequenceVector& hsv) {

    bam_hdr_t * hdr = bam_hdr_init();
    hdr->n_targets = hsv.size();
    hdr->target_len = (uint32_t*)malloc(hdr->n_targets * sizeof(uint32_t));
    hdr->target_name = (char**)malloc(hdr->n_targets * sizeof(char*));
    
    // fill the names and make the text
    std::stringstream text;
    text << "@HD\tVN:1.4" << std::endl;
    for (size_t i = 0; i < hsv.size(); ++i) {
      hdr->target_len[i] = hsv[i].Length;
      hdr->target_name[i] = strdup(hsv[i].Name.c_str());
      text << "@SQ\tSN:" << hsv[i].Name << "\tLN:" << hsv[i].Length << std::endl;
    }
    hdr->text = strdup(text.str().c_str());

    // give to object
    h = SeqPointer<bam_hdr_t>(hdr, bam_hdr_delete()); 
    ConstructName2IDTable();
  }

  int BamHeader::GetSequenceLength(int id) const {
    if (h && id < NumSequences())
      return h->target_len[id];
    return -1;
      
  }

  int BamHeader::GetSequenceLength(const std::string& id) const {
    
    int nid = Name2ID(id);
    if (nid == -1)
      return -1;

    if (h && nid < NumSequences())
      return h->target_len[nid];

    return -1;
  }

BamHeader::BamHeader(const std::string& hdr)  {

  h = SeqPointer<bam_hdr_t>(sam_hdr_read2(hdr), bam_hdr_delete()); 

  ConstructName2IDTable();

}

  std::string BamHeader::AsString() const {
    
    std::stringstream ss;
    
    ss << h->text;
    return ss.str();
    
  }

  BamHeader::BamHeader(const bam_hdr_t * hdr) {

    h = SeqPointer<bam_hdr_t>(bam_hdr_dup(hdr), bam_hdr_delete());

    ConstructName2IDTable();

  }

  void BamHeader::ConstructName2IDTable() {

    // create the lookup table if not already made
    if (!n2i) {
      n2i = SeqPointer<SeqHashMap<std::string, int> >(new SeqHashMap<std::string, int>());
      for (int i = 0; i < h->n_targets; ++i)
	n2i->insert(std::pair<std::string, int>(std::string(h->target_name[i]), i));
    }

  }

  int BamHeader::Name2ID(const std::string& name) const {

    SeqHashMap<std::string, int>::const_iterator ff = n2i->find(name);
    if (ff != n2i->end())
      return ff->second;
    else
      return -1;

  }

bam_hdr_t* BamHeader::sam_hdr_read2(const std::string& hdr) const {

  kstring_t str;
  bam_hdr_t *hhh;
  str.l = str.m = 0; str.s = 0;
  
  std::istringstream iss(hdr);
  std::string line;
  while (std::getline(iss, line, '\n')) {
    //while (hts_getline(fp, KS_SEP_LINE, &fp->line) >= 0) {
    if (line.length() == 0 || line.at(0) != '@') break;
    
    //if (line.length() > 3 && line.substr(0,3) == "@SQ") has_SQ = 1;
    //if (fp->line.l > 3 && strncmp(fp->line.s,"@SQ",3) == 0) has_SQ = 1;
    //kputsn(fp->line.s, fp->line.l, &str);
    kputsn(line.c_str(), line.length(), &str);
    kputc('\n', &str);
  }

  /*
    if (! has_SQ && fp->fn_aux) {
    char line[2048];
    FILE *f = fopen(fp->fn_aux, "r");
    if (f == NULL) return NULL;
    while (fgets(line, sizeof line, f)) {
    const char *name = strtok(line, "\t");
    const char *length = strtok(NULL, "\t");
    ksprintf(&str, "@SQ\tSN:%s\tLN:%s\n", name, length);
    }
    fclose(f);
    }
  */

  if (str.l == 0) 
    kputsn("", 0, &str);
  hhh = sam_hdr_parse(str.l, str.s);
  hhh->l_text = str.l; hhh->text = str.s; // hhh->text needs to be freed
  return hhh;
}

  /** Return the reference sequences as vector of HeaderSequence objects */
  HeaderSequenceVector BamHeader::GetHeaderSequenceVector() const {

    std::vector<HeaderSequence> out;
    for (int i = 0; i < h->n_targets; ++i)
      out.push_back(HeaderSequence(std::string(h->target_name[i]), h->target_len[i]));
    return out;
  }


int BamHeader::NumSequences() const {
  
  if (!h)
    return 0;
  return h->n_targets;

}

std::string BamHeader::IDtoName(int id) const {

  if (id < 0)
    throw std::invalid_argument("BamHeader::IDtoName - ID must be >= 0");

  if (!h)
    throw std::out_of_range("BamHeader::IDtoName - Header is uninitialized");

  if (id >= h->n_targets)
    throw std::out_of_range("BamHeader::IDtoName - Requested ID is higher than number of sequences");
  
  return std::string(h->target_name[id]);

}

  // copied from htslib - sam.c
  /*
  std::string BamHeader::sam_hdr_write2(htsFile *fp, const bam_hdr_t *h)
  {
    switch (fp->format.format) {
    case binary_format:
      fp->format.category = sequence_data;
      fp->format.format = bam;
      // fall-through 
    case bam:
      if (bam_hdr_write(fp->fp.bgzf, h) < 0) return -1;
      break;

    case cram: {
      cram_fd *fd = fp->fp.cram;
      SAM_hdr *hdr = bam_header_to_cram((bam_hdr_t *)h);
      if (! hdr) return -1;
      if (cram_set_header(fd, hdr) < 0) return -1;
      if (fp->fn_aux)
	cram_load_reference(fd, fp->fn_aux);
      if (cram_write_SAM_hdr(fd, fd->header) < 0) return -1;
    }
      break;

    case text_format:
      fp->format.category = sequence_data;
      fp->format.format = sam;
      // fall-through
    case sam: {
      char *p;
      hputs(h->text, fp->fp.hfile);
      p = strstr(h->text, "@SQ\t"); // FIXME: we need a loop to make sure "@SQ\t" does not match something unwanted!!!
      if (p == 0) {
	int i;
	for (i = 0; i < h->n_targets; ++i) {
	  fp->line.l = 0;
	  kputsn("@SQ\tSN:", 7, &fp->line); kputs(h->target_name[i], &fp->line);
	  kputsn("\tLN:", 4, &fp->line); kputw(h->target_len[i], &fp->line); kputc('\n', &fp->line);
	  if ( hwrite(fp->fp.hfile, fp->line.s, fp->line.l) != fp->line.l ) return -1;
	}
      }
      if ( hflush(fp->fp.hfile) != 0 ) return -1;
    }
      break;

    default:
      abort();
    }
    return 0;
  }
  */
}
