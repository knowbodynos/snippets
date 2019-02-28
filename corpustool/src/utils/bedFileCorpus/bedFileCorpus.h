/*****************************************************************************
  bedFileCorpus.h

  2019 - Ross Altman
  Inari Agriculture
  raltman@inari.com
******************************************************************************/
#ifndef BEDFILECORPUS_H
#define BEDFILECORPUS_H

// "local" includes
#include "bedFile.h"
#include "Fasta.h"

using namespace std;


//************************************************
// BedFileCorpus Class methods and elements
//************************************************
class BedFileCorpus : public BedFile {
    using BedFile::BedFile;

    public:

        int GetChunkSize(FastaReference *fr, unsigned int threads);

        void GetNextSeq(FastaReference *fr, BED &bedEntry, string &sequence);

        int ChunkFile(FastaReference *fr, string chunkName, int threads, int chunkSize);

};

#endif /* BEDFILECORPUS_H */
