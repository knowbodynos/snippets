/*****************************************************************************
  bedFileCorpus.cpp

  2019 - Ross Altman
  Inari Agriculture
  raltman@inari.com
******************************************************************************/
#include "bedFileCorpus.h"
#include "Fasta.h"
#include <math.h>


int BedFileCorpus::GetChunkSize(FastaReference *fr, unsigned int threads) {
    BED bedEntry, nullBed;
    Open();
    int fullLength = 0;
    while (GetNextBed(bedEntry)) {
        if (_status == BED_VALID) {
            // make sure we are extracting >= 1 bp
            if (bedEntry.zeroLength == false) {

                size_t seqLength = fr->sequenceLength(bedEntry.chrom);
                // seqLength > 0 means chrom was found in index.
                // seqLength == 0 otherwise.
                if (seqLength) {
                    // make sure this feature will not exceed 
                    // the end of the chromosome.
                    if ((bedEntry.start <= seqLength) && (bedEntry.end <= seqLength)) {
                        int length = bedEntry.end - bedEntry.start;
                        fullLength += length;
                    } else {
                        cerr << "Feature (" << bedEntry.chrom << ":" 
                             << bedEntry.start << "-" << bedEntry.end 
                             << ") beyond the length of "
                             << bedEntry.chrom 
                             << " size (" << seqLength << " bp).  Skipping." 
                             << endl;
                    }
                } else {
                    cerr << "WARNING. chromosome (" 
                         << bedEntry.chrom 
                         << ") was not found in the FASTA file. Skipping."
                         << endl;
                }
            } else { // handle zeroLength 
                cerr << "Feature (" << bedEntry.chrom << ":" 
                     << bedEntry.start+1 << "-" << bedEntry.end-1
                     << ") has length = 0, Skipping." 
                     << endl;
            }
            bedEntry = nullBed;
        }
    }
    Close();
    double chunkSizeDouble = ceil((double)fullLength / threads);
    return (int)chunkSizeDouble;
}


void BedFileCorpus::GetNextSeq(FastaReference *fr, BED &bedEntry, string &sequence) {
    BED nullBed;
    GetNextBed(bedEntry);
    if (_status == BED_VALID) {
        // make sure we are extracting >= 1 bp
        if (bedEntry.zeroLength == false) {

            size_t seqLength = fr->sequenceLength(bedEntry.chrom);
            // seqLength > 0 means chrom was found in index.
            // seqLength == 0 otherwise.
            if (seqLength) {
                // make sure this feature will not exceed 
                // the end of the chromosome.
                if ((bedEntry.start <= seqLength) && (bedEntry.end <= seqLength)) {
                    int length = bedEntry.end - bedEntry.start;
                    sequence = fr->getSubSequence(bedEntry.chrom, bedEntry.start, length);
                } else {
                    cerr << "Feature (" << bedEntry.chrom << ":" 
                         << bedEntry.start << "-" << bedEntry.end 
                         << ") beyond the length of "
                         << bedEntry.chrom 
                         << " size (" << seqLength << " bp).  Skipping." 
                         << endl;
                }
            } else {
                cerr << "WARNING. chromosome (" 
                     << bedEntry.chrom 
                     << ") was not found in the FASTA file. Skipping."
                     << endl;
            }
        } else { // handle zeroLength 
            cerr << "Feature (" << bedEntry.chrom << ":" 
                 << bedEntry.start+1 << "-" << bedEntry.end-1
                 << ") has length = 0, Skipping." 
                 << endl;
        }
        // bed = nullBed;
    }
}


// Chunks a file by breaking it up into chunks of "chunkSize" bytes.
int BedFileCorpus::ChunkFile(FastaReference *fr, string chunkName, int threads, int chunkSize) {
    string fullChunkName;
    BED bedEntry;

    Open();

    string sequence;
    int zeros = to_string(threads).length();
    int counter = 1;
    GetNextSeq(fr, bedEntry, sequence);
    while (_status == BED_VALID) { 
        // Build the chunk file name. Usually /tmp/corpustool/chunkName.ext.N
        // N represents the Nth chunk
        fullChunkName.clear();
        fullChunkName.append(chunkName);
        fullChunkName.append(".");

        // Convert counter integer into string and append to name.
        string ext = to_string(counter);
        ext = string(zeros - ext.length(), '0') + ext;
        fullChunkName.append(ext);

        // Open new chunk file name for output
        ofstream output;
        output.open(fullChunkName.c_str(), ios::out | ios::trunc | ios::binary);

        // If chunk file opened successfully, read from input and 
        // write to output chunk. Then close.
        if (output.is_open()) {
            int currChunkSize = 0;
            while ((currChunkSize < chunkSize) && (_status == BED_VALID)) {
                int seqLength = sequence.length();
                if (seqLength <= chunkSize - currChunkSize) {
                    output << sequence << endl;
                    currChunkSize += seqLength;
                    sequence = "";
                    GetNextSeq(fr, bedEntry, sequence);
                } else {
                    output << sequence.substr(0, chunkSize - currChunkSize);
                    sequence = sequence.substr(chunkSize - currChunkSize, seqLength - chunkSize + currChunkSize);
                    currChunkSize += chunkSize - currChunkSize;
                }
            }
            output.close();
            counter++;
        } else {
            cerr << "Error: Unable to open file '" << fullChunkName << "' for writing." << endl;
        }
    }
    return counter - 1;
}
