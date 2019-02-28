/*****************************************************************************
  getcorpus.h

  2019 - Ross Altman
  Inari Agriculture
  raltman@inari.com
******************************************************************************/
#ifndef GETCORPUS_H
#define GETCORPUS_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <map>
#include <math.h>
#include <algorithm>
#include <functional>
#include <thread>
#include <stdio.h>
#include <stdlib.h>
#include <libgen.h>
#include <linux/limits.h>
#include <sys/stat.h>
#include <unistd.h>
#include "gzstream.h"
#include "Fasta.h"
#include "bedFileCorpus.h"

using namespace std;


static unsigned int MAX_THREADS = thread::hardware_concurrency();
static string tmp_dir = "/tmp/corpustool/";


// Helper methods

string padZeros(int val, int threads);

int getFileSize(ifstream *file);

int chunkFile(string inFileName, string chunkName, int threads, unsigned long chunkSize);
void joinFile(string chunkName, string outFileName, int threads, bool gzip_flag);

string translate(string token, string alphabet, string complement);
string tokenizeSeq(string seq, map<string,string> &rev_comp, int k, int stride,
                   string alphabet, string complement, char sep, bool rcomp_flag);
void tokenizeFile(string inFileName, int k, int stride, string alphabet, string complement,
                  char sep, bool rcomp_flag);

void popFreqs(string inFileName, map<string,int> &freq, char sep);

// Indexes the tokens in a corpus tsv file.
map<string,int> buildIndex(string outFileName, map<string,int> freq);
string indexTokens(string tokens, map<string,int> index, char sep);
void indexFile(string inFileName, map<string,int> index, char sep);
map<string,int> readIndex(string inFileName);

template <class T>
void broadcast(T func, string chunkName, int threads);

template <class T>
void broadcast(T func, string chunkName, vector<map<string,int>> &freqs, int threads);

void removeTmpfiles(string chunkName, int threads);


// Class definitions


class Bed2Corpus {
    private:
        // input files
        string _fastaInFile;
        string _indexInFile;
        string _bedFile;
        string _corpusInFile;

        // output files
        string _corpusOutFile;
        string _indexOutFile;

        // options
        string _chunkName;
        string _alphabet;
        string _complement;

        bool _verbose_flag;
        bool _rcomp_flag;
        bool _gzip_flag;

        char _sep;

        int _k;
        int _stride;
        int _threads;

        vector<map<string,int>> _freqs;
        map<string,int> _index;


    public:
        Bed2Corpus(int k, string fastaInFile, string bedFile, int stride,
                   string alphabet, string complement, char sep, int threads,
                   bool verbose_flag, bool rcomp_flag, bool gzip_flag) : _fastaInFile(fastaInFile),
                                                                         _bedFile(bedFile),
                                                                         _alphabet(alphabet),
                                                                         _complement(complement),
                                                                         _verbose_flag(verbose_flag),
                                                                         _rcomp_flag(rcomp_flag),
                                                                         _gzip_flag(gzip_flag), _sep(sep),
                                                                         _k(k), _stride(stride),
                                                                         _threads(threads) {}

        Bed2Corpus(string corpusInFile, string indexInFile, string corpusOutFile,
                   char sep, int threads, bool verbose_flag, bool rcomp_flag,
                   bool gzip_flag) : _indexInFile(indexInFile), _corpusInFile(corpusInFile),
                                                                _corpusOutFile(corpusOutFile),
                                                                _verbose_flag(verbose_flag), _rcomp_flag(rcomp_flag),
                                                                _gzip_flag(gzip_flag), _sep(sep),
                                                                _threads(threads) {}

        ~Bed2Corpus(void);

        void NewCorpusOnly(string corpusOutFile);

        void NewIndexNewCorpus(string corpusOutFile, string indexOutFile);

        void ExistingIndexNewCorpus(string indexInFile, string corpusOutFile);

        void NewIndexOnly(string indexOutFile);

        void ExistingIndexExistingCorpus(void);

        // Intermediate methods
        void doChunkTokenize(void);
        void doNewIndex(string indexOutFile);
        void doIndexCorpus(string corpusOutFile);
};


class fileTokenizer {
    private:
        string _alphabet;
        string _complement;
        char _sep;
        int _k;
        int _stride;
        bool _rcomp_flag;

    public:
        fileTokenizer(int k, int stride, string alphabet, string complement,
                      char sep, bool rcomp_flag) : _alphabet(alphabet), _complement(complement),
                                                     _sep(sep), _k(k), _stride(stride),
                                                     _rcomp_flag(rcomp_flag) {}
        void operator()(string inFileName) { 
            tokenizeFile(inFileName, _k, _stride, _alphabet, _complement, _sep, _rcomp_flag);
        }
};


class freqPopulator {
    private:
        char _sep;
    public:
        freqPopulator(char sep) : _sep(sep) {}
        void operator()(string inFileName, map<string,int> &freq) { 
            popFreqs(inFileName, freq, _sep);
        }
};


class fileIndexer {
    private:
        map<string,int> _index;
        char _sep;
    public:
        fileIndexer(map<string,int> index, char sep) : _index(index), _sep(sep) {}
        void operator()(string inFileName) {
            indexFile(inFileName, _index, _sep);
        }
};


#endif