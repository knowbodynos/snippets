/*****************************************************************************
  getcorpusMain.cpp

  2019 - Ross Altman
  Inari Agriculture
  raltman@inari.com
******************************************************************************/
#include "getcorpus.h"

using namespace std;

// define our program name
#define PROGRAM_NAME "corpustool getcorpus"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)


void getcorpus_help(int code);


int getcorpus_main(int argc, char* argv[]) {
    bool showHelp = false;

    // input files
    string fastaInFile;
    string indexInFile;
    string bedFile;
    string corpusInFile;

    // output files
    string corpusOutFile;
    string indexOutFile;

    // checks for existence of parameters
    bool haveFastaIn = false;
    bool haveIndexIn = false;
    bool haveCorpusIn = false;
    bool haveBed = false;
    bool haveCorpusOut = false;
    bool haveIndexOut = false;
    bool haveK = false;

    bool verbose_flag = true;
    bool rcomp_flag = true;
    bool gzip_flag = false;

    int k = 0;
    int stride = 1;
    int threads = MAX_THREADS;

    char sep = ' ';

    string alphabet = "ACGTNBDHVRYKM";
    string complement = "TGCANVHDBYRMK";
    
    // check to see if we should print out some help
    if (argc <= 1) {
        showHelp = true;
    }

    for (int i=2;i<argc;i++) {
        int parameterLength = (int)strlen(argv[i]);
        if ((PARAMETER_CHECK("-h", 2, parameterLength)) || (PARAMETER_CHECK("--help", 6, parameterLength))) {
            showHelp = true;
        }
    }

    if (showHelp) {
        getcorpus_help(0);
    }

    for (int i=2;i<argc;i++) {
        int parameterLength = (int)strlen(argv[i]);
        if (PARAMETER_CHECK("-fi", 3, parameterLength)) {
            if ((i+1) < argc) {
                haveFastaIn = true;
                fastaInFile = argv[i + 1];
                i++;
            }
        } else if (PARAMETER_CHECK("-ii", 3, parameterLength)) {
            if ((i+1) < argc) {
                haveIndexIn = true;
                indexInFile = argv[i + 1];
                i++;
            }
        } else if (PARAMETER_CHECK("-ci", 3, parameterLength)) {
            if ((i+1) < argc) {
                haveCorpusIn = true;
                corpusInFile = argv[i + 1];
                i++;
            }
        } else if (PARAMETER_CHECK("-bed", 4, parameterLength)) {
            if ((i+1) < argc) {
                haveBed = true;
                bedFile = argv[i + 1];
                i++;
            }
        } else if (PARAMETER_CHECK("-co", 3, parameterLength)) {
            if ((i+1) < argc) {
                haveCorpusOut = true;
                corpusOutFile = argv[i + 1];
                i++;
            }
        } else if (PARAMETER_CHECK("-io", 3, parameterLength)) {
            if ((i+1) < argc) {
                haveIndexOut = true;
                indexOutFile = argv[i + 1];
                i++;
            }
        } else if (PARAMETER_CHECK("-k", 2, parameterLength)) {
            if ((i+1) < argc) {
                k = stoi(argv[i + 1]);
                haveK=true;
                i++;
            }
        } else if ((PARAMETER_CHECK("--stride", 8, parameterLength)) || (PARAMETER_CHECK("-s", 2, parameterLength))) {
            if ((i+1) < argc) {
                stride = stoi(argv[i + 1]);
                i++;
            }
        } else if ((PARAMETER_CHECK("--alphabet", 10, parameterLength)) || (PARAMETER_CHECK("-a", 2, parameterLength))) {
            if ((i+1) < argc) {
                alphabet = argv[i + 1];
                i++;
            }
        } else if ((PARAMETER_CHECK("--complement", 12, parameterLength)) || (PARAMETER_CHECK("-c", 2, parameterLength))) {
            if ((i+1) < argc) {
                complement = argv[i + 1];
                i++;
            }
        } else if ((PARAMETER_CHECK("--threads", 9, parameterLength)) || (PARAMETER_CHECK("-t", 2, parameterLength))) {
            if ((i+1) < argc) {
                threads = stoul(argv[i + 1]);
                i++;
            }
        } else if (PARAMETER_CHECK("--quiet", 7, parameterLength)) {
            verbose_flag = false;
        } else if (PARAMETER_CHECK("--no-rcomp", 10, parameterLength)) {
            rcomp_flag = false;
        } else if (PARAMETER_CHECK("--gzip", 6, parameterLength)) {
            gzip_flag = true;
        } else {
            getcorpus_help(1);
        }
    }

    if (haveK && haveFastaIn && haveBed) {
        Bed2Corpus *bc = new Bed2Corpus(k, fastaInFile, bedFile, stride, alphabet, complement,
                                        sep, threads, verbose_flag, rcomp_flag, gzip_flag);
        if (haveCorpusOut) {
            if (haveIndexOut) {
                bc->NewIndexNewCorpus(corpusOutFile, indexOutFile);
            } else if (haveIndexIn) {
                bc->ExistingIndexNewCorpus(indexInFile, corpusOutFile);
            } else {
                bc->NewCorpusOnly(corpusOutFile);
            }
        } else if (haveIndexOut) {
            bc->NewIndexOnly(indexOutFile);
        }
    } else if (haveCorpusIn && haveIndexIn && haveCorpusOut) {
        Bed2Corpus *bc = new Bed2Corpus(corpusInFile, indexInFile, corpusOutFile, sep,
                                        threads, verbose_flag, rcomp_flag, gzip_flag);
        bc->ExistingIndexExistingCorpus();
    } else {
        getcorpus_help(1);
    }
    return 0;
}


void getcorpus_help(int code) {
    cerr << "\nTool:    corpustool getcorpus " << endl;
    cerr << "Summary: Extract kmer tokens from a FASTA file based on feature coordinates." << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME 
         << " [OPTIONS] -k <kmer_size> -fi <fasta> -bed <bed/gff/vcf> -co <text_file>" << endl
         << "    or   " << PROGRAM_NAME
         << " [OPTIONS] -k <kmer_size> -fi <fasta> -bed <bed/gff/vcf> -co <text_file> -io <index_file>" << endl
         << "    or   " << PROGRAM_NAME
         << " [OPTIONS] -k <kmer_size> -fi <fasta> -bed <bed/gff/vcf> -ii <index_file> -co <text_file>" << endl
         << "    or   " << PROGRAM_NAME
         << " [OPTIONS] -k <kmer_size> -fi <fasta> -bed <bed/gff/vcf> -io <index_file>" << endl
         << "    or   " << PROGRAM_NAME
         << " [OPTIONS] -ci <corpus> -ii <index> -co <indexed_corpus>" << endl
         << endl;

    cerr << "Options:   " << endl;
    cerr << "\t-k                 Size of kmer tokens." << endl;
    cerr << "\t-fi                Input FASTA file." << endl;
    cerr << "\t-bed               BED/GFF/VCF file of ranges to extract from -fi." << endl;
    cerr << "\t-ii                Input index file." << endl;
    cerr << "\t-ci                Input corpus file (plain text)." << endl;
    cerr << "\t-io                Output index file." << endl;
    cerr << "\t-co                Output (indexed) corpus file (plain text)." << endl << endl;
    
    cerr << "\t-s, --stride       Stride for breaking sequence into kmers." << endl;
    cerr << "\t-a, --alphabet     Sequence alphabet." << endl;
    cerr << "\t-c, --complement   Sequence complement alphabet." << endl;
    cerr << "\t-t, --threads      Number of threads to run on." << endl;
    cerr << "\t--no-rcomp         Do not treat kmers and their reverse complements the same" << endl;
    cerr << "\t                   (otherwise, kmers are represented as max(kmer, rev_comp))." << endl;
    cerr << "\t--gzip             Compress output corpus file in gzip format." << endl;
    cerr << "\t--quiet            Quiet mode." << endl << endl;
    

    exit(code);
}