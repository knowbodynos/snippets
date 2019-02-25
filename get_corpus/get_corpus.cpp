#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <map>
#include <math.h>
#include <algorithm>
#include <thread>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <libgen.h>
#include <linux/limits.h>
#include <sys/stat.h>
#include <unistd.h>

using namespace std;

static int build_index_flag = 0;
static int index_only_flag = 0;
static int header_flag = 1;
static int verbose_flag = 1;
static int rcomp_flag = 1;
static unsigned int MAX_THREADS = thread::hardware_concurrency();


// Simply gets the file size of file.
int getFileSize(ifstream *file) {
    file->seekg(0,ios::end);
    int filesize = file->tellg();
    file->seekg(ios::beg);
    return filesize;
}


// Chunks a file by breaking it up into chunks of "chunkSize" bytes.
void chunkFile(string fullFilePath, string chunkName, unsigned long chunkSize) {
    ifstream fileStream;
    fileStream.open(fullFilePath.c_str(), ios::in | ios::binary);

    // File open a success
    if (fileStream.is_open()) {
        ofstream output;
        int counter = 1;

        string fullChunkName;
        string line;

        if (header_flag) {
            string header;
            getline(fileStream, header);
        }

        getline(fileStream, line);

        // Keep reading until end of file
        while (!fileStream.eof()) {

            // Build the chunk file name. Usually drive:\\chunkName.ext.N
            // N represents the Nth chunk
            fullChunkName.clear();
            fullChunkName.append(chunkName);
            fullChunkName.append(".");

            // Convert counter integer into string and append to name.
            fullChunkName.append(to_string(counter));

            // Open new chunk file name for output
            output.open(fullChunkName.c_str(), ios::out | ios::trunc | ios::binary);

            // If chunk file opened successfully, read from input and 
            // write to output chunk. Then close.
            if (output.is_open()) {
                for (int i=0;(i<chunkSize)&&(!fileStream.eof());i++) {
                    // gcount() returns number of bytes read from stream.
                    output << line << endl;
                    getline(fileStream, line);
                }
                output.close();
                counter++;
            } else {
                cerr << "Error: Unable to open file '" << fullChunkName << "' for writing." << endl;
            }
        }

        // Close input file stream.
        fileStream.close();
        if (verbose_flag) {
            cout << "-> Chunking complete! " << counter - 1 << " files created." << endl;
        }
    } else {
        cerr << "Error: Unable to open file '" << fullFilePath << "' for reading." << endl;
    }
}


// Finds chunks by "chunkName" and creates file specified in outputFileName
void joinFile(string chunkName, string outputFileName, string header, int column, char delim, string ext="") {
    string inputFileName;

    // Create our output file
    ofstream fileOutput;
    fileOutput.open(outputFileName.c_str(), ios::out | ios::binary);

    // If successful, loop through chunks matching chunkName
    if (fileOutput.is_open()) {
        bool filefound = true;
        int counter = 1;
        int fileSize = 0;

        if (header_flag) {
            istringstream iss(header);
            string field;
            for (int i=0;i<column;i++) {
                getline(iss, field, delim);
                fileOutput << field << delim;
            }
            getline(iss, field, delim);
            fileOutput << "Tokens";
            while (!iss.eof()) {
                getline(iss, field, delim);
                fileOutput << delim << field;
            }
            fileOutput << endl;
        }

        while (filefound) {
            filefound = false;

            // Build the filename
            inputFileName.clear();
            inputFileName.append(chunkName);
            inputFileName.append(".");
            inputFileName.append(to_string(counter));
            if (ext != "") {
                inputFileName.append(".");
                inputFileName.append(ext);
            }
            // Open chunk to read
            ifstream fileInput;
            fileInput.open(inputFileName.c_str(), ios::in | ios::binary);

            // If chunk opened successfully, read it and write it to 
            // output file.
            if (fileInput.is_open()) {
                filefound = true;
                fileSize = getFileSize(&fileInput);
                char *inputBuffer = new char[fileSize];

                fileInput.read(inputBuffer, fileSize);
                fileOutput.write(inputBuffer, fileSize);
                delete(inputBuffer);

                fileInput.close();
            }
            counter++;
        }

        // Close output file.
        fileOutput.close();

        if (verbose_flag) {
            char resolved_path[PATH_MAX];
            cout << "-> Joining complete! Corpus saved to '" << realpath(&outputFileName[0], resolved_path) << "'." << endl;
        }
    } else {
        cerr << "Error: Unable to open file '" << outputFileName << "' for writing." << endl;
    }
}


string translate(string token, string alphabet, string complement) {
    string origToken = token;
    for (int i=0;i<alphabet.length();i++) {
        size_t pos = origToken.find(alphabet[i], 0);
        while(pos != string::npos)
        {
            token[pos] = complement[i];
            pos = origToken.find(alphabet[i], pos + 1);
        }
    }
    reverse(token.begin(),token.end());
    if (token > origToken) {
        return token;
    } else {
        return origToken;
    }
}


string tokenizeSeq(string seq, map<string,string> &rev_comp, int k, int stride, string alphabet, string complement, char sep) {
    transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    string origSeq = seq;
    string kmer = seq.substr(0, k);
    string tokens = translate(kmer, alphabet, complement);
    for (int i=stride;i<seq.length();i+=stride) {
        string kmer = seq.substr(i, k);
        if (kmer.length() == k) {
            string token;
            if (rcomp_flag) {
                if (rev_comp.find(kmer) == rev_comp.end()) {
                    token = translate(kmer, alphabet, complement);
                    rev_comp[kmer] = token;
                } else {
                    token = rev_comp[kmer];
                }
            } else {
                token = kmer;
            }
            tokens += sep;
            tokens.append(token);
        }
    }
    return tokens;
}


// Tokenizes a sequence tsv file.
void tokenizeFile(string fullFilePath, int k, int stride, string alphabet, string complement, int column, char delim, char sep) {
    map<string,string> rev_comp;
    ifstream fileStream;
    fileStream.open(fullFilePath.c_str(), ios::in | ios::binary);

    // File open a success
    if (fileStream.is_open()) {
        // Open new file name for output
        string outFilePath = fullFilePath;
        outFilePath.append(".corpus");

        ofstream output;
        output.open(outFilePath.c_str(), ios::out | ios::trunc | ios::binary);
        if (output.is_open()) {
            string line;
            string field;
            // Keep reading until end of file
            while (getline(fileStream, line)) {
                istringstream iss(line);
                for (int i=0;i<column;i++) {
                    getline(iss, field, delim);
                    output << field << delim;
                }
                getline(iss, field, delim);
                output << tokenizeSeq(field, rev_comp, k, stride, alphabet, complement, sep) << endl;
            }
            output.close();
        } else {
            cerr << "Error: Unable to open file '" << outFilePath << "' for writing." << endl;
        }

        // Close input file stream.
        fileStream.close();
    } else {
        cerr << "Error: Unable to open file '" << fullFilePath << "' for reading." << endl;
    }
}


class fileTokenizer {
    private:
        string alphabet, complement;
        char delim, sep;
        int k, stride, column;
    public:
        fileTokenizer(int k1, int stride1, string alphabet1, string complement1, int column1, char delim1, char sep1) {
            k = k1;
            stride = stride1;
            alphabet = alphabet1;
            complement = complement1;
            column = column1;
            delim = delim1;
            sep = sep1;
        }
        void operator()(string fullFilePath) { 
            tokenizeFile(fullFilePath, k, stride, alphabet, complement, column, delim, sep);
        }
};


// Populate token frequencies from an input corpus.
void popFreqs(map<string,int> &freq, string fullFilePath, int column, char delim, char sep) {
    ifstream fileStream;
    fileStream.open(fullFilePath.c_str(), ios::in | ios::binary);

    // File open a success
    if (fileStream.is_open()) {
        string line;
        string field;
        // Keep reading until end of file
        while (getline(fileStream, line)) {
            istringstream iss1(line);
            for (int i=0;i<column;i++) {
                getline(iss1, field, delim);
            }
            getline(iss1, field, delim);
            string token;
            istringstream iss2(field);
            while (getline(iss2, token, sep)) {
                if (freq.find(token) == freq.end()) {
                    freq[token] = 1;
                } else {
                    freq[token]++;
                }
            }
        }
        fileStream.close();
    } else {
        cerr << "Error: Unable to open file '" << fullFilePath << "' for reading." << endl;
    }
}


class freqPopulator {
    private:
        int column;
        char delim, sep;
    public:
        freqPopulator(int column1, char delim1, char sep1) {
            column = column1;
            delim = delim1;
            sep = sep1;
        }
        void operator()(map<string,int> &freq, string fullFilePath) { 
            popFreqs(freq, fullFilePath, column, delim, sep);
        }
};


// Indexes the tokens in a corpus tsv file.
map<string,int> buildIndex(string fullIndexPath, map<string, int> freq, char delim) {
    vector<pair<string, int>> pairs;
    map<string,int> index;
    for (auto itr = freq.begin(); itr != freq.end(); ++itr)
        pairs.push_back(*itr);

    sort(pairs.begin(), pairs.end(), [](pair<string, int> a, pair<string, int> b) {
        return a.second > b.second;
    });

    // Open new file name for output
    ofstream output;
    output.open(fullIndexPath.c_str(), ios::out | ios::trunc | ios::binary);
    if (output.is_open()) {
        output << "Token" << delim << "Counts" << endl; 
        int i = 1;
        for (auto itr = pairs.begin(); itr != pairs.end(); ++itr) {
            index[itr->first] = i;
            output << itr->first << "\t" << itr->second << endl;
            i++;
        }
        output.close();
    } else {
        cerr << "Error: Unable to open file '" << fullIndexPath << "' for writing." << endl;
    }
    if (verbose_flag) {

        char resolved_path[PATH_MAX];
        cout << "-> Building index complete! Indexes saved to '" << realpath(&fullIndexPath[0], resolved_path) << "'." << endl;
    }
    return index;
}


string indexTokens(string tokens, map<string,int> index, char sep) {
    string token;
    istringstream iss(tokens);
    getline(iss, token, sep);
    string indexedTokens = to_string(index[token]);
    while(getline(iss, token, sep)) {
        indexedTokens += sep;
        indexedTokens.append(to_string(index[token]));
    }
    return indexedTokens;
}


void indexFile(string fullFilePath, map<string,int> index, int column, char delim, char sep) {
    ifstream fileStream;
    fileStream.open(fullFilePath.c_str(), ios::in | ios::binary);
    // File open a success
    if (fileStream.is_open()) {
        // Open new file name for output
        string outFilePath = fullFilePath;
        outFilePath.append(".index");
        ofstream output;
        output.open(outFilePath.c_str(), ios::out | ios::trunc | ios::binary);
        if (output.is_open()) {
            string line;
            string field;
            // Keep reading until end of file
            while (getline(fileStream, line)) {
                istringstream iss(line);
                for (int i=0;i<column;i++) {
                    getline(iss, field, delim);
                    output << field << delim;
                }
                getline(iss, field, delim);
                output << indexTokens(field, index, sep) << endl;
            }
            output.close();
        } else {
            cerr << "Error: Unable to open file '" << outFilePath << "' for writing." << endl;
        }
        fileStream.close();
    } else {
        cerr << "Error: Unable to open file '" << fullFilePath << "' for reading." << endl;
    }
}


class fileIndexer {
    private:
        map<string, int> index;
        int column;
        char delim, sep;
    public:
        fileIndexer(map<string, int> index1, int column1, char delim1, char sep1) {
            index = index1;
            column = column1;
            delim = delim1;
            sep = sep1;
        }
        void operator()(string fullFilePath) {
            indexFile(fullFilePath, index, column, delim, sep);
        }
};


map<string, int> readIndex(string indexFile, char delim) {
    map<string, int> index;
    ifstream fileStream;
    fileStream.open(indexFile.c_str(), ios::in | ios::binary);
    // File open a success
    if (fileStream.is_open()) {
        string line;
        getline(fileStream, line);
        int i = 1;
        while (getline(fileStream, line)) {
            istringstream iss(line);
            string token;
            string counts;
            getline(iss, token, delim);
            getline(iss, counts, delim);
            index[token] = i;
            i++;
        }
        fileStream.close();
    }
    return index;
}


template <class T>
void doThreading(T func, string chunkName, int threads, string msg, string ext="") {
    string fullChunkName;
    vector<thread> threadList;
    for (int i=0;i<threads;i++) {
        fullChunkName.clear();
        fullChunkName.append(chunkName);
        fullChunkName.append(".");
        fullChunkName.append(to_string(i + 1));
        if (ext != "") {
            fullChunkName.append(".");
            fullChunkName.append(ext);
        }
        // cout << fullChunkName << endl;
        threadList.push_back(thread(func, fullChunkName));
        // threadList.push_back(thread(tokenizeFile, fullChunkName, 5, 1, "ACGTNBDHU", "TGCANUHDB", " "));
    }
    for (auto& th : threadList) {
        th.join();
    }
    if (verbose_flag) {
        cout << "-> " << msg << endl;
    }
}


template <class T>
void doThreading(T func, vector<map<string,int>> &freqs, string chunkName, int threads, string msg, string ext="") {
    string fullChunkName;
    vector<thread> threadList;
    for (int i=0;i<threads;i++) {
        fullChunkName.clear();
        fullChunkName.append(chunkName);
        fullChunkName.append(".");
        fullChunkName.append(to_string(i + 1));
        if (ext != "") {
            fullChunkName.append(".");
            fullChunkName.append(ext);
        }
        // cout << fullChunkName << endl;
        threadList.push_back(thread(func, ref(freqs[i]), fullChunkName));
        // threadList.push_back(thread(tokenizeFile, fullChunkName, 5, 1, "ACGTNBDHU", "TGCANUHDB", " "));
    }
    for (auto& th : threadList) {
        th.join();
    }
    if (verbose_flag) {
        cout << "-> " << msg << endl;
    }
}


void removeTmpfiles(string chunkName, int threads, string ext="") {
    string fullChunkName;
    for (int i=0;i<threads;i++) {
        fullChunkName.clear();
        fullChunkName.append(chunkName);
        fullChunkName.append(".");
        fullChunkName.append(to_string(i + 1));
        if (ext != "") {
            fullChunkName.append(".");
            fullChunkName.append(ext);
        }
        if (remove(fullChunkName.c_str()) != 0) {
            cerr << "Error: Could not remove intermediate file '" << fullChunkName << "'.";
        }
    }
}


void usage(int status) {
    const char *text = "usage: get_corpus [-acdhostz] [--no-header] [--quiet] [--no-rcomp] ... [Sequence Corpus File] [kmer size] [stride length]\n"
    "\n"
    "       Tokenize sequence corpus.\n"
    "\n"
    "       Options:\n"
    "           --alphabet, -a      Sequence alphabet.\n"
    "           --col, -c           Complement alphabet.\n"
    "           --delim, -d         Input table delimitor.\n"
    "           --help, -h          Display this help and exit.\n"
    "           --index, -i         Index file path.\n"
    "           --out, -o           Output file path.\n"
    "           --sep, -s           Kmer token seperator.\n"
    "           --threads, -t       Number of threads to run on.\n"
    "           --complement, -z    Complement alphabet.\n"
    "\n"
    "           --build-index       Build new kmer token index from input sequences.\n"
    "           --no-header         Input file has no header.\n"
    "           --quiet             Quiet mode.\n"
    "           --no-rcomp          Do not treat kmers and their reverse complements the same\n"
    "                               (otherwise, kmers are represented as max(kmer, rev_comp)).\n\n";
    
    
    cout << text;
    exit(status);
}


int main(int argc, char **argv) {
    string alphabet = "ACGTNBDHVRYKM";
    string complement = "TGCANVHDBYRMK";
    char sep = ' ';
    int column = 6;
    char delim = '\t';
    int threads = MAX_THREADS;
    string outFile = "corpus.tsv";
    string indexFile = "";
    int opt;
    while (1) {
        static struct option long_options[] = {
            /* These options set a flag. */
            {"build-index", no_argument, &build_index_flag, 1},
            {"index-only", no_argument, &index_only_flag, 1},
            {"no-header", no_argument, &header_flag, 0},
            {"quiet", no_argument, &verbose_flag, 0},
            // Match rcomp as well
            {"no-rcomp", no_argument, &rcomp_flag, 0},
            /* These options don’t set a flag.
            We distinguish them by their indices. */
            {"alphabet", required_argument, 0, 'a'},
            
            {"col", required_argument, 0, 'c'},
            {"delim", required_argument, 0, 'd'},
            {"help", no_argument, 0, 'h'},
            {"index", no_argument, 0, 'i'},
            {"out", required_argument, 0, 'o'},
            {"sep", required_argument, 0, 's'},
            {"threads", required_argument, 0, 't'},
            {"complement", required_argument, 0, 'z'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        opt = getopt_long(argc, argv, "a:c:d:hi:o:s:t:z:",
                          long_options, &option_index);

        /* Detect the end of the options. */
        if (opt == -1) {
            break;
        }

        switch (opt) {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0) {
                    break;
                }
                printf ("option %s", long_options[option_index].name);
                if (optarg) {
                    printf (" with arg %s", optarg);
                }
                printf ("\n");
                break;

            case 'a':
                if (optarg) {
                    alphabet = optarg;
                } else {
                    usage(1);
                }
                break;

            case 'c':
                if (optarg) {
                    column = stoi(optarg);
                } else {
                    usage(1);
                }
                break;

            case 'd':
                if (optarg && optarg[0]) {
                    delim = optarg[0];
                } else {
                    usage(1);
                }
                break;

            case 'h':
                usage(0);
                break;

            case 'i':
                if (optarg) {
                    indexFile = optarg;
                } else {
                    usage(1);
                }
                break;

            case 'o':
                if (optarg) {
                    outFile = optarg;
                } else {
                    usage(1);
                }
                break;

            case 's':
                if (optarg && optarg[0]) {
                    sep = optarg[0];
                } else {
                    usage(1);
                }
                break;

            case 't':
                if (optarg) {
                    int new_threads = stoi(optarg);
                    if (new_threads <= MAX_THREADS) {
                        threads = new_threads;
                    }
                } else {
                    usage(1);
                }
                break;

            case 'z':
                if (optarg) {
                    complement = optarg;
                } else {
                    usage(1);
                }
                break;

            case '?':
                /* getopt_long already printed an error message. */
                break;

            default:
                usage(1);
        }
    }

    string fullFilePath = argv[optind];
    int k = stoi(argv[optind + 1]);
    int stride = stoi(argv[optind + 2]);
    // string chunkName = fullFilePath.substr(fullFilePath.find_last_of("/\\") + 1);
    string tmp_dir = "/tmp/get_corpus/";
    string chunkName = tmp_dir;
    chunkName.append(basename(&fullFilePath[0]));
    string fullChunkName;
    string fullTokenizedChunkName;
    string header;
    map<string, int> index;
    ifstream fileStream;
    fileStream.open(fullFilePath.c_str(), ios::in | ios::binary);

    if (build_index_flag && (indexFile == "")) {
        indexFile = dirname(&outFile[0]);
        indexFile.append("/index.tsv");
    }

    int lines = 0;
    if (fileStream.is_open()) {
        string line;
        if (header_flag) {
            getline(fileStream, header);
        }
        while (getline(fileStream, line)) {
            lines++;
        }
        fileStream.close();
    }

    if (verbose_flag) {
        char resolved_path[PATH_MAX];
        cout << "The file '" << realpath(&fullFilePath[0], resolved_path) << "' has " << lines << " lines." << endl;
        cout << "There are " << threads << " out of " << MAX_THREADS << " threads in use." << endl << endl;
    }

    int chunkSize = ceil((float)lines / threads);

    if (mkdir(tmp_dir.c_str(), 0777) && (errno != EEXIST)) {
        cerr << "Error: Unable to create directory '" << tmp_dir << "'.";
    }

    chunkFile(fullFilePath, chunkName, chunkSize);

    fileTokenizer tokenizer(k, stride, alphabet, complement, column, delim, sep);
    doThreading(tokenizer, chunkName, threads, "Tokenizing complete!");

    removeTmpfiles(chunkName, threads);
    
    if (build_index_flag) {
        freqPopulator populator(column, delim, sep);
        vector<map<string,int>> freqs(threads);
        doThreading(populator, freqs, chunkName, threads, "Counting complete!", "corpus");
        map<string,int> freq;
        for (map<string,int> f : freqs) {
            for (auto itr = f.begin(); itr != f.end(); ++itr) {
                if (freq.find(itr->first) == freq.end()) {
                    freq[itr->first] = itr->second;
                } else {
                    freq[itr->first] += itr->second;
                }
            }
        }
        index = buildIndex(indexFile, freq, delim);
    }

    if (!index_only_flag) {
        string joinExt = "corpus";

        if (indexFile != "") {
            if (!build_index_flag) {
                index = readIndex(indexFile, delim);
            }
            fileIndexer indexer(index, column, delim, sep);
            doThreading(indexer, chunkName, threads, "Indexing complete!", "corpus");
            joinExt.append(".index");
        }

        joinFile(chunkName, outFile, header, column, delim, joinExt);
    }

    removeTmpfiles(chunkName, threads, "corpus");
    if (indexFile != "") {
        removeTmpfiles(chunkName, threads, "corpus.index");
    }

    if (rmdir(tmp_dir.c_str())) {
        cerr << "Error: Unable to remove directory '" << tmp_dir << "'.";
    }
}