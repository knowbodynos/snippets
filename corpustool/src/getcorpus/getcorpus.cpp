/*****************************************************************************
  getcorpus.cpp

  2019 - Ross Altman
  Inari Agriculture
  raltman@inari.com
******************************************************************************/
#include "getcorpus.h"

using namespace std;


// High level functions


void Bed2Corpus::NewCorpusOnly(string corpusOutFile) {
    doChunkTokenize();
    joinFile(_chunkName + ".corpus", corpusOutFile, _threads, _gzip_flag);    
    if (_verbose_flag) {
        char resolved_path[PATH_MAX];
        cout << "-> Joining complete! Corpus saved to '" << realpath(&corpusOutFile[0], resolved_path) << "'." << endl;
    }
    removeTmpfiles(_chunkName + ".corpus", _threads); 
    if (rmdir(tmp_dir.c_str())) {
        cerr << "Error: Unable to remove directory '" << tmp_dir << "'.";
    }
}


void Bed2Corpus::NewIndexNewCorpus(string corpusOutFile, string indexOutFile) {
    doChunkTokenize();
    doNewIndex(indexOutFile);
    doIndexCorpus(corpusOutFile);
    removeTmpfiles(_chunkName + ".corpus", _threads);
    if (rmdir(tmp_dir.c_str())) {
        cerr << "Error: Unable to remove directory '" << tmp_dir << "'.";
    }
}


void Bed2Corpus::ExistingIndexNewCorpus(string indexInFile, string corpusOutFile) {
    doChunkTokenize();
    _index = readIndex(indexInFile);
    doIndexCorpus(corpusOutFile);
    removeTmpfiles(_chunkName + ".corpus", _threads);
    if (rmdir(tmp_dir.c_str())) {
        cerr << "Error: Unable to remove directory '" << tmp_dir << "'.";
    }
}


void Bed2Corpus::NewIndexOnly(string indexOutFile) {
    doChunkTokenize();
    doNewIndex(indexOutFile);
    removeTmpfiles(_chunkName + ".corpus", _threads);
    if (rmdir(tmp_dir.c_str())) {
        cerr << "Error: Unable to remove directory '" << tmp_dir << "'.";
    }
}


void Bed2Corpus::ExistingIndexExistingCorpus(void) {
    _chunkName = tmp_dir + basename(&_corpusInFile[0]);
    _chunkName = _chunkName.substr(0, _chunkName.find_last_of("."));
    int lines = 0;
    int ext_start = _corpusInFile.find_last_of(".");
    string ext = _corpusInFile.substr(ext_start + 1, _corpusInFile.length() - ext_start);
    if (ext == ".gz") {
        igzstream inFileStream;
        inFileStream.open(_corpusInFile.c_str(), ios::in | ios::binary);
        if (!inFileStream.fail()) {
            string line;
            while (getline(inFileStream, line)) {
                lines++;
            }
            inFileStream.close();
        }
    } else {
        ifstream inFileStream;
        inFileStream.open(_corpusInFile.c_str(), ios::in | ios::binary);
        if (!inFileStream.fail()) {
            string line;
            while (getline(inFileStream, line)) {
                lines++;
            }
            inFileStream.close();
        }
    }

    int chunkSize = (int)ceil((double)lines / _threads);
    int files = chunkFile(_corpusInFile, _chunkName, _threads, chunkSize);
    if (_verbose_flag) {
        cout << "-> Chunking complete! " << files << " files created." << endl;
    }

    _index = readIndex(_indexInFile);

    doIndexCorpus(_corpusOutFile);
    removeTmpfiles(_chunkName + ".corpus", _threads);
    if (rmdir(tmp_dir.c_str())) {
        cerr << "Error: Unable to remove directory '" << tmp_dir << "'.";
    }
}


// Intermediate methods


void Bed2Corpus::doChunkTokenize(void) {
    _chunkName = tmp_dir + basename(&_bedFile[0]);
    _chunkName = _chunkName.substr(0, _chunkName.find_last_of("."));
    // string faFilePath = "/home/ec2-user/data/ref/Zea_mays.B73_RefGen_v4.fa";
    // string gtfFilePath = "/home/ec2-user/data/ref/five_prime_utr_u1000d0.gtf";
    BedFileCorpus *bf = new BedFileCorpus(_bedFile);
    bf->Open();

    FastaReference *fr = new FastaReference;
    bool memmap = true;
    fr->open(_fastaInFile, memmap);

    if (mkdir(tmp_dir.c_str(), 0777) && (errno != EEXIST)) {
        cerr << "Error: Unable to create directory '" << tmp_dir << "'." << endl;
    }

    int chunkSize = bf->GetChunkSize(fr, _threads);
    cout << "Threads: " << _threads << endl;
    cout << "Chunk Size: " << chunkSize << endl;
    int files = bf->ChunkFile(fr, _chunkName, _threads, chunkSize);
    if (_verbose_flag) {
        cout << "-> Chunking complete! " << files << " files created." << endl;
    }

    fileTokenizer tokenizer(_k, _stride, _alphabet, _complement, _sep, _rcomp_flag);
    broadcast(tokenizer, _chunkName, _threads);
    if (_verbose_flag) {
        cout << "-> Tokenizing complete!" << endl;
    }

    removeTmpfiles(_chunkName, _threads);
}


void Bed2Corpus::doNewIndex(string indexOutFile) {
    _freqs.resize(_threads);
    freqPopulator populator(_sep);
    broadcast(populator, _chunkName, _freqs, _threads);
    if (_verbose_flag) {
        cout << "-> Counting complete!" << endl;
    }
    map<string,int> freq;
    for (map<string,int> f : _freqs) {
        for (auto itr = f.begin(); itr != f.end(); ++itr) {
            if (freq.find(itr->first) == freq.end()) {
                freq[itr->first] = itr->second;
            } else {
                freq[itr->first] += itr->second;
            }
        }
    }
    _index = buildIndex(indexOutFile, freq);
    if (_verbose_flag) {
        char resolved_path[PATH_MAX];
        cout << "-> Building index complete! Indexes saved to '" << realpath(&indexOutFile[0], resolved_path) << "'." << endl;
    }
}


void Bed2Corpus::doIndexCorpus(string corpusOutFile) {
    fileIndexer indexer(_index, _sep);
    broadcast(indexer, _chunkName, _threads);
    if (_verbose_flag) {
        cout << "-> Indexing complete!" << endl;
    }
    joinFile(_chunkName + ".icorpus", corpusOutFile, _threads, _gzip_flag);
    if (_verbose_flag) {
        char resolved_path[PATH_MAX];
        cout << "-> Joining complete! Corpus saved to '" << realpath(&corpusOutFile[0], resolved_path) << "'." << endl;
    }
    removeTmpfiles(_chunkName + ".icorpus", _threads);
}


// Helper functions

string padZeros(int val, int threads) {
    int zeros = to_string(threads).length();
    string val_padded = to_string(val);
    return string(zeros - val_padded.length(), '0') + val_padded;
}

// Simply gets the file size of file.
int getFileSize(ifstream *file) {
    file->seekg(0,ios::end);
    int filesize = file->tellg();
    file->seekg(ios::beg);
    return filesize;
}


// Chunks a file by breaking it up into chunks of "chunkSize" bytes.
int chunkFile(string inFileName, string chunkName, int threads, unsigned long chunkSize) {
    ifstream inFileStream;
    inFileStream.open(inFileName.c_str(), ios::in | ios::binary);
    // File open a success
    if (inFileStream.is_open()) {
        string outFileName;
        ofstream outFileStream;
        string line;
        getline(inFileStream, line);
        int counter = 1;
        // Keep reading until end of file
        while (!inFileStream.eof()) {
            // Build the chunk file name. Usually drive:\\chunkName.ext.N
            // N represents the Nth chunk
            outFileName.clear();
            outFileName.append(chunkName);
            outFileName.append(".");
            outFileName.append("corpus");
            outFileName.append(".");
            outFileName.append(padZeros(counter, threads));

            // Open new chunk file name for output
            outFileStream.open(outFileName.c_str(), ios::out | ios::trunc | ios::binary);

            // If chunk file opened successfully, read from input and 
            // write to output chunk. Then close.
            if (outFileStream.is_open()) {
                for (size_t i=0;(i<chunkSize)&&(!inFileStream.eof());i++) {
                    // gcount() returns number of bytes read from stream.
                    outFileStream << line << endl;
                    getline(inFileStream, line);
                }
                outFileStream.close();
                counter++;
            } else {
                cerr << "Error: Unable to open file '" << outFileName << "' for writing." << endl;
            }
        }
        // Close input file stream.
        inFileStream.close();

        return counter - 1;
    } else {
        cerr << "Error: Unable to open file '" << inFileName << "' for reading (chunkFile)." << endl;
        return 0;
    }
}


// Finds chunks by "chunkName" and creates file specified in outFileName
void joinFile(string chunkName, string outFileName, int threads, bool gzip_flag) {
    string inFileName;

    // Create our output file
    if (gzip_flag) {
        ogzstream outFileStream;
        outFileName.append(".gz");
        outFileStream.open(outFileName.c_str(), ios::out | ios::binary);

        // If successful, loop through chunks matching chunkName
        if (!outFileStream.fail()) {
            bool fileFound = true;
            int counter = 1;
            int fileSize = 0;

            while (fileFound) {
                fileFound = false;

                // Build the filename
                inFileName.clear();
                inFileName.append(chunkName);
                inFileName.append(".");
                inFileName.append(padZeros(counter, threads));

                // Open chunk to read
                ifstream inFileStream;
                inFileStream.open(inFileName.c_str(), ios::in | ios::binary);

                // If chunk opened successfully, read it and write it to 
                // output file.
                if (inFileStream.is_open()) {
                    fileFound = true;
                    fileSize = getFileSize(&inFileStream);
                    char *inputBuffer = new char[fileSize];

                    inFileStream.read(inputBuffer, fileSize);
                    outFileStream.write(inputBuffer, fileSize);
                    delete(inputBuffer);

                    inFileStream.close();
                }
                counter++;
            }

            // Close output file.
            outFileStream.close();
        } else {
            cerr << "Error: Unable to open file '" << outFileName << "' for writing." << endl;
        }
    } else {
        ofstream outFileStream;
        outFileStream.open(outFileName.c_str(), ios::out | ios::binary);

        // If successful, loop through chunks matching chunkName
        if (!outFileStream.fail()) {
            bool fileFound = true;
            int counter = 1;
            int fileSize = 0;

            while (fileFound) {
                fileFound = false;

                // Build the filename
                inFileName.clear();
                inFileName.append(chunkName);
                inFileName.append(".");
                inFileName.append(padZeros(counter, threads));

                // Open chunk to read
                ifstream inFileStream;
                inFileStream.open(inFileName.c_str(), ios::in | ios::binary);

                // If chunk opened successfully, read it and write it to 
                // output file.
                if (inFileStream.is_open()) {
                    fileFound = true;
                    fileSize = getFileSize(&inFileStream);
                    char *inputBuffer = new char[fileSize];

                    inFileStream.read(inputBuffer, fileSize);
                    outFileStream.write(inputBuffer, fileSize);
                    delete(inputBuffer);

                    inFileStream.close();
                }
                counter++;
            }

            // Close output file.
            outFileStream.close();
        } else {
            cerr << "Error: Unable to open file '" << outFileName << "' for writing." << endl;
        }
    }
}


string translate(string token, string alphabet, string complement) {
    string origToken = token;
    for (size_t i=0;i<alphabet.length();i++) {
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


string tokenizeSeq(string seq, map<string,string> &rev_comp, int k, int stride,
                   string alphabet, string complement, char sep, bool rcomp_flag) {
    transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    string origSeq = seq;
    string kmer = seq.substr(0, k);
    string tokens = "";
    if ((int)kmer.length() == k) {
        tokens = translate(kmer, alphabet, complement);
    }
    for (size_t i=stride;i<seq.length();i+=stride) {
        string kmer = seq.substr(i, k);
        if ((int)kmer.length() == k) {
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
void tokenizeFile(string inFileName, int k, int stride, string alphabet, string complement,
                  char sep, bool rcomp_flag) {
    map<string,string> rev_comp;
    ifstream inFileStream;
    inFileStream.open(inFileName.c_str(), ios::in | ios::binary);

    // File open a success
    if (inFileStream.is_open()) {
        // Open new file name for output
        string outFileName = inFileName.insert(inFileName.find_last_of("."), ".corpus");

        ofstream outFileStream;
        outFileStream.open(outFileName.c_str(), ios::out | ios::trunc | ios::binary);
        if (outFileStream.is_open()) {
            string line;
            // Keep reading until end of file
            while (getline(inFileStream, line)) {
                outFileStream << tokenizeSeq(line, rev_comp, k, stride, alphabet, complement, sep, rcomp_flag) << endl;
            }
            outFileStream.close();
        } else {
            cerr << "Error: Unable to open file '" << outFileName << "' for writing." << endl;
        }

        // Close input file stream.
        inFileStream.close();
    } else {
        cerr << "Error: Unable to open file '" << inFileName << "' for reading (tokenizeFile)." << endl;
    }
}


// Populate token frequencies from an input corpus.
void popFreqs(string inFileName, map<string,int> &freq, char sep) {
    inFileName.insert(inFileName.find_last_of("."), ".corpus");
    ifstream inFileStream;
    inFileStream.open(inFileName.c_str(), ios::in | ios::binary);

    // File open a success
    if (inFileStream.is_open()) {
        string line;
        // Keep reading until end of file
        while (getline(inFileStream, line)) {
            string token;
            istringstream iss(line);
            while (getline(iss, token, sep)) {
                if (freq.find(token) == freq.end()) {
                    freq[token] = 1;
                } else {
                    freq[token]++;
                }
            }
        }
        inFileStream.close();
    } else {
        cerr << "Error: Unable to open file '" << inFileName << "' for reading (popFreqs)." << endl;
    }
}


// Indexes the tokens in a corpus tsv file.
map<string,int> buildIndex(string outFileName, map<string,int> freq) {
    vector<pair<string, int>> pairs;
    map<string,int> index;
    for (auto itr = freq.begin(); itr != freq.end(); ++itr)
        pairs.push_back(*itr);

    sort(pairs.begin(), pairs.end(), [](pair<string, int> a, pair<string, int> b) {
        return a.second > b.second;
    });
    // Open new file name for output
    ofstream outFileStream;
    outFileStream.open(outFileName.c_str(), ios::out | ios::trunc | ios::binary);
    if (outFileStream.is_open()) {
        int i = 1;
        for (auto itr = pairs.begin(); itr != pairs.end(); ++itr) {
            index[itr->first] = i;
            outFileStream << itr->first << "\t" << itr->second << endl;
            i++;
        }
        outFileStream.close();
    } else {
        cerr << "Error: Unable to open file '" << outFileName << "' for writing." << endl;
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


void indexFile(string inFileName, map<string,int> index, char sep) {
    string outFileName = inFileName;
    outFileName.insert(inFileName.find_last_of("."), ".icorpus");
    inFileName.insert(inFileName.find_last_of("."), ".corpus");
    ifstream inFileStream;
    inFileStream.open(inFileName.c_str(), ios::in | ios::binary);
    // File open a success
    if (inFileStream.is_open()) {
        // Open new file name for output
        // size_t start = inFileName.find_first_of(".") + 1;
        // size_t end = inFileName.find_last_of(".");
        // string outFileName = inFileName.replace(start, end - start, "icorpus");
        ofstream outFileStream;
        outFileStream.open(outFileName.c_str(), ios::out | ios::trunc | ios::binary);
        if (outFileStream.is_open()) {
            string line;
            string field;
            // Keep reading until end of file
            while (getline(inFileStream, line)) {
                outFileStream << indexTokens(line, index, sep) << endl;
            }
            outFileStream.close();
        } else {
            cerr << "Error: Unable to open file '" << outFileName << "' for writing." << endl;
        }
        inFileStream.close();
    } else {
        cerr << "Error: Unable to open file '" << inFileName << "' for reading (indexFile)." << endl;
    }
}


map<string,int> readIndex(string inFileName) {
    map<string,int> index;
    ifstream inFileStream;
    inFileStream.open(inFileName.c_str(), ios::in | ios::binary);
    // File open a success
    if (inFileStream.is_open()) {
        string line;
        getline(inFileStream, line);
        int i = 1;
        while (getline(inFileStream, line)) {
            istringstream iss(line);
            string token;
            string counts;
            getline(iss, token, '\t');
            getline(iss, counts, '\t');
            index[token] = i;
            i++;
        }
        inFileStream.close();
    }
    return index;
}


template <class T>
void broadcast(T func, string chunkName, int threads) {
    string fullChunkName;
    vector<thread> threadList;
    for (int i=0;i<threads;i++) {
        // Build the filename
        fullChunkName.clear();
        fullChunkName.append(chunkName);
        fullChunkName.append(".");
        fullChunkName.append(padZeros(i + 1, threads));
        threadList.push_back(thread(func, fullChunkName));
    }
    for (auto& th : threadList) {
        th.join();
    }
}


template <class T>
void broadcast(T func, string chunkName, vector<map<string,int>> &freqs, int threads) {
    string fullChunkName;
    vector<thread> threadList;
    for (int i=0;i<threads;i++) {
        // Build the filename
        fullChunkName.clear();
        fullChunkName.append(chunkName);
        fullChunkName.append(".");
        fullChunkName.append(padZeros(i + 1, threads));
        threadList.push_back(thread(func, fullChunkName, ref(freqs[i])));
    }
    for (auto& th : threadList) {
        th.join();
    }
}


void removeTmpfiles(string chunkName, int threads) {
    string fullChunkName;
    for (int i=0;i<threads;i++) {
        fullChunkName.clear();
        fullChunkName.append(chunkName);
        fullChunkName.append(".");
        fullChunkName.append(padZeros(i + 1, threads));
        if (remove(fullChunkName.c_str()) != 0) {
            cerr << "Error: Could not remove intermediate file '" << fullChunkName << "'.";
        }
    }
}