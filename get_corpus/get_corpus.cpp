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

static int quiet_flag, rcomp_flag;
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
        string header;

        getline(fileStream, header);
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
        if (!quiet_flag) {
            cout << endl << "Chunking complete! " << counter - 1 << " files created." << endl;
        }
    } else {
        cerr << "Error: Unable to open file '" << fullFilePath << "' for reading." << endl;
    }
}


// Finds chunks by "chunkName" and creates file specified in outputFileName
void joinFile(string chunkName, string outputFileName, string header) {
    string inputFileName;

    // Create our output file
    ofstream fileOutput;
    fileOutput.open(outputFileName.c_str(), ios::out | ios::binary);

    // If successful, loop through chunks matching chunkName
    if (fileOutput.is_open()) {
        bool filefound = true;
        int counter = 1;
        int fileSize = 0;

        fileOutput << header << endl;

        while (filefound) {
            filefound = false;

            // Build the filename
            inputFileName.clear();
            inputFileName.append(chunkName);
            inputFileName.append(".");
            inputFileName.append(to_string(counter));
            inputFileName.append(".corpus");

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

        if (!quiet_flag) {
            char resolved_path[PATH_MAX];
            cout << endl << "File '" << realpath(&outputFileName[0], resolved_path) << "' assembly complete!" << endl;
        }
    } else {
        cerr << "Error: Unable to open file '" << outputFileName << "' for writing." << endl;
    }

}


string translate(string token, string alphabet, string complement) {
    transform(token.begin(), token.end(), token.begin(), ::toupper);
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


string tokenizeSeq(string seq, map<string,string> &rev_comp, int k, int stride, string alphabet, string complement, string delim) {
    string origSeq = seq;
    string kmer = seq.substr(0, k);
    string tokens = translate(kmer, alphabet, complement);
    for (int i=stride;i<seq.length();i+=stride) {
        string kmer = seq.substr(i, k);
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
        tokens.append(delim);
        tokens.append(token);
    }
    return tokens;
}


// Tokenizes a sequence tsv file.
void tokenizeFile(string fullFilePath, int k, int stride, string alphabet, string complement, string delim) {
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
                for (int i=0;i<6;i++) {
                    getline(iss, field, '\t');
                    output << field << '\t';
                }
                getline(iss, field, '\t');
                output << tokenizeSeq(field, rev_comp, k, stride, alphabet, complement, delim) << endl;
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
        string alphabet, complement, delim;
        int k, stride;
    public:
        fileTokenizer(int k1, int stride1, string alphabet1, string complement1, string delim1) {
            k = k1;
            stride = stride1;
            alphabet = alphabet1;
            complement = complement1;
            delim = delim1;
        }
        void operator()(string fullFilePath) { 
            tokenizeFile(fullFilePath, k, stride, alphabet, complement, delim);
        }
};


void usage() {
    const char *text = "usage: get_corpus [-acdtoh] [--verbose] [--rcomp] ... [Sequence Corpus File] [kmer size] [stride length]\n"
    "\n"
    "       Tokenize sequence corpus.\n"
    "\n"
    "       Options:\n"
    "           -a    Sequence alphabet.\n"
    "           -c    Complement alphabet.\n"
    "           -d    Kmer token delimiter.\n"
    "           -t    Number of threads to run on.\n"
    "           -o    Output file path.\n"
    "           -h    Display this help and exit.\n\n";
    cout << text;
    exit(1);
}


int main(int argc, char **argv) {
    string alphabet = "ACGTNBDHU";
    string complement = "TGCANUHDB";
    string delim = " ";
    int threads = MAX_THREADS;
    string outFile = "corpus.tsv";
    int c;
    while (1) {
        static struct option long_options[] = {
            /* These options set a flag. */
            {"quiet", no_argument, &quiet_flag, 1},
            // Match rcomp as well
            {"rcomp", no_argument, &rcomp_flag, 1},
            /* These options don’t set a flag.
            We distinguish them by their indices. */
            {"alphabet", required_argument, 0, 'a'},
            {"complement", required_argument, 0, 'c'},
            {"delim", required_argument, 0, 'd'},
            {"threads", required_argument, 0, 't'},
            {"out", required_argument, 0, 'o'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long(argc, argv, "a:c:d:t:o:h",
                       long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1) {
            break;
        }

        switch (c) {
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
                    usage();
                }
                break;

            case 'c':
                if (optarg) {
                    complement = optarg;
                } else {
                    usage();
                }
                break;

            case 'd':
                if (optarg) {
                    delim = optarg;
                } else {
                    usage();
                }
                break;

            case 't':
                if (optarg) {
                    int new_threads = stoi(optarg);
                    if (new_threads <= MAX_THREADS) {
                        threads = new_threads;
                    }
                } else {
                    usage();
                }
                break;

            case 'o':
                if (optarg) {
                    outFile = optarg;
                } else {
                    usage();
                }
                break;

            case 'h':
                usage();
                break;

            case '?':
                /* getopt_long already printed an error message. */
                break;

            default:
                usage();
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
    ifstream fileStream;
    fileStream.open(fullFilePath.c_str(), ios::in | ios::binary);

    int lines = 0;
    if (fileStream.is_open()) {
        string line;
        getline(fileStream, header);
        while (getline(fileStream, line)) {
            lines++;
        }
        fileStream.close();
    }

    if (!quiet_flag) {
        char resolved_path[PATH_MAX];
        cout << "The file '" << realpath(&fullFilePath[0], resolved_path) << "' has " << lines << " lines.\n";
        cout << "There are " << threads << " out of " << MAX_THREADS << " threads in use.\n";
    }

    int chunkSize = ceil((float)lines / threads);

    if (mkdir(tmp_dir.c_str(), 0777) && (errno != EEXIST)) {
        cerr << "Error: Unable to create directory '" << tmp_dir << "'.";
    }

    chunkFile(fullFilePath, chunkName, chunkSize);

    fileTokenizer tokenizer(k, stride, alphabet, complement, delim);
    vector<thread> threadList;

    for (int i=0;i<threads;i++) {
        fullChunkName.clear();
        fullChunkName.append(chunkName);
        fullChunkName.append(".");
        fullChunkName.append(to_string(i + 1));
        // cout << fullChunkName << endl;
        threadList.push_back(thread(tokenizer, fullChunkName));
        // threadList.push_back(thread(tokenizeFile, fullChunkName, 5, 1, "ACGTNBDHU", "TGCANUHDB", " "));
    }
    for (auto& th : threadList) {
        th.join();
    }
    if (!quiet_flag) {
        cout << endl << "Tokenizing complete!" << endl;
    }

    for (int i=0;i<threads;i++) {
        fullChunkName.clear();
        fullChunkName.append(chunkName);
        fullChunkName.append(".");
        fullChunkName.append(to_string(i + 1));
        if (remove(fullChunkName.c_str()) != 0) {
            cerr << "Error: Could not remove intermediate file '" << fullChunkName << "'.";
        }
    }

    joinFile(chunkName, outFile, header);

    for (int i=0;i<threads;i++) {
        fullTokenizedChunkName.clear();
        fullTokenizedChunkName.append(chunkName);
        fullTokenizedChunkName.append(".");
        fullTokenizedChunkName.append(to_string(i + 1));
        fullTokenizedChunkName.append(".corpus");
        if (remove(fullTokenizedChunkName.c_str()) != 0) {
            cerr << "Error: Could not remove intermediate file '" << fullTokenizedChunkName << "'.";
        }
    }

    if (rmdir(tmp_dir.c_str())) {
        cerr << "Error: Unable to remove directory '" << tmp_dir << "'.";
    }
}