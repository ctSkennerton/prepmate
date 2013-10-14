//
//  main.cpp
//  prepmate
//
//  Created by Connor Skennerton on 30/07/13.
//  Copyright (c) 2013 Connor Skennerton. All rights reserved.
//

#include <iostream>
#include <sstream>
#include <vector>
#include <stack>
#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>
#include <chrono>
#include <cstdlib>

#include <unistd.h>
#include <cstdio>
#include <string>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <zlib.h>

#define VERSION "0.2"

typedef enum PROCESS_DECISION {
    MP,  // read suitable for mate-pair
    PE,  // read suitable for paired-end
    R,   // remove (don't output) read
    NF   // adaptor not found
} processDecision_t;


typedef struct ALI_POS {
    unsigned int adaptorBegin;
    unsigned int adaptorEnd;
    unsigned int readBegin;
    unsigned int readEnd;
    int score;
    processDecision_t decision;
} alignmentResult_t;

struct Fx
{
    seqan::CharString id;
    seqan::CharString seq;
    seqan::CharString qual;
    
    Fx()
    {};
    
    Fx(seqan::CharString _id, seqan::CharString _seq, seqan::CharString _qual = "") :
    id(_id), seq(_seq), qual(_qual){}
};


struct Options
{
    int minLength;
    std::string outPrefix;
    seqan::CharString adaptor;
    int bandLength;
    int minScore;
    int threads;
    Options() :
    minLength(50), outPrefix("prepmate"), adaptor("CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG"), bandLength(9), minScore(15), threads(2)
    {}

};

class OutputFiles
{
private:
    std::recursive_mutex outMutexMP;
    std::recursive_mutex outMutexPE;
    std::recursive_mutex outMutexS;
    std::recursive_mutex outMutexL;
    
    gzFile matePairOutFile1;
    gzFile matePairOutFile2;
    gzFile pairedEndOutFile1;
    gzFile pairedEndOutFile2;
    gzFile singletonsOutFile;
    FILE* logFile;


public:
    OutputFiles(std::string outPrefix) {
        std::string mpfn1 = outPrefix+"_MP.1.fq.gz";
        std::string mpfn2 = outPrefix+"_MP.2.fq.gz";
        std::string pefn1 = outPrefix+"_PE.1.fq.gz";
        std::string pefn2 = outPrefix+"_PE.2.fq.gz";
        std::string sfn = outPrefix+"_singletons.fq.gz";
        std::string log_file_name = outPrefix + "_log.txt";


        matePairOutFile1 = gzopen(mpfn1.c_str(), "wb");
        matePairOutFile2 = gzopen(mpfn2.c_str(),"wb");
        pairedEndOutFile1 = gzopen(pefn1.c_str(), "wb");
        pairedEndOutFile2 = gzopen(pefn2.c_str(), "wb");
        singletonsOutFile = gzopen(sfn.c_str(), "wb");
        logFile = fopen(log_file_name.c_str(), "w");
    }

    ~OutputFiles() {
        gzclose(matePairOutFile1);
        gzclose(matePairOutFile2);
        gzclose(pairedEndOutFile1);
        gzclose(pairedEndOutFile2);
        gzclose(singletonsOutFile);
        fclose(logFile);
    }

    void printMatePair(Fx& mate1, Fx& mate2){
        std::lock_guard<std::recursive_mutex> lock(outMutexMP);
        gzprintf(matePairOutFile1, "@%s\n%s\n+\n%s\n", 
                seqan::toCString(mate1.id),
                seqan::toCString(mate1.seq),
                seqan::toCString(mate1.qual));

        gzprintf(matePairOutFile2, "@%s\n%s\n+\n%s\n", 
                seqan::toCString(mate2.id),
                seqan::toCString(mate2.seq),
                seqan::toCString(mate2.qual));
    }

    void printPairEnd(Fx& mate1, Fx& mate2){
        std::lock_guard<std::recursive_mutex> lock(outMutexPE);
        gzprintf(pairedEndOutFile1, "@%s\n%s\n+\n%s\n", 
                seqan::toCString(mate1.id),
                seqan::toCString(mate1.seq),
                seqan::toCString(mate1.qual));

        gzprintf(pairedEndOutFile2, "@%s\n%s\n+\n%s\n", 
                seqan::toCString(mate2.id),
                seqan::toCString(mate2.seq),
                seqan::toCString(mate2.qual));
    }

    void printSingletons(Fx& read){
        std::lock_guard<std::recursive_mutex> lock(outMutexS);
        gzprintf(singletonsOutFile, "@%s\n%s\n+\n%s\n", 
                seqan::toCString(read.id),
                seqan::toCString(read.seq),
                seqan::toCString(read.qual));
    
    }

    void printLog(const char * msg){
        std::lock_guard<std::recursive_mutex> lock(outMutexL);
        fprintf(logFile, "%s", msg);
    }

};
//      Constants
//
const int producer_delay_to_produce = 1;   // in milliseconds
const int consumer_delay_to_consume = 30;   // in milliseconds

const int consumer_max_wait_time = 200;     // in milliseconds - max time that a
// consumer can wait for a product to be produced.

const int max_production = 1000;              // When producers has produced this quantity
// they will stop to produce
const int max_products = 10;                // Maximum number of products that can be stored

//
//      Variables
//
std::atomic<bool> started_production(false);

std::atomic<bool> finished_production(false);       // When there's no producer working the consumers
// will stop, and the program will stop.
std::stack<std::stack<std::pair<Fx, Fx> > > products;
std::mutex xmutex;                               // Our mutex, without this mutex our program will cry
std::mutex print_mutex;

std::condition_variable is_not_full;             // to indicate that our stack is not full between
// the thread operations
std::condition_variable is_not_empty;            // to indicate that our stack is not empty between
// the thread operations


void usage(Options& opts) {
    printf("prepmate [-hVmsao] <read1.fx> <read2.fx>\n");
    //std::cout<<"\t-j           Force bzip2 formatting\n";
    //std::cout<<"\t-q           Force fastq formatting\n";
    printf("\t-h           Print this help\n");
    printf("\t-V           Print version\n");
    printf("\t-t <int>     number of threads (processors) to use. default: %d\n", opts.threads);
    printf("\t-m <int>     minimum read length after adaptor trimming. default: %d\n", opts.minLength);
    printf("\t-s <int>     minimum score for the overall alignment. default: %d\n", opts.minScore);
    printf("\t-a <DNA>     DNA sequence of the matepair double adaptor. default: %s\n", seqan::toCString(opts.adaptor));
    printf("\t-o <text>    prefix for output files. default: %s\n", opts.outPrefix.c_str());
    exit(1);
}

int parseOptions(int argc,  char * argv[], Options& opts) {
    int c;
    while ((c = getopt(argc, argv, "hVm:b:a:o:s:t:")) != -1 ) {
        switch (c) {
            case 't':
                opts.threads = atoi(optarg);
                if(opts.threads > 8) {
                    fprintf(stderr, "NOTE: The guy that wrote this doesn't have a great grasp of threading\n");
                    fprintf(stderr, "so actually any more than 8 threads gives you no speedup... changing\n");
                    fprintf(stderr, "number of threads to 8 now!\n");
                    opts.threads = 8;
                }
                break;
            case 'm':
                opts.minLength = atoi(optarg);
                break;
            case 's':
                opts.minScore = atoi(optarg);
                break;
            case 'b':
                opts.bandLength = atoi(optarg);
                break;
            case 'a':
                opts.adaptor = optarg;
                break;
            case 'o':
                opts.outPrefix = optarg;
                break;
            case 'V':
                printf("%s\n",VERSION);
                exit(1);
                break;
            case 'h':
            default:
                usage(opts);
                break;
        }
    }
    return optind;
}

void trim_seqs(Fx& mate1, 
        Fx& mate2, 
        alignmentResult_t& m1a, 
        alignmentResult_t& m2a ){

    switch(m1a.decision) {
        case MP:
            seqan::erase(mate1.seq, m1a.readBegin, seqan::length(mate1.seq));
            seqan::erase(mate1.qual, m1a.readBegin, seqan::length(mate1.qual));  
            break;
        case PE:
            seqan::erase(mate1.seq, 0, m1a.readEnd);
            seqan::erase(mate1.qual, 0, m1a.readEnd);
            break;
        default:
            break;
    }

    switch (m2a.decision) {
        case MP:
            seqan::erase(mate2.seq, m2a.readBegin, seqan::length(mate2.seq));
            seqan::erase(mate2.qual, m2a.readBegin, seqan::length(mate2.qual));
            break;
        case PE:
            seqan::erase(mate2.seq, 0, m2a.readEnd);
            seqan::erase(mate2.qual, 0, m2a.readEnd);
            break;
        default:
            break;
    }
}

void output_alignments(Fx& mate1, 
        Fx& mate2, 
        alignmentResult_t& m1a, 
        alignmentResult_t& m2a, 
        OutputFiles& outFs
        ){

    switch(m1a.decision) {
        case MP: {
            switch(m2a.decision) {
                case PE:
                    outFs.printLog("mate1 is mate-pair, mate2 is paired-end. PE trumps MP\n");
                    outFs.printPairEnd(mate1, mate2);
                    break;
                case R:
                    outFs.printLog("mate1 to singletons\n");
                    outFs.printSingletons(mate1);
                    break;
                case MP:
                default:
                    outFs.printLog("mate-pair\n");
                    outFs.printMatePair(mate1, mate2);
                    break;
            }
            break;
        }
        case PE: {
            switch(m2a.decision) {
                case PE:
                    outFs.printLog("paired-end\n");
                    outFs.printPairEnd(mate1, mate2);
                    break;
                case R:
                    outFs.printLog( "mate1 to singletons\n");
                    outFs.printSingletons(mate1);
                    break;
                case MP:
                    outFs.printLog( "mate1 is paired-end, mate2 is mate-pair. PE trumps MP\n");
                    outFs.printPairEnd(mate1, mate2);
                    break;
                default:
                    outFs.printLog( "paired-end\n");
                    outFs.printPairEnd(mate1, mate2);
                    break;
            }
            break;
        }
        case R: {
            switch(m2a.decision) {
                case PE:
                    outFs.printLog( "mate2 to singletons\n");
                    outFs.printSingletons(mate2);
                    break;
                case R:
                    outFs.printLog( "both removed\n");
                    break;
                case MP:
                    outFs.printLog( "mate2 to singletons\n");
                    outFs.printSingletons(mate2);
                    break;
                default:
                    outFs.printLog( "mate2 to singletons\n");
                    outFs.printSingletons(mate2);
                    break;
            }
            break;
        }
        default: {
            switch(m2a.decision) {
                case PE:
                    outFs.printLog( "paired-end\n");
                    outFs.printPairEnd(mate1, mate2);
                    break;
                case R:
                    outFs.printLog( "mate1 to singletons\n");
                    outFs.printSingletons(mate1);
                    break;
                case MP:
                    outFs.printLog( "mate-pair\n");
                    outFs.printMatePair(mate1, mate2);
                    break;
                default:
                    outFs.printLog( "mate-pair\n");
                    outFs.printMatePair(mate1, mate2);
                    break;
            }
            break;
        }

    }
}

void log_alignment(seqan::CharString& name1, 
        seqan::CharString& name2, 
        alignmentResult_t& mate1_alignment, 
        alignmentResult_t& mate2_alignment,
        OutputFiles& outFs ) {
    
    // print a summary of the alignments
    char buffer[1024];
    sprintf(buffer, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t", 
            seqan::toCString(name1), 
            mate1_alignment.score, 
            mate1_alignment.adaptorBegin, 
            mate1_alignment.adaptorEnd, 
            mate1_alignment.readBegin, 
            mate1_alignment.readEnd, 
            mate1_alignment.decision,
            seqan::toCString(name2), 
            mate2_alignment.score, 
            mate2_alignment.adaptorBegin, 
            mate2_alignment.adaptorEnd, 
            mate2_alignment.readBegin, 
            mate2_alignment.readEnd, 
            mate2_alignment.decision);

    outFs.printLog(buffer);
}

void align_pair(Fx& mate1, 
        Fx& mate2, 
        seqan::Align< seqan::String<char> >& ali, 
        const int minLength, 
        const int minScore,
        OutputFiles& outFs
        ) {

    
    alignmentResult_t mate1_alignment;
    alignmentResult_t mate2_alignment;
    
    // add in mate 1 to the alignment
    seqan::assignSource(seqan::row(ali, 1), mate1.seq);
    mate1_alignment.score = seqan::localAlignment(ali, seqan::Score<int>(1,-1,-5, -4));

    mate1_alignment.adaptorBegin = seqan::clippedBeginPosition(seqan::row(ali, 0));
    mate1_alignment.adaptorEnd = seqan::clippedEndPosition(seqan::row(ali, 0)) - 1;
    mate1_alignment.readBegin = seqan::clippedBeginPosition(seqan::row(ali, 1));
    mate1_alignment.readEnd = seqan::clippedEndPosition(seqan::row(ali, 1)) - 1;


    seqan::assignSource(seqan::row(ali, 1), mate2.seq);
    mate2_alignment.score = seqan::localAlignment(ali, seqan::Score<int>(1,-1,-5, -4));

    mate2_alignment.adaptorBegin = seqan::clippedBeginPosition(seqan::row(ali, 0));
    mate2_alignment.adaptorEnd = seqan::clippedEndPosition(seqan::row(ali, 0)) - 1;
    mate2_alignment.readBegin = seqan::clippedBeginPosition(seqan::row(ali, 1));
    mate2_alignment.readEnd = seqan::clippedEndPosition(seqan::row(ali, 1)) - 1;

    
    if(mate1_alignment.score >= minScore) {
        if(mate1_alignment.readBegin >= minLength) {
            // this read is long enough to be a mate-pair
            mate1_alignment.decision = MP;
        } else if(seqan::length(mate1.seq) - 1 - mate1_alignment.readEnd >= minLength) {
            // this read can make a good paired-end
            mate1_alignment.decision = PE;
        } else {
            // remove this read altogether
            mate1_alignment.decision = R;
        }
    } else {
        mate1_alignment.decision = NF;
    }
    
    if(mate2_alignment.score >= minScore) {
        if(mate2_alignment.readBegin >= minLength) {
            // this read is long enough to be a mate-pair
            mate2_alignment.decision = MP;
        } else if(seqan::length(mate2.seq) - 1 - mate2_alignment.readEnd >= minLength) {
            // this read can make a good paired-end
            mate2_alignment.decision = PE;
        } else {
            // remove this read altogether
            mate2_alignment.decision = R;
        }
    } else {
        mate2_alignment.decision = NF;
    }
    log_alignment(mate1.id, mate2.id, mate1_alignment, mate2_alignment, outFs);
    trim_seqs(mate1, mate2, mate1_alignment, mate2_alignment);
    output_alignments(mate1, 
        mate2, 
        mate1_alignment, 
        mate2_alignment,
        outFs
        );
}
void consume(int consumer_id, seqan::Align< seqan::String<char> >& alignment, const int minLength, const int minScore, OutputFiles& outFs)
{
    std::stack<std::pair<Fx, Fx> > product;
    {
        std::unique_lock<std::mutex> lock(xmutex);
        
        if(is_not_empty.wait_for(lock, std::chrono::milliseconds(consumer_max_wait_time),
                                 [] { return products.size() > 0; }))
        {
            product = products.top();
            products.pop();
            is_not_full.notify_all();
        }
    }
    std::pair<Fx,Fx> readPair;
    while(! product.empty()) {
        readPair = product.top();
        product.pop();
        align_pair(readPair.first, 
            readPair.second, 
            alignment, 
            minLength,
            minScore,
            outFs
            );
    }
}

void worker(int threadId, seqan::CharString adaptorSeq, int minLength, int minScore, OutputFiles& outFs)
{
    // Wait until there is any producer working
    while(! started_production ) std::this_thread::yield();

    // make our alignment object
    seqan::Align< seqan::String<char> > ali;
    seqan::resize(rows(ali), 2);
    seqan::assignSource(seqan::row(ali, 0), adaptorSeq);
    
    while(!finished_production || products.size() > 0)
    {
        consume(threadId, ali, minLength, minScore, outFs);
        std::this_thread::sleep_for(std::chrono::milliseconds(consumer_delay_to_consume));
    }

}
std::thread::id main_thread_id;
int main(int argc, char * argv[])
{
    main_thread_id = std::this_thread::get_id();
    Options opts;

    int opt_idx = parseOptions(argc, argv, opts);

    if (opt_idx > argc - 1) {
        fprintf(stderr, "%s\n", "Please provide input files");
        usage(opts);
    }
    
    seqan::SequenceStream read1(argv[opt_idx++]);
    if (!seqan::isGood(read1))
        fprintf(stderr, "%s\n", "Could not open read1 file");
    
    // Read one record.
    if (opt_idx >= argc) {
        fprintf(stderr, "%s\n", "only one input file has been provided, please give two");
        usage(opts);
    }
    // we have a mate file
    seqan::SequenceStream read2(argv[opt_idx]);
    if (!seqan::isGood(read2))
        fprintf(stderr, "%s\n", "Could not open read2 file");


    OutputFiles outFs(opts.outPrefix);
    char buffer[300];
    sprintf(buffer, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
        "mate1", "score1", "adaptor_begin1", "adaptor_end1", "read_begin1", "read_end1", "type1",
        "mate2", "score2", "adaptor_begin2", "adaptor_end2", "read_begin2", "read_end2", "type2", "output");

    outFs.printLog(buffer);

    std::vector<std::thread> consumers;

    // Create consumers
    for(int i = 0; i < opts.threads; ++i)
        consumers.push_back(std::thread(worker, i, opts.adaptor, opts.minLength, opts.minScore, std::ref(outFs)));
    
    if (std::this_thread::get_id() == main_thread_id) {

        seqan::StringSet<seqan::CharString> ids1;
        seqan::StringSet<seqan::CharString> seqs1;
        seqan::StringSet<seqan::CharString> quals1;
        seqan::StringSet<seqan::CharString> ids2;
        seqan::StringSet<seqan::CharString> seqs2;
        seqan::StringSet<seqan::CharString> quals2;


        while (! seqan::atEnd(read1)) {

            if (seqan::atEnd(read2)) {
                fprintf(stderr, "%s\n","files have different number of reads");
                break;
            }

            if (readBatch(ids1, seqs1, quals1, read1, max_production) != 0) {
                std::unique_lock<std::mutex> lock(print_mutex);
                fprintf(stderr, "Malformed record in file 1\n");
                break;
            }

            if (seqan::readBatch(ids2, seqs2, quals2, read2, max_production) != 0) {
                std::unique_lock<std::mutex> lock(print_mutex);
                fprintf(stderr, "Malformed record in file 2");
                break;
            }
            if(seqan::length(ids1) != seqan::length(ids2)) {
                std::unique_lock<std::mutex> lock(print_mutex);
                fprintf(stderr, "different number of reads in pair files");
                break;
            }
            std::stack<std::pair<Fx, Fx> > readPairs;

            typedef seqan::Iterator<seqan::StringSet<seqan::CharString> >::Type TStringSetIterator;
            TStringSetIterator it_id = seqan::begin(ids1);
            TStringSetIterator it_seq = seqan::begin(seqs1);
            TStringSetIterator it_qual = seqan::begin(quals1);
            TStringSetIterator it_id2 = seqan::begin(ids2);
            TStringSetIterator it_seq2 = seqan::begin(seqs2);
            TStringSetIterator it_qual2 = seqan::begin(quals2);
            for (; it_id != seqan::end(ids1); ++it_id, ++it_seq, ++it_qual, ++it_id2, ++it_seq2, ++it_qual2) {
                Fx current_read1;
                Fx current_read2;
                current_read1.id = value(it_id);
                current_read1.seq = value(it_seq);
                current_read1.qual = value(it_qual);
                current_read2.id = value(it_id2);
                current_read2.seq = value(it_seq2);
                current_read2.qual = value(it_qual2);
                readPairs.push(std::make_pair(current_read1, current_read2));
            }
            {
                std::unique_lock<std::mutex> lock(xmutex);
                is_not_full.wait(lock, [] { return products.size() <= max_products; });
                products.push(readPairs);
                started_production = true;
            }
            seqan::clear(ids1);
            seqan::clear(ids2);
            seqan::clear(seqs1);
            seqan::clear(seqs2);
            seqan::clear(quals1);
            seqan::clear(quals2);
        }
        finished_production = true;
    }
    
    // Wait for consumers and producers to finish
    for(auto& t : consumers)
        t.join();

    
    return 0;
}
