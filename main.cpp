//
//  main.cpp
//  prepmate
//
//  Created by Connor Skennerton on 30/07/13.
//  Copyright (c) 2013 Connor Skennerton. All rights reserved.
//


#include <unistd.h>
#include <cstdlib>
#include <cstdio>
#include <string>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <zlib.h>

#define VERSION "0.1"
using namespace seqan;

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
    Options() :
    minLength(50), outPrefix("prepmate"), adaptor("CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG"), bandLength(9), minScore(15)
    {}

};


void printSingle(Fx& mate1, gzFile out ) {
    gzprintf(out, "@%s\n%s\n+\n%s\n", 
            seqan::toCString(mate1.id),
            seqan::toCString(mate1.seq),
            seqan::toCString(mate1.qual));
    //seqan::writeRecord(out, mate1.id, mate1.seq, mate1.qual, seqan::Fastq());
}

void printPair(Fx& mate1, Fx& mate2, gzFile out1, gzFile out2) {
    // match in the first read print out pair
    printSingle(mate1, out1);
    printSingle(mate2, out2);
}

void usage(Options& opts) {
    printf("prepmate [-hVmsao] <read1.fx> <read2.fx>\n");
    //std::cout<<"\t-j           Force bzip2 formatting\n";
    //std::cout<<"\t-q           Force fastq formatting\n";
    printf("\t-h           Print this help\n");
    printf("\t-V           Print version\n");
    printf("\t-m <int>     minimum read length after adaptor trimming. default: %d\n", opts.minLength);
    printf("\t-b <int>     minimum length for the initial match to the adaptor. default: %d\n", opts.bandLength);
    printf("\t-s <int>     minimum score for the overall alignment. default: %d\n", opts.minScore);
    printf("\t-a <DNA>     DNA sequence of the matepair double adaptor. default: %s\n", seqan::toCString(opts.adaptor));
    printf("\t-o <text>    prefix for output files. default: %s\n", opts.outPrefix.c_str());
    exit(1);
}

int parseOptions(int argc,  char * argv[], Options& opts) {
    int c;
    while ((c = getopt(argc, argv, "hVm:b:a:o:s:")) != -1 ) {
        switch (c) {
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
        gzFile matePairOutFile1, 
        gzFile matePairOutFile2, 
        gzFile pairedEndOutFile1,
        gzFile pairedEndOutFile2,
        gzFile singletonsOutFile,
        FILE* logFile){

    switch(m1a.decision) {
        case MP: {
            switch(m2a.decision) {
                case PE:
                    fprintf(logFile, "%s\n", "mate1 is mate-pair, mate2 is paired-end. PE trumps MP");
                    printPair(mate1, mate2, pairedEndOutFile1, pairedEndOutFile2);
                    break;
                case R:
                    fprintf(logFile, "%s\n", "mate1 to singletons");
                    printSingle(mate1, singletonsOutFile);
                    break;
                case MP:
                default:
                    fprintf(logFile, "%s\n", "mate-pair");
                    printPair(mate1, mate2, matePairOutFile1, matePairOutFile2);
                    break;
            }
            break;
        }
        case PE: {
            switch(m2a.decision) {
                case PE:
                    fprintf(logFile, "%s\n", "paired-end");
                    printPair(mate1, mate2, pairedEndOutFile1, pairedEndOutFile2);
                    break;
                case R:
                    fprintf(logFile, "%s\n", "mate1 to singletons");
                    printSingle(mate1, singletonsOutFile);
                    break;
                case MP:
                    fprintf(logFile, "%s\n", "mate1 is paired-end, mate2 is mate-pair. PE trumps MP");
                    printPair(mate1, mate2, pairedEndOutFile1, pairedEndOutFile2);
                    break;
                default:
                    fprintf(logFile, "%s\n", "paired-end");
                    printPair(mate1, mate2, pairedEndOutFile1, pairedEndOutFile2);
                    break;
            }
            break;
        }
        case R: {
            switch(m2a.decision) {
                case PE:
                    fprintf(logFile, "%s\n", "mate2 to singletons");
                    printSingle(mate2, singletonsOutFile);
                    break;
                case R:
                    fprintf(logFile, "%s\n", "both removed");
                    break;
                case MP:
                    fprintf(logFile, "%s\n", "mate2 to singletons");
                    printSingle(mate2, singletonsOutFile);
                    break;
                default:
                    fprintf(logFile, "%s\n", "mate2 to singletons");
                    printSingle(mate2, singletonsOutFile);
                    break;
            }
            break;
        }
        default: {
            switch(m2a.decision) {
                case PE:
                    fprintf(logFile, "%s\n", "paired-end");
                    printPair(mate1, mate2, pairedEndOutFile1, pairedEndOutFile2);
                    break;
                case R:
                    fprintf(logFile, "%s\n", "mate1 to singletons");
                    printSingle(mate1, singletonsOutFile);
                    break;
                case MP:
                    fprintf(logFile, "%s\n", "mate-pair");
                    printPair(mate1, mate2, matePairOutFile1, matePairOutFile2);
                    break;
                default:
                    fprintf(logFile, "%s\n", "mate-pair");
                    printPair(mate1, mate2, matePairOutFile1, matePairOutFile2);
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
        FILE* logFile ) {
    
    // print a summary of the alignments
    fprintf(logFile, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t", 
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
    
}

void align_pair(Fx& mate1, 
        Fx& mate2, 
        Align< seqan::String<char> >& ali, 
        int minScore, 
        int minLength,
        gzFile matePairOutFile1, 
        gzFile matePairOutFile2, 
        gzFile pairedEndOutFile1,
        gzFile pairedEndOutFile2,
        gzFile singletonsOutFile,
        FILE* logFile) {

    
    alignmentResult_t mate1_alignment;
    alignmentResult_t mate2_alignment;
    
    // add in mate 1 to the alignment
    assignSource(row(ali, 1), mate1.seq);
    mate1_alignment.score = localAlignment(ali, Score<int>(1,-1,-5, -4));

    mate1_alignment.adaptorBegin = clippedBeginPosition(row(ali, 0));
    mate1_alignment.adaptorEnd = clippedEndPosition(row(ali, 0)) - 1;
    mate1_alignment.readBegin = clippedBeginPosition(row(ali, 1));
    mate1_alignment.readEnd = clippedEndPosition(row(ali, 1)) - 1;


    assignSource(row(ali, 1), mate2.seq);
    mate2_alignment.score = localAlignment(ali, Score<int>(1,-1,-5, -4));

    mate2_alignment.adaptorBegin = clippedBeginPosition(row(ali, 0));
    mate2_alignment.adaptorEnd = clippedEndPosition(row(ali, 0)) - 1;
    mate2_alignment.readBegin = clippedBeginPosition(row(ali, 1));
    mate2_alignment.readEnd = clippedEndPosition(row(ali, 1)) - 1;

    
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
    log_alignment(mate1.id, mate2.id, mate1_alignment, mate2_alignment, logFile);
    trim_seqs(mate1, mate2, mate1_alignment, mate2_alignment);
    
    output_alignments(mate1, 
        mate2, 
        mate1_alignment, 
        mate2_alignment,
        matePairOutFile1, 
        matePairOutFile2, 
        pairedEndOutFile1,
        pairedEndOutFile2,
        singletonsOutFile,
        logFile);
}

int main(int argc, char * argv[])
{
    
    Options opts;

    int opt_idx = parseOptions(argc, argv, opts);

    if (opt_idx > argc - 1) {
        fprintf(stderr, "%s\n", "Please provide input files");
        usage(opts);
    }
    
    seqan::SequenceStream read1(argv[opt_idx++]);
    if (!isGood(read1))
        fprintf(stderr, "%s\n", "Could not open read1 file");
    // Read one record.
    Fx mate1 = Fx();
    if (opt_idx >= argc) {
        fprintf(stderr, "%s\n", "only one input file has been provided, please give two");
        usage(opts);
    }
    // we have a mate file
    seqan::SequenceStream read2(argv[opt_idx]);
    if (!isGood(read2))
        fprintf(stderr, "%s\n", "Could not open read2 file");
    Fx mate2 = Fx();

    // make our alignment object
    Align< seqan::String<char> > ali;
    resize(rows(ali), 2);
    assignSource(row(ali, 0), opts.adaptor);

    seqan::CharString mpfn1 = opts.outPrefix+"_MP.1.fq.gz";
    seqan::CharString mpfn2 = opts.outPrefix+"_MP.2.fq.gz";
    seqan::CharString pefn1 = opts.outPrefix+"_PE.1.fq.gz";
    seqan::CharString pefn2 = opts.outPrefix+"_PE.2.fq.gz";
    seqan::CharString sfn = opts.outPrefix+"_singletons.fq.gz";
    seqan::CharString log_file_name = opts.outPrefix + "_log.txt";


    gzFile matePairOutFile1 = gzopen(seqan::toCString(mpfn1), "wb");
    gzFile matePairOutFile2 = gzopen(seqan::toCString(mpfn2),"wb");
    gzFile pairedEndOutFile1 = gzopen(seqan::toCString(pefn1), "wb");
    gzFile pairedEndOutFile2 = gzopen(seqan::toCString(pefn2), "wb");
    gzFile singletonsOutFile = gzopen(seqan::toCString(sfn), "wb");
    FILE* log_file_fp = fopen(seqan::toCString(log_file_name), "w");

    fprintf(log_file_fp, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
        "mate1", "score1", "adaptor_begin1", "adaptor_end1", "read_begin1", "read_end1", "type1",
        "mate2", "score2", "adaptor_begin2", "adaptor_end2", "read_begin2", "read_end2", "type2", "output");
    
    int counter1 = 0;
    while (!atEnd(read1))
    {
        if (atEnd(read2)) {
            fprintf(stderr, "%s\n","files have different number of reads");
            break;
        }

        
        if (readRecord(mate1.id, mate1.seq, mate1.qual, read1) != 0) {
            fprintf(stderr, "%s: read %d\n", "Malformed record in file 1", counter1);
        }

        if (readRecord(mate2.id, mate2.seq, mate2.qual, read2) != 0) {
            fprintf(stderr, "%s: read %d\n", "Malformed record in file 2", counter1);
        }
        counter1++;
        if((counter1 % 1000) == 0) {
            printf("%d\n", counter1);
        }
        align_pair(mate1, 
            mate2, 
            ali, 
            opts.minScore, 
            opts.minLength,
            matePairOutFile1,
            matePairOutFile2,
            pairedEndOutFile1,
            pairedEndOutFile2,
            singletonsOutFile,
            log_file_fp);
 
    }
    gzclose(matePairOutFile1);
    gzclose(matePairOutFile2);
    gzclose(pairedEndOutFile1);
    gzclose(pairedEndOutFile2);
    gzclose(singletonsOutFile);
    fclose(log_file_fp);
    return 0;
}

