//
//  main.cpp
//  prepmate
//
//  Created by Connor Skennerton on 30/07/13.
//  Copyright (c) 2013 Connor Skennerton. All rights reserved.
//

#include <fstream>
#include <iostream>
#include <set>
#include <vector>
#include <string>
#include <unistd.h>
#include <cstdlib>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/align.h>


#define VERSION "0.1"
using namespace seqan;

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
    seqan::CharString outPrefix;
    seqan::CharString adaptor;
    int bandLength;
    int minScore;
    Options() :
    minLength(50), outPrefix("prepmate"), adaptor("CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG"), bandLength(9), minScore(9)
    {}

};



void printSingle(Fx& mate1, std::ostream& out ) {
    if (seqan::empty(mate1.qual)) {
        //out<<">"<<mate1.id<<std::endl;
        //out<<mate1.seq<<std::endl;
        seqan::writeRecord(out, mate1.id, mate1.seq, seqan::Fasta());
        
    } else {
        //std::cout<<"@"<<mate1.id<<'\n'<<mate1.seq<<"\n+\n"<<mate1.qual<<std::endl;
        seqan::writeRecord(out, mate1.id, mate1.seq, mate1.qual, seqan::Fastq());
    }
}

void printPair(Fx& mate1, Fx& mate2, std::ostream& out1, std::ostream& out2) {
    // match in the first read print out pair
    printSingle(mate1, out1);
    printSingle(mate2, out2);
}

void usage(Options& opts) {
    std::cout<< "prepmate [-hVmsao] <read1.fx> <read2.fx>\n";
    //std::cout<<"\t-j           Force bzip2 formatting\n";
    //std::cout<<"\t-q           Force fastq formatting\n";
    std::cout<<"\t-h           Print this help"<<std::endl;
    std::cout<<"\t-V           Print version"<<std::endl;
    std::cout<<"\t-m <int>     minimum read length after adaptor trimming. default: "<<opts.minLength <<std::endl;
    std::cout<<"\t-b <int>     minimum length for the initial match to the adaptor. default: "<<opts.bandLength <<std::endl;
    std::cout<<"\t-s <int>     minimum score for the overall alignment. default: "<<opts.minScore <<std::endl;
    std::cout<<"\t-a <DNA>     DNA sequence of the matepair double adaptor. default: "<<opts.adaptor<<std::endl;
    std::cout<<"\t-o <text>    prefix for output files. default: "<<opts.outPrefix <<std::endl;
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
                std::cout <<VERSION<<std::endl;
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

int main(int argc, char * argv[])
{
    
    Options opts;

    int opt_idx = parseOptions(argc, argv, opts);

    if (opt_idx > argc - 1) {
        std::cout << "Please provide input files"<<std::endl;
        usage(opts);
    }
    
    seqan::SequenceStream read1(argv[opt_idx++]);
    if (!isGood(read1))
        std::cerr << "Could not open read1 file\n";
    // Read one record.
    Fx mate1 = Fx();
    if (opt_idx >= argc) {
        std::cerr << "only one input file has been provided, please give two"<<std::endl;
        usage(opts);
    }
    // we have a mate file
    seqan::SequenceStream read2(argv[opt_idx]);
    if (!isGood(read2))
        std::cerr << "Could not open read2 file\n";
    Fx mate2 = Fx();

    // make our alignment object
    Align< seqan::String<char> > ali;
    resize(rows(ali), 2);
    assignSource(row(ali, 0), opts.adaptor);

    int counter = 0;
    int counter1 = 0;
    while (!atEnd(read1))
    {
        if (atEnd(read2)) {
            std::cerr<< "files have different number of reads"<<std::endl;
            break;
        }

        
        if (readRecord(mate1.id, mate1.seq, mate1.qual, read1) != 0) {
            std::cerr<<"Malformed record"<<std::endl;
        }

        if (readRecord(mate2.id, mate2.seq, mate2.qual, read2) != 0) {
            std::cerr<<"Malformed record"<<std::endl;
        }

        /*if((counter1 % 100) == 0) {
            std::cout<<counter1<<std::endl;
        }*/

        unsigned int cBeginPos0;
        unsigned int cEndPos0;
        unsigned int cBeginPos1;
        unsigned int cEndPos1;
        // add in mate 1 to the alignment
        assignSource(row(ali, 1), mate1.seq);
        //int ali_score = localAlignment(ali, Score<int>(1,-1,-5, -4), -opts.bandLength, opts.bandLength);
        int ali_score = localAlignment(ali, Score<int>(1,-1,-5, -4));
        // check to see if either read has got a match above the seed length
        if(ali_score >= opts.minScore) {
            cBeginPos0 = clippedBeginPosition(row(ali, 0));
            cEndPos0 = clippedEndPosition(row(ali, 0)) - 1;
            cBeginPos1 = clippedBeginPosition(row(ali, 1));
            cEndPos1 = clippedEndPosition(row(ali, 1)) - 1;
            std::cout<< cBeginPos0<<" : "<<cBeginPos1<<std::endl;
            std::cout <<mate1.id<<std::endl;
            std::cout<<"adaptor = "<<cBeginPos0<<" : "<<cEndPos0<<std::endl;
            std::cout<<"seq     = "<<cBeginPos1<<" : "<<cEndPos1<<std::endl;
            std::cout <<ali<<std::endl;
        }

        assignSource(row(ali, 1), mate2.seq);
        //ali_score = localAlignment(ali, Score<int>(1,-1,-5, -4), -opts.bandLength, opts.bandLength);
        ali_score = localAlignment(ali, Score<int>(1,-1,-5, -4));
        if(ali_score >= opts.minScore) {
        // read 1 has a seed match
            cBeginPos0 = clippedBeginPosition(row(ali, 0));
            cEndPos0 = clippedEndPosition(row(ali, 0)) - 1;
            cBeginPos1 = clippedBeginPosition(row(ali, 1));
            cEndPos1 = clippedEndPosition(row(ali, 1)) - 1;
            std::cout<< cBeginPos0<<" : "<<cBeginPos1<<std::endl;
            std::cout <<mate2.id<<std::endl;
            std::cout<<"adaptor = "<<cBeginPos0<<" : "<<cEndPos0<<std::endl;
            std::cout<<"seq     = "<<cBeginPos1<<" : "<<cEndPos1<<std::endl;
            std::cout <<ali<<std::endl;
        }
        counter1++;
    }

    return 0;
}

