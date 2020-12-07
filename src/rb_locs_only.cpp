#include <cstdio>
#include <getopt.h>
#include <string>
extern "C" {
#include <zlib.h>
#ifndef AC_KSEQ_H
#include "kseq.h"
KSEQ_INIT(gzFile, gzread);
#endif
}
#include "rowbowt.hpp"
#include "rowbowt_io.hpp"

struct RbAlignArgs {
    std::string inpre = "";
    std::string fastq_fname = "";
    std::string outpre = "";
    size_t max_hits = -1;
    size_t wsize = 10;
};

void print_help() {
    fprintf(stderr, "rb_locs_only");
    fprintf(stderr, "Usage: rb_locs_only [options] <index_prefix> <input_fastq_name>\n");
    fprintf(stderr, "    --output_prefix/-o <basename>    output prefix\n");
    fprintf(stderr, "    <input_prefix>                   index prefix\n");
    fprintf(stderr, "    <input_fastq>                    input fastq\n");
}

RbAlignArgs parse_args(int argc, char** argv) {
    int c;
    char* end;
    RbAlignArgs args;
    static struct option long_options[] {
        {"wsize", required_argument, 0, 'w'},
        {"output_prefix", required_argument, 0, 'o'},
        {"max-hits", required_argument, 0, 'm'},
    };
    int long_index = 0;
    while((c = getopt_long(argc, argv, "o:w:m:h", long_options, &long_index)) != -1) {
        switch (c) {
            case 'w':
                args.wsize = std::atol(optarg);
                break;
            case 'm':
                args.max_hits = std::atol(optarg);
                break;
            case 'o':
                args.outpre = optarg;
                break;
            case 'h':
                print_help();
                exit(0);
                break;
            default:
                print_help();
                exit(1);
                break;
        }
    }

    if (argc - optind < 2) {
        fprintf(stderr, "no argument provided\n");
        exit(1);
    }

    args.inpre = argv[optind++];
    args.fastq_fname = argv[optind++];
    if (args.outpre == "") {
        args.outpre = args.inpre;
    }
    return args;
}


void rb_report(const rbwt::RowBowt& rbwt, const rle_window_arr<>& midx, const RbAlignArgs args, kseq_t* seq, std::vector<uint64_t>& locs) {
    locs.clear();
    locs = rbwt.find_locs_greedy_seeding(seq->seq.s, args.wsize, args.max_hits, locs);
    std::cout << seq->name.s;
    for (auto l: locs) {
        auto pair = rbwt.resolve_offset(l);
        auto ms = midx.at_range(l, l+seq->seq.l-1);
        for (auto m: ms) {
            std::cout << " " << get_pos(m) << "/" << static_cast<uint64_t>(get_allele(m));
        }
    } 
    std::cout << "\n";
}

rbwt::RowBowt load_rbwt(const RbAlignArgs args) {
    rbwt::LoadRbwtFlag flag;
    flag = flag | rbwt::LoadRbwtFlag::DL | rbwt::LoadRbwtFlag::SA;
    rbwt::RowBowt rbwt(rbwt::load_rowbowt(args.inpre, flag));
    return rbwt;
}

void rb_locs_all(RbAlignArgs args) {
    rbwt::RowBowt rbwt(load_rbwt(args));
    rle_window_arr<> midx;
    std::ifstream midx_ifs(args.inpre + ".midx"); 
    midx.load(midx_ifs);
    int err, nreads = 0, noccs = 0;
    gzFile fq_fp(gzopen(args.fastq_fname.data(), "r"));
    if (fq_fp == NULL) {
        fprintf(stderr, "invalid file\n");
        exit(1);
    }
    kseq_t* seq(kseq_init(fq_fp));
    std::vector<uint64_t> locs;
    while ((err = kseq_read(seq)) >= 0) {
        rb_report(rbwt, midx, args, seq, locs);
    }
    // error checking here
    switch(err) {
        case -2:
            fprintf(stderr, "ERROR: truncated quality string\n");
            exit(1); break;
        case -3:
            fprintf(stderr, "ERROR: error reading stream\n");
            exit(1); break;
        default:
            break;
    }
}

int main(int argc, char** argv) {
    rb_locs_all(parse_args(argc, argv));
}

