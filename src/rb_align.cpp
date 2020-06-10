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
    int sam = 0;
    int markers = 0;
};

void print_help() {
    fprintf(stderr, "rb_align");
    fprintf(stderr, "Usage: rb_align [options] <index_prefix> <input_fastq_name>\n");
    fprintf(stderr, "    --output_prefix/-o <basename>    output prefix\n");
    fprintf(stderr, "    --markers/-m <basename>          print markers to <output_prefix>.markers\n");
    fprintf(stderr, "    --sam/-s <basename>              print locations to <output_prefix>.sam\n");
    fprintf(stderr, "    <input_prefix>                   index prefix\n");
    fprintf(stderr, "    <input_fastq>                    input fastq\n");
}

RbAlignArgs parse_args(int argc, char** argv) {
    int c;
    char* end;
    RbAlignArgs args;
    static struct option long_options[] {
        {"output_prefix", required_argument, 0, 'o'},
        {"markers", required_argument, 0, 'm'}, // outputs markers
        {"sam", required_argument, 0, 's'} // output locations in SAM format
    };
    int long_index = 0;
    while((c = getopt_long(argc, argv, "o:smh", long_options, &long_index)) != -1) {
        switch (c) {
            case 'o':
                args.outpre = optarg;
                break;
            case 'h':
                print_help();
                exit(0);
                break;
            case 's':
                args.sam = 1;
                break;
            case 'm':
                args.markers = 1;
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


struct RowBowtRet {
    rbwt::RowBowt::range_t r;
    uint64_t ts; // toehold sample
    std::vector<uint64_t> locs;
    std::vector<MarkerT> markers;
};

RowBowtRet rb_get_range(const rbwt::RowBowt& rbwt, std::string query, bool sa) {
    RowBowtRet ret;
    if (sa)  {
        auto range = rbwt.find_range_w_toehold(query);
        ret.r = range.first;
        ret.ts = range.second;
    }
    else {
        ret.r = rbwt.find_range(query);
    }
    return ret;
}

void rb_report(const rbwt::RowBowt& rbwt, const RbAlignArgs args, kseq_t* seq) {
    std::cout << seq->name.s << " ";
    RowBowtRet ret = rb_get_range(rbwt, seq->seq.s, args.sam); // TODO: other options could also potentially trigger SA
    std::cout << "(" << ret.r.first << "," << ret.r.second << "):\n";
    if (args.sam) {
        // TODO: output SAM line to args.output + SAM
        rbwt.locs_at(ret.r, ret.ts, static_cast<uint64_t>(-1), ret.locs);
        std:: cout << "\tlocs: ";
        for (auto l: ret.locs) {
            std:: cout << l << " ";
        } std::cout << "\n";
    }
    if (args.markers) {
        rbwt.markers_at(ret.r, ret.markers);
        std:: cout << "\tmarkers: ";
        if (!ret.markers.size()) cout << "no markers (consider building the marker array with a larger window size)";
        for (auto m:  ret.markers) {
            std:: cout << get_pos(m) << "/" << static_cast<int>(get_allele(m)) << " ";
        } std::cout << "\n";
    }
}

rbwt::RowBowt load_rbwt(const RbAlignArgs args) {
    rbwt::LoadRbwtFlag flag;
    if (args.sam) flag = flag | rbwt::LoadRbwtFlag::SA;
    if (args.markers) flag = flag | rbwt::LoadRbwtFlag::MA;
    rbwt::RowBowt rbwt(rbwt::load_rowbowt(args.inpre, flag));
    return rbwt;
}

void rb_align_all(RbAlignArgs args) {
    rbwt::RowBowt rbwt(load_rbwt(args));
    int err, nreads = 0, noccs = 0;
    gzFile fq_fp(gzopen(args.fastq_fname.data(), "r"));
    if (fq_fp == NULL) {
        fprintf(stderr, "invalid file\n");
        exit(1);
    }
    kseq_t* seq(kseq_init(fq_fp));
    while ((err = kseq_read(seq)) >= 0) {
        rb_report(rbwt, args, seq);
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
    rb_align_all(parse_args(argc, argv));
}
