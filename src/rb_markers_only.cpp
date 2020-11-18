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
    int inexact = 0;
    size_t wsize = 10;
    size_t max_range = 1000;
};

void print_help() {
    fprintf(stderr, "rb_markers_only");
    fprintf(stderr, "Usage: rb_markers_only [options] <index_prefix> <input_fastq_name>\n");
    fprintf(stderr, "    --output_prefix/-o <basename>    output prefix\n");
    fprintf(stderr, "    --inexact                                     \n");
    fprintf(stderr, "    --wsize            <int>         window size for performing marker queries along read\n");
    fprintf(stderr, "    --max-range        <int>         range-size upper threshold for performing marker queries\n");
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
        {"inexact", no_argument, &args.inexact, 1},
        {"max-range", required_argument, 0, 'r'}
    };
    int long_index = 0;
    while((c = getopt_long(argc, argv, "o:w:r:h", long_options, &long_index)) != -1) {
        switch (c) {
            case 'o':
                args.outpre = optarg;
                break;
            case 'r':
                args.max_range = std::atol(optarg);
                break;
            case 'w':
                args.wsize = std::atol(optarg);
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


struct RowBowtRet {
    rbwt::RowBowt::range_t r;
    uint64_t ts; // toehold sample
    std::vector<uint64_t> locs;
    std::vector<MarkerT> markers;
};

void rb_report(const rbwt::RowBowt& rbwt, const RbAlignArgs args, kseq_t* seq, std::vector<MarkerT>& markers) {
    markers.clear();
    auto fn = [&](rbwt::RowBowt::range_t p, std::pair<size_t, size_t> q, std::vector<MarkerT> mbuf) {
        std::cerr << " [" << p.second - p.first + 1 << "]" << ":[" << q.first << "-" << q.second << "]";
        if (mbuf.size()) {
            for (auto m: mbuf) {
                std::cerr << ":" << get_pos(m) << "/" << static_cast<int>(get_allele(m));
                markers.push_back(m);
            }
        } 
    };
    std::cout << seq->name.s;
    std::cerr << seq->name.s;
    rbwt.get_markers_greedy_seeding(seq->seq.s, args.wsize, args.max_range, fn);
    std::sort(markers.begin(), markers.end(), [](MarkerT a, MarkerT b) {
            return get_pos(a) < get_pos(b);
    });
    markers.erase(std::unique(markers.begin(), markers.end()), markers.end());
    for (auto m: markers) {
        std::cout << " " << get_pos(m) << "/" << static_cast<int>(get_allele(m));
    } 
    std::cerr << "\n";
    std::cout << "\n";
}

rbwt::RowBowt load_rbwt(const RbAlignArgs args) {
    rbwt::LoadRbwtFlag flag;
    flag = flag | rbwt::LoadRbwtFlag::MA;
    rbwt::RowBowt rbwt(rbwt::load_rowbowt(args.inpre, flag));
    return rbwt;
}

void rb_markers_all(RbAlignArgs args) {
    rbwt::RowBowt rbwt(load_rbwt(args));
    int err, nreads = 0, noccs = 0;
    gzFile fq_fp(gzopen(args.fastq_fname.data(), "r"));
    if (fq_fp == NULL) {
        fprintf(stderr, "invalid file\n");
        exit(1);
    }
    kseq_t* seq(kseq_init(fq_fp));
    std::vector<MarkerT> markers;
    while ((err = kseq_read(seq)) >= 0) {
        if (!args.inexact) {
            rb_report(rbwt, args, seq, markers);
        } else {
            rb_report(rbwt, args, seq, markers);
        }
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
    rb_markers_all(parse_args(argc, argv));
}

