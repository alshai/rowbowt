#include <cstdio>
#include <getopt.h>
#include <string>
#include <chrono>
extern "C" {
#include <zlib.h>
#ifndef AC_KSEQ_H
#include "kseq.h"
KSEQ_INIT(gzFile, gzread);
#endif
}
#include "rowbowt.hpp"
#include "rowbowt_io.hpp"
#include "fbb_string.hpp"
#include "rle_string.hpp"

struct RbAlignArgs {
    std::string inpre = "";
    std::string fastq_fname = "";
    std::string outpre = "";
    int sam = 0;
    int markers = 0;
    int fbb = 0;
};

void print_help() {
    fprintf(stderr, "rb_align");
    fprintf(stderr, "Usage: rb_align [options] <index_prefix> <input_fastq_name>\n");
    fprintf(stderr, "    --output_prefix/-o <basename>    output prefix\n");
    fprintf(stderr, "    --markers/-m <basename>          print markers to <output_prefix>.markers\n");
    fprintf(stderr, "    --sam/-s <basename>              print locations to <output_prefix>.sam\n");
    fprintf(stderr, "    --fbb                            index is based on wt-fbb\n");
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
        {"sam", required_argument, 0, 's'}, // output locations in SAM format
        {"fbb", no_argument, 0, 'f'} // output locations in SAM format
    };
    int long_index = 0;
    while((c = getopt_long(argc, argv, "o:smh", long_options, &long_index)) != -1) {
        switch (c) {
            case 'f':
                args.fbb = 1;
                break;
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


template<typename StringT>
struct RowBowtRet {
    typename rbwt::RowBowt<StringT>::range_t r;
    uint64_t ts; // toehold sample
    std::vector<uint64_t> locs;
    std::vector<MarkerT> markers;
};

template<typename StringT>
typename std::enable_if<std::is_same<StringT,ri::rle_string_sd>::value, RowBowtRet<StringT> >::type
rb_get_range(const rbwt::RowBowt<StringT>& rbwt, std::string query, bool sa) {
    RowBowtRet<StringT> ret;
    if (sa)  {
        auto range = rbwt.find_range_w_toehold(query);
        ret.r = range.rn;
        ret.ts = range.ssamp;
    }
    else {
        ret.r = rbwt.find_range(query);
    }
    return ret;
}

template<typename StringT>
typename std::enable_if<std::is_same<StringT,ri::fbb_string>::value, RowBowtRet<StringT> >::type
rb_get_range(const rbwt::RowBowt<StringT>& rbwt, std::string query, bool sa) {
    RowBowtRet<StringT> ret;
    ret.r = rbwt.find_range(query);
    return ret;
}

template<typename StringT>
void rb_report(const rbwt::RowBowt<StringT>& rbwt, const RbAlignArgs args, kseq_t* seq) {
    std::cout << seq->name.s << " ";
    RowBowtRet<StringT> ret = rb_get_range<StringT>(rbwt, seq->seq.s, args.sam); // TODO: other options could also potentially trigger SA
    std::cout << "(" << ret.r.first << "," << ret.r.second << "), count=" << ret.r.second-ret.r.first + 1 << "\n";
    if (args.sam) {
        // TODO: output SAM line to args.output + SAM
        rbwt.locs_at(ret.r, ret.ts, static_cast<uint64_t>(-1), ret.locs);
        std:: cout << "\tlocs: ";
        for (auto l: ret.locs) {
            std::cout << l << "/";
            auto x = rbwt.resolve_offset(l);
            std::cout << x.first << ":" << x.second << " ";
        } std::cout << "\n";
    }
    if (args.markers) {
        rbwt.markers_at(ret.r.first, ret.markers);
        ret.markers.clear();
        rbwt.markers_at(ret.r.second, ret.markers);
        ret.markers.clear();
        rbwt.markers_at(ret.r, ret.markers);
        std:: cout << "\tmarkers: ";
        if (!ret.markers.size()) cout << "no markers (consider building the marker array with a larger window size)";
        for (auto m:  ret.markers) {
            std:: cout << get_pos(m) << "/" << static_cast<int>(get_allele(m)) << " ";
        } std::cout << "\n";
    }
}

template<typename StringT>
rbwt::RowBowt<StringT> load_rbwt(const RbAlignArgs args) {
    rbwt::LoadRbwtFlag flag = static_cast<rbwt::LoadRbwtFlag>(0);
    if (args.sam) {
        std::cerr << "will load SA and DA" << std::endl;
        flag = flag | rbwt::LoadRbwtFlag::SA | rbwt::LoadRbwtFlag::DL;
    }
    if (args.markers) {
        std::cerr << "will load SA and DA" << std::endl;
        flag = flag | rbwt::LoadRbwtFlag::MA;
    }
    rbwt::RowBowt<StringT> rbwt(rbwt::load_rowbowt<StringT>(args.inpre, flag));
    return rbwt;
}

template<typename StringT>
void rb_align_all(RbAlignArgs args) {
    auto start = std::chrono::high_resolution_clock::now();
    rbwt::RowBowt<StringT> rbwt(load_rbwt<StringT>(args));
    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> index_load_time = stop - start;
    int err, nreads = 0, noccs = 0;
    gzFile fq_fp(gzopen(args.fastq_fname.data(), "r"));
    if (fq_fp == NULL) {
        fprintf(stderr, "invalid file\n");
        exit(1);
    }
    kseq_t* seq(kseq_init(fq_fp));
    start = std::chrono::high_resolution_clock::now();
    while ((err = kseq_read(seq)) >= 0) {
        rb_report<StringT>(rbwt, args, seq);
    }
    stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> total_query_time = stop - start;
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
    std::cerr << index_load_time.count() << " " << total_query_time.count() << std::endl;
}

int main(int argc, char** argv) {
    auto args = parse_args(argc, argv);
    if (args.fbb) {
        rb_align_all<ri::fbb_string>(args);
    } else {
        rb_align_all<ri::rle_string_sd>(args);
    }
}
