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
    size_t min_range = 0;
};

void print_help() {
    fprintf(stderr, "rb_markers_only");
    fprintf(stderr, "Usage: rb_markers_only [options] <index_prefix> <input_fastq_name>\n");
    fprintf(stderr, "    --output_prefix/-o <basename>    output prefix\n");
    fprintf(stderr, "    --inexact                                     \n");
    fprintf(stderr, "    --wsize            <int>         window size for performing marker queries along read\n");
    fprintf(stderr, "    --max-range        <int>         range-size upper threshold for performing marker queries\n");
    fprintf(stderr, "    --min-range        <int>         range-size upper threshold for performing marker queries\n");
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
        {"max-range", required_argument, 0, 'r'},
        {"min-range", required_argument, 0, 'm'}
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
            case 'm':
                args.min_range = std::atol(optarg);
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

char comp_tab[] = {
    0,    1,    2,    3,    4,    5,    6,    7,    8,    9,    10,   11,   12,   13,   14,   15,
    16,   17,   18,   19,   20,   21,   22,   23,   24,   25,   26,   27,   28,   29,   30,   31,
    32,   33,   34,   35,   36,   37,   38,   39,   40,   41,   42,   43,   44,   45,   46,   47,
    48,   49,   50,   51,   52,   53,   54,   55,   56,   57,   58,   59,   60,   61,   62,   63,
    64,   'T',  'V',  'G',  'H',  'E',  'F',  'C',  'D',  'I',  'J',  'M',  'L',  'K',  'N',  'O',
    'P',  'Q',  'Y',  'S',  'A',  'A',  'B',  'W',  'X',  'R',  'Z',  91,   92,   93,   94,   95,
    64,   't',  'v',  'g',  'h',  'e',  'f',  'c',  'd',  'i',  'j',  'm',  'l',  'k',  'n',  'o',
    'p',  'q',  'y',  's',  'a',  'a',  'b',  'w',  'x',  'r',  'z',  123,  124,  125,  126,  127
};

/* code copied from seqtk
 * https://github.com/lh3/seqtk/blob/7c04ce7898ad5909bd309c6ba3cd9c3bd0651f0e/seqtk.c#L1464
 */
void revc_in_place(kseq_t* seq) {
    int c0, c1;
    for (size_t i = 0; i < seq->seq.l>>1; ++i) { // reverse complement sequence
        c0 = comp_tab[(int)seq->seq.s[i]];
        c1 = comp_tab[(int)seq->seq.s[seq->seq.l - 1 - i]];
        seq->seq.s[i] = c1;
        seq->seq.s[seq->seq.l - 1 - i] = c0;
    }
    if (seq->seq.l & 1) // complement the remaining base
        seq->seq.s[seq->seq.l>>1] = comp_tab[(int)seq->seq.s[seq->seq.l>>1]];
    if (seq->qual.l) {
        for (size_t i = 0; i < seq->seq.l>>1; ++i) // reverse quality
            c0 = seq->qual.s[i], seq->qual.s[i] = seq->qual.s[seq->qual.l - 1 - i], seq->qual.s[seq->qual.l - 1 - i] = c0;
    }
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
        if (p.second - p.first + 1 >= args.min_range && mbuf.size()) {
            for (auto m: mbuf) {
                std::cerr << ":" << get_seq(m) << "/" <<  get_pos(m) << "/" << static_cast<int>(get_allele(m));
                markers.push_back(m);
            }
        }
    };
    std::cout << seq->name.s;
    std::cerr << seq->name.s;
    rbwt.get_markers_greedy_seeding(seq->seq.s, args.wsize, args.max_range, fn);
    revc_in_place(seq);
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
    auto flag = rbwt::LoadRbwtFlag::MA;
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
        rb_report(rbwt, args, seq, markers);
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

