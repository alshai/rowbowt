#include <cstdio>
#include <getopt.h>
#include <string>
#include "rowbowt_io.hpp"
#include "rle_string.hpp"
#include "fbb_string.hpp"

/* build rl-bwt from input bwt and saves to disk */


void print_help() {
    fprintf(stderr, "rb_build\n");
    fprintf(stderr, "Usage: rb_build [options] <index_prefix>\n");
    fprintf(stderr, "    --output_prefix/-o <basename>    output prefix\n");
    fprintf(stderr, "    --tsa/-s <basename>                 build toehold suffix array\n");
    fprintf(stderr, "    --ma/-m <basename>                  build marker array\n");
    fprintf(stderr, "    --ftab/-f <basename>                construct offset table (for faster querying)\n");
    fprintf(stderr, "    --fbb                               use fbb_string instead of rle_string\n");
    fprintf(stderr, "    <input_prefix>                   index prefix\n");
}

rbwt::RowBowtConstructArgs parse_args(int argc, char** argv) {
    int c;
    rbwt::RowBowtConstructArgs args;
    static struct option long_options[] {
        {"output-prefix", required_argument, 0, 'o'},
        {"tsa", no_argument, 0, 's'},
        {"dl", no_argument, 0, 'l'},
        {"ftab-only", no_argument, 0, 'a'},
        {"ma", no_argument, 0, 'm'},
        {"ft", no_argument, 0, 'f'},
        {"fbb", no_argument, 0, 'x'},
        {0, 0, 0, 0}
    };
    int long_index = 0;
    while((c = getopt_long(argc, argv, "xo:k:lfsmha", long_options, &long_index)) != -1) {
        switch (c) {
            case 'x':
                args.fbb = 1;
                break;
            case 'o':
                args.prefix = optarg;
                break;
            case 's':
                args.tsa = 1;
                break;
            case 'm':
                args.ma = 1;
                break;
            case 'l':
                args.dl = 1;
                break;
            case 'f':
                args.ft = 1;
                break;
            case 'a':
                args.ft_only = 1;
                break;
            case 'k':
                args.k = std::stoull(optarg); break;
            case 'h':
                print_help();
                exit(0);
                break;
            case '?':
                break;
            default:
                print_help();
                exit(1);
                break;
        }
    }

    if (argc - optind < 1) {
        fprintf(stderr, "no argument provided\n");
        exit(1);
    }

    std::string inpre = argv[optind++];
    if (args.prefix == "") {
        args.prefix = inpre;
    }
    args.bwt_fname = inpre + ".bwt";
    if (args.tsa) {
        args.ssa_fname = inpre + ".ssa";
        args.esa_fname  = inpre + ".esa";
    }
    if (args.ma) {
        args.ma_fname = inpre + ".ma";
    }
    if (args.dl) {
        args.dl_fname = inpre + ".docs";
    }
    return args;
}

void rb_ftab(rbwt::RowBowtConstructArgs args) {
    rbwt::construct_and_serialize_ftab(args);
}

template<typename T=ri::rle_string_sd>
void rb_build(rbwt::RowBowtConstructArgs args) {
    rbwt::construct_and_serialize_rowbowt<T>(args);
}

int main(int argc, char** argv) {
    auto args = parse_args(argc, argv);
    if (args.ft_only) {
        rb_ftab(args);
    } else {
        if (args.fbb) {
            rb_build<ri::fbb_string>(args);
        } else {
            rb_build<ri::rle_string_sd>(args);
        }
    }
    return 0;
}
