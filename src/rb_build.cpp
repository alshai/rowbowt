#include <cstdio>
#include <getopt.h>
#include <string>
#include "rowbowt_io.hpp"

/* build rl-bwt from input bwt and saves to disk */


void print_help() {
    fprintf(stderr, "rb_build");
    fprintf(stderr, "Usage: rb_build [options] <index_prefix>\n");
    fprintf(stderr, "    --output_prefix/-o <basename>    output prefix\n");
    fprintf(stderr, "    --tsa <basename>                 build toehold suffix array\n");
    fprintf(stderr, "    --ma <basename>                  build marker array\n");
    fprintf(stderr, "    <input_prefix>                   index prefix\n");
}

rbwt::RowBowtConstructArgs parse_args(int argc, char** argv) {
    int c;
    rbwt::RowBowtConstructArgs args;
    static struct option long_options[] {
        {"output_prefix", required_argument, 0, 'o'},
        {"tsa", required_argument, 0, 's'},
        {"ma", required_argument, 0, 'm'},
        {"dl", required_argument, 0, 'l'},
    };
    int long_index = 0;
    while((c = getopt_long(argc, argv, "o:lsmh", long_options, &long_index)) != -1) {
        switch (c) {
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

void rb_build(rbwt::RowBowtConstructArgs args) {
    rbwt::construct_and_serialize_rowbowt(args);
}

int main(int argc, char** argv) {
    rb_build(parse_args(argc, argv));
    return 0;
}

