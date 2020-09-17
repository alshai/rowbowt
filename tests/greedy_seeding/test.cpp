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

int main() {
    std::string fa("/home/taher/rowbowt/tests/greedy_seeding/ref.fa");
    std::string fq("/home/taher/rowbowt/tests/greedy_seeding/query.fq");
    rbwt::LoadRbwtFlag flag;
    flag = rbwt::LoadRbwtFlag::SA | rbwt::LoadRbwtFlag::DL;
    rbwt::RowBowt rbwt(rbwt::load_rowbowt(fa, flag));
    int err, nreads = 0, noccs = 0;
    gzFile fq_fp(gzopen(fq.data(), "r"));
    if (fq_fp == NULL) {
        fprintf(stderr, "invalid file\n");
        exit(1);
    }
    kseq_t* seq(kseq_init(fq_fp));
    std::vector<uint64_t> locs;
    while ((err = kseq_read(seq)) >= 0) {
        auto locs = rbwt.locate_pattern_greedy_seeding(seq->seq.s, 5, -1); // ret: {range, k}
        for (auto l: locs) {
            auto x = rbwt.resolve_offset(l);
            if (l != 10000) {
                std::cerr << "failed to locate: " << seq->name.s << std::endl;
                return 1;
            }
        }
        locs.clear();
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
    return 0;
}
