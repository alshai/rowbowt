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
#include "gtest/gtest.h"

namespace {

// Tests simple seeding scheme which assumes no errors in reads
// ie. does full LF mapping of the read from rightmost end to leftmost end.
class SimpleSeedingTester : public testing::Test {
    protected:
    void SetUp() override {
        std::string fa("/home/taher/rowbowt/tests/data/small.fa");
        std::string fq("/home/taher/rowbowt/tests/data/simple_query.fq");
        rbwt::LoadRbwtFlag flag;
        flag = rbwt::LoadRbwtFlag::SA | rbwt::LoadRbwtFlag::DL;
        rbwt = rbwt::load_rowbowt(fa, flag);
        int err, nreads = 0, noccs = 0;
        gzFile fq_fp(gzopen(fq.data(), "r"));
        if (fq_fp == NULL) {
            fprintf(stderr, "invalid file\n");
            exit(1);
        }
        kseq_t* seq(kseq_init(fq_fp));
        rets.clear();
        while ((err = kseq_read(seq)) >= 0) {
            auto ret = rbwt.find_range_w_toehold(seq->seq.s);
            rets.push_back(ret);
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

    void LocateTester() {
        std::vector<uint64_t> locs;
        std::vector<uint64_t> all_locs;
        for (auto ret: rets) {
            locs.clear();
            rbwt.locs_at(ret.first, ret.second, static_cast<uint64_t>(-1), locs);
            all_locs.insert(all_locs.end(), locs.begin(), locs.end());
        }
        // EXPECT_EQ(all_locs.size(), 6);
        EXPECT_EQ(all_locs[0], 20306);
        EXPECT_EQ(all_locs[1], 286);
        EXPECT_EQ(all_locs[2], 10296);
        EXPECT_EQ(all_locs[3], 11897);
        EXPECT_EQ(all_locs[4], 21907);
        EXPECT_EQ(all_locs[5], 1887);
        EXPECT_EQ(all_locs[6], 11897);
        EXPECT_EQ(all_locs[7], 21907);
        EXPECT_EQ(all_locs[8], 1887);
        EXPECT_EQ(all_locs[9], 4644);
        EXPECT_EQ(all_locs[10], 14654);
        EXPECT_EQ(all_locs[11], 24664);
    }

   void CountTester() {
       EXPECT_EQ(rets[0].first, rbwt::RowBowt::range_t(24279,24280));
       EXPECT_EQ(rets[1].first, rbwt::RowBowt::range_t(24175,24175));
       EXPECT_EQ(rets[2].first, rbwt::RowBowt::range_t(27430,27432));
       EXPECT_EQ(rets[3].first, rbwt::RowBowt::range_t(27430,27432));
       EXPECT_EQ(rets[4].first, rbwt::RowBowt::range_t(17409,17409));
       EXPECT_EQ(rets[5].first, rbwt::RowBowt::range_t(17416,17417));
   }

    rbwt::RowBowt rbwt;
    std::vector<std::pair<rbwt::RowBowt::range_t, uint64_t>> rets;
};
//
// Tests simple seeding scheme which assumes no errors in reads
// ie. does full LF mapping of the read from rightmost end to leftmost end.
class SimpleMarkerTester : public testing::Test {
    protected:
    void SetUp() override {
        std::string fa("/home/taher/rowbowt/tests/data/small.fa");
        std::string fq("/home/taher/rowbowt/tests/data/simple_query.fq");
        rbwt::LoadRbwtFlag flag;
        flag = rbwt::LoadRbwtFlag::MA;
        rbwt = rbwt::load_rowbowt(fa, flag);
        int err, nreads = 0, noccs = 0;
        gzFile fq_fp(gzopen(fq.data(), "r"));
        if (fq_fp == NULL) {
            fprintf(stderr, "invalid file\n");
            exit(1);
        }
        kseq_t* seq(kseq_init(fq_fp));
        rets.clear();
        while ((err = kseq_read(seq)) >= 0) {
            auto markers = rbwt.get_markers_simple_seeding(seq->seq.s, 4);
            rets.push_back(markers);
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

    void MarkerTester() {
        EXPECT_EQ(get_pos(rets[0][0]), 289);
        EXPECT_EQ(get_allele(rets[0][0]), 0);
        EXPECT_EQ(get_pos(rets[1][0]), 289);
        EXPECT_EQ(get_allele(rets[1][0]), 1);
        EXPECT_EQ(rets[2].size(), 0); // second two reads do not have markers
        EXPECT_EQ(rets[3].size(), 0);
        EXPECT_EQ(get_pos(rets[4][0]), 4650);
        EXPECT_EQ(get_allele(rets[4][0]), 0);
        EXPECT_EQ(get_pos(rets[5][0]), 4650);
        EXPECT_EQ(get_allele(rets[5][0]), 1);
    }

    rbwt::RowBowt rbwt;
    std::vector<std::vector<MarkerT>> rets;
};

}

TEST_F(SimpleSeedingTester, Count) {
    CountTester();
}

TEST_F(SimpleSeedingTester, Locate) {
    LocateTester();
}

TEST_F(SimpleMarkerTester, Marker) {
    MarkerTester();
}
