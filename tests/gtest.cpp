#include <cstdio>
#include <getopt.h>
#include <string>
#include <tuple>
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
            rbwt.locs_at(ret.rn, ret.ssamp, static_cast<uint64_t>(-1), locs);
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
       EXPECT_EQ(rets[0].rn, rbwt::RowBowt::range_t(24279,24280));
       EXPECT_EQ(rets[1].rn, rbwt::RowBowt::range_t(24175,24175));
       EXPECT_EQ(rets[2].rn, rbwt::RowBowt::range_t(27430,27432));
       EXPECT_EQ(rets[3].rn, rbwt::RowBowt::range_t(27430,27432));
       EXPECT_EQ(rets[4].rn, rbwt::RowBowt::range_t(17409,17409));
       EXPECT_EQ(rets[5].rn, rbwt::RowBowt::range_t(17416,17417));
   }

    rbwt::RowBowt rbwt;
    std::vector<rbwt::RowBowt::LFData> rets;
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
        lfs.clear();
        int i = 0;
        while ((err = kseq_read(seq)) >= 0) {
            auto lf = rbwt.find_range_w_markers(seq->seq.s, 10, -1);
            lfs.push_back(lf);
            ++i;
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
        EXPECT_EQ(get_pos(lfs[0].markers[0]), 289);
        EXPECT_EQ(get_allele(lfs[0].markers[0]), 0);
        EXPECT_EQ(get_pos(lfs[1].markers[0]), 289);
        EXPECT_EQ(get_allele(lfs[1].markers[0]), 1);
        EXPECT_EQ(lfs[2].markers.size(), 0); // second two reads do not have markers
        EXPECT_EQ(lfs[3].markers.size(), 0);
        EXPECT_EQ(get_pos(lfs[4].markers[0]), 4650);
        EXPECT_EQ(get_allele(lfs[4].markers[0]), 0);
        EXPECT_EQ(get_pos(lfs[5].markers[0]), 4650);
        EXPECT_EQ(get_allele(lfs[5].markers[0]), 1);
    }

    rbwt::RowBowt rbwt;
    std::vector<rbwt::RowBowt::LFData> lfs;
};

class GreedySeedingTester : public testing::Test {
    protected:
    void SetUp() override {
        std::string fa("/home/taher/rowbowt/tests/data/small.fa");
        std::string fq("/home/taher/rowbowt/tests/data/error_query.fq");
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
        int i = 0;
        while ((err = kseq_read(seq)) >= 0) {
            auto ret = rbwt.get_seeds_greedy_w_sample(seq->seq.s, 10);
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
        std::vector<std::vector<uint64_t>> all_locs;
        for (auto r: rets) {
            auto locs = rbwt.locate_from_longest_seed(static_cast<uint64_t>(-1), r);
            all_locs.push_back(locs);
        }
        EXPECT_EQ(all_locs[0][0], 10296);
        EXPECT_EQ(all_locs[0][1], 20306);
        EXPECT_EQ(all_locs[0][2], 286);
        EXPECT_EQ(all_locs[1][0], 10296);
        EXPECT_EQ(all_locs[2][0], 11897);
        EXPECT_EQ(all_locs[2][1], 21907);
        EXPECT_EQ(all_locs[2][2], 1887);
        EXPECT_EQ(all_locs[3][0], 11897);
        EXPECT_EQ(all_locs[3][1], 21907);
        EXPECT_EQ(all_locs[3][2], 1887);
        EXPECT_EQ(all_locs[4].size(), 0);
        EXPECT_EQ(all_locs[5][0], 14654);
        EXPECT_EQ(all_locs[5][1], 4644);
    }

   void CountTester() {
       // EXPECT_EQ(rets[0].first, rbwt::RowBowt::range_t(24279,24280));
       // EXPECT_EQ(rets[1].first, rbwt::RowBowt::range_t(24175,24175));
       // EXPECT_EQ(rets[2].first, rbwt::RowBowt::range_t(27430,27432));
       // EXPECT_EQ(rets[3].first, rbwt::RowBowt::range_t(27430,27432));
       // EXPECT_EQ(rets[4].first, rbwt::RowBowt::range_t(17409,17409));
       // EXPECT_EQ(rets[5].first, rbwt::RowBowt::range_t(17416,17417));
   }

    rbwt::RowBowt rbwt;
    std::vector<std::vector<rbwt::RowBowt::LFData>> rets;
};

class GreedySeedingMarkerTester : public testing::Test {
    protected:
    void SetUp() override {
        std::string fa("/home/taher/rowbowt/tests/data/small.fa");
        std::string fq("/home/taher/rowbowt/tests/data/error_query.fq");
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
        int i = 0;
        while ((err = kseq_read(seq)) >= 0) {
            std::vector<MarkerT> ret;
            ret.clear();
            auto fn = [&](rbwt::RowBowt::range_t p, std::pair<size_t, size_t> q, std::vector<MarkerT> mbuf) {
                std::cerr << " [" << p.second - p.first + 1 << "]" << ":[" << q.first << "-" << q.second << "]";
                if (mbuf.size()) {
                    for (auto m: mbuf) {
                        ret.push_back(m);
                    }
                } 
            };
            rbwt.get_markers_greedy_seeding(seq->seq.s, 10, 100, fn);
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

    void MarkerTester() {
        std::cerr << "todo: properly generate truth set for this test\n" << std::endl;
    }

    rbwt::RowBowt rbwt;
    std::vector<std::vector<MarkerT>> rets;
};

class LMemMarkerTester : public testing::Test {
    protected:
    void SetUp() override {
        std::string fa("/home/taher/rowbowt/tests/data/small.fa");
        std::string fq("/home/taher/rowbowt/tests/data/error_query.fq");
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
        int i = 0;
        while ((err = kseq_read(seq)) >= 0) {
            std::vector<MarkerT> ret;
            ret.clear();
            std::cerr << seq->name.s;
            auto fn = [&](rbwt::RowBowt::range_t p, std::pair<size_t, size_t> q, std::vector<MarkerT> mbuf) {
                std::cerr << " [" << p.second - p.first + 1 << "]" << ":[" << q.first << "-" << q.second << "]";
                if (mbuf.size()) {
                    for (auto m: mbuf) {
                        ret.push_back(m);
                    }
                }
            };
            rbwt.get_markers_lmems(seq->seq.s, 10, 100, fn);
            std::cerr << std::endl;
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

    void MarkerTester() {
        std::cerr << "todo: properly generate truth set for this test\n" << std::endl;
    }

    rbwt::RowBowt rbwt;
    std::vector<std::vector<MarkerT>> rets;
};

class FTabTester  : public testing::Test {
    protected:
    void SetUp() override {
        std::string fa("/home/taher/rowbowt/tests/data/small.fa");
        rbwt::LoadRbwtFlag flag;
        flag = rbwt::LoadRbwtFlag::FT;
        rbwt = rbwt::load_rowbowt(fa, flag);

    }

    void FTabTest10mer() {
        std::pair<size_t, size_t> rn;
        rn = rbwt.find_range("TTCGTCGTAA");
        EXPECT_EQ(rn.first, 28942);
        EXPECT_EQ( rn.second, 28944);
        rn = rbwt.find_range("CCGCGGACAT");
        EXPECT_EQ(rn.first, 10673);
        EXPECT_EQ( rn.second, 10675);
        rn = rbwt.find_range("GGCAGGCGGA");
        EXPECT_EQ(rn.first,  19418);
        EXPECT_EQ( rn.second, 19423);
    }

    void FTabTest11merfrom10merFtab() {
        std::pair<size_t, size_t> rn;
        rn = rbwt.find_range("TATCGTGGAA");
        EXPECT_EQ(rn.first,  24272);
        EXPECT_EQ( rn.second, 24274);
        rn = rbwt.find_range("GTATCGTGGAA");
        EXPECT_EQ(rn.first,  21142);
        EXPECT_EQ( rn.second, 21144);
        rn = rbwt.find_range("GGAGATATTG");
        EXPECT_EQ(rn.first,  19097);
        EXPECT_EQ( rn.second, 19099);
        rn = rbwt.find_range("TGGAGATATTG");
        EXPECT_EQ(rn.first,  27180);
        EXPECT_EQ( rn.second, 27182);
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

TEST_F(GreedySeedingTester, Locate) {
    LocateTester();
}

TEST_F(GreedySeedingMarkerTester, Marker) {
    MarkerTester();
}

TEST_F(LMemMarkerTester, Marker) {
    MarkerTester();
}

TEST_F(FTabTester, FTab) {
    FTabTest10mer();
}

TEST_F(FTabTester, FTabExtend) {
    FTabTest11merfrom10merFtab();
}