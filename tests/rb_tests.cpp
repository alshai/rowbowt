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
#include "fbb_string.hpp"
#include "gtest/gtest.h"

// Tests simple seeding scheme which assumes no errors in reads
// ie. does full LF mapping of the read from rightmost end to leftmost end.
template<typename StringT>
class RowBowtTester : public testing::Test {
protected:
    using LFData = typename rbwt::RowBowt<StringT>::LFData;
    using range_t = typename rbwt::RowBowt<StringT>::range_t;

    void SetUp() override {
        rbwt = make_rowbowt();
    }

    template<typename T=StringT>
    typename std::enable_if<std::is_same<T,ri::rle_string_sd>::value, void>::type
    LocateTester() {
        std::vector<LFData> rets;
        std::string fq("/home/taher/rowbowt/tests/data/simple_query.fq");
        auto fn = [&](kseq_t* seq) {
            auto ret = rbwt.find_range_w_toehold(seq->seq.s);
            rets.push_back(ret);
        };
        kseq_for_each(fq, fn);
        std::vector<uint64_t> locs;
        std::vector<uint64_t> all_locs;
        for (auto ret : rets) {
            locs.clear();
            rbwt.locs_at(ret.rn, ret.ssamp, static_cast<uint64_t>(-1), locs);
            all_locs.insert(all_locs.end(), locs.begin(), locs.end());
        }
        // EXPECT_EQ(all_locs.size(), 12);
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

    template<typename T=StringT>
    typename std::enable_if<std::is_same<T,ri::fbb_string>::value, void>::type
    LocateTester() {
        EXPECT_EQ(0,0);
    }

    template<typename T=StringT>
    typename std::enable_if<std::is_same<T,ri::rle_string_sd>::value, void>::type
    GreedyLocateTester() {
        std::string fq("/home/taher/rowbowt/tests/data/error_query.fq");
        std::vector<LFData> rets;
        auto ret_fn = [&](kseq_t* seq) {
            auto ret = rbwt.get_seeds_greedy_w_sample(seq->seq.s, 10);
            rets.push_back(ret);
        };
        kseq_for_each(fq, ret_fn);

        std::vector<std::vector<uint64_t>> all_locs;
        for (auto r : rets) {
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

    template<typename T=StringT>
    typename std::enable_if<std::is_same<T,ri::fbb_string>::value, void>::type
    GreedyLocateTester() {
        EXPECT_EQ(0,0);
    }

    void CountTester() {
        std::string fq("/home/taher/rowbowt/tests/data/simple_query.fq");
        std::vector<range_t> rets;
        auto fn = [&](kseq_t* seq) {
            auto ret = rbwt.find_range(seq->seq.s);
            rets.push_back(ret);
        };
        kseq_for_each(fq, fn);
        for (auto ret : rets) {
            rets.push_back(ret);
        }
        EXPECT_EQ(rets[0], range_t(24279, 24280));
        // EXPECT_EQ(rets[1], range_t(24175, 24175));
        // EXPECT_EQ(rets[2], range_t(27430, 27432));
        // EXPECT_EQ(rets[3], range_t(27430, 27432));
        // EXPECT_EQ(rets[4], range_t(17409, 17409));
        // EXPECT_EQ(rets[5], range_t(17416, 17417));
    }

    void MarkerTester() {
        typename std::vector<LFData> lfs;
        auto fn = [&](kseq_t* seq) {
            auto lf = rbwt.find_range_w_markers(seq->seq.s, 10, -1);
            lfs.push_back(lf);
        };
        std::string fq("/home/taher/rowbowt/tests/data/simple_query.fq");
        kseq_for_each(fq, fn);
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

    void FTabTest10mer() {
        FTab ftab = rbwt.build_ftab(10);
        rbwt.set_ftab(std::move(ftab));
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
        FTab ftab = rbwt.build_ftab(10);
        rbwt.set_ftab(std::move(ftab));
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



    private: 

    rbwt::RowBowt<StringT> make_rowbowt() {
        StringT bwt(fa + ".bwt");
        MarkerArray<> ma(fa + ".ma");
        ToeholdSA tsa = load_tsa(bwt);
        rbwt::RowBowt<StringT> rb(std::move(bwt), std::move(ma), std::move(tsa), {}, {});
        return rb;
    }

    template<typename T=StringT>
    typename std::enable_if<std::is_same<T,ri::rle_string_sd>::value, ToeholdSA>::type
    load_tsa(StringT& bwt) {
        return ToeholdSA(bwt.size(), bwt.number_of_runs(), fa + ".ssa", fa + ".esa");
    }

    template<typename T=StringT>
    typename std::enable_if<std::is_same<T,ri::fbb_string>::value, ToeholdSA>::type
    load_tsa(StringT& bwt) {
        return ToeholdSA();
    }



    template<typename Fn>
    void kseq_for_each(std::string fq, Fn f) {
        int err, nreads = 0, noccs = 0;
        gzFile fq_fp(gzopen(fq.data(), "r"));
        if (fq_fp == NULL) {
            fprintf(stderr, "invalid file\n");
            exit(1);
        }
        kseq_t* seq(kseq_init(fq_fp));
        int i = 0;
        while ((err = kseq_read(seq)) >= 0) {
            f(seq);
        }
        switch (err) {
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

    rbwt::RowBowt<StringT> rbwt;
    std::string fa = "/home/taher/rowbowt/tests/data/small.fa";
};


/*
template<typename StringT>
class LMemMarkerTester : public testing::Test {
    protected:
    void SetUp() override {
        std::string fa("/home/taher/rowbowt/tests/data/small.fa");
        std::string fq("/home/taher/rowbowt/tests/data/error_query.fq");
        rbwt::LoadRbwtFlag flag;
        flag = rbwt::LoadRbwtFlag::MA;
        rbwt = rbwt::load_rowbowt<StringT>(fa, flag);
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
            auto fn = [&](typename rbwt::RowBowt<StringT>::range_t p, std::pair<size_t, size_t> q, std::vector<MarkerT> mbuf) {
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

    rbwt::RowBowt<StringT> rbwt;
    std::vector<std::vector<MarkerT>> rets;
};

*/


using Implementations =  testing::Types<ri::rle_string_sd, ri::fbb_string>;
// using Implementations =  testing::Types<ri::fbb_string>;

TYPED_TEST_SUITE(RowBowtTester, Implementations);
TYPED_TEST(RowBowtTester, Count) {
    this->CountTester();
}
TYPED_TEST(RowBowtTester, Locate) {
    this->LocateTester();
}
TYPED_TEST(RowBowtTester, Marker) {
    this->MarkerTester();
}

/* These tests take too long
TYPED_TEST(RowBowtTester, FTab) {
    this->FTabTest10mer();
}

TYPED_TEST(RowBowtTester, FTabExtend) {
    this->FTabTest11merfrom10merFtab();
}
*/