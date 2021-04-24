#ifndef ROWBOWT_IO_HPP
#define ROWBOWT_IO_HPP

#include <fstream>
#include <iostream>
#include <string>
#include <optional>
#include "rowbowt.hpp"
#include "rle_string.hpp"
#include "toehold_sa.hpp"
#include "marker_array.hpp"
#include "ftab.hpp"

namespace rbwt {

std::string rbwt_suffix = ".rbwt";
std::string tsa_suffix = ".tsa";
std::string ma_suffix = ".mab";
std::string dl_suffix = ".docs";
std::string ft_suffix = ".ftab";

bool file_exists(std::string fname) {
    std::ifstream ifs(fname);
    return ifs.good();
}

void file_ne_error(std::string fname) {
    std::cerr << fname << " does not exist!\n";
    exit(1);
}

struct RowBowtConstructArgs {
    std::string bwt_fname;
    std::string ssa_fname;
    std::string esa_fname;
    std::string ma_fname;
    std::string dl_fname;
    std::string prefix;
    int ma = 0;
    int tsa = 0;
    int dl = 0;
    int ft = 0;
    size_t k = 10;
};

void construct_and_serialize_rowbowt(RowBowtConstructArgs args) {
    std::ifstream bwt_ifs(args.bwt_fname);
    rle_string_t bwt(bwt_ifs);
    std::ofstream bwt_ofs(args.prefix + rbwt_suffix);
    bwt.serialize(bwt_ofs);
    bwt_ofs.close();
    if (args.ma && args.ma_fname != "") {
        if (!file_exists(args.ma_fname)) file_ne_error(args.ma_fname);
        MarkerArray<> ma(args.ma_fname);
        std::ofstream ma_ofs(args.prefix + ma_suffix);
        ma.serialize(ma_ofs);
        ma_ofs.close();
    }
    if (args.tsa && args.ssa_fname != "" && args.esa_fname != "") {
        if (!file_exists(args.ssa_fname)) file_ne_error(args.ssa_fname);
        if (!file_exists(args.esa_fname)) file_ne_error(args.esa_fname);
        ToeholdSA tsa(bwt.size(), bwt.number_of_runs(), args.ssa_fname, args.esa_fname);
        std::ofstream tsa_ofs(args.prefix + tsa_suffix);
        tsa.serialize(tsa_ofs);
        tsa_ofs.close();
    }
    if (args.dl && args.dl_fname != "") {
        std::ifstream ifs(args.dl_fname);
        if (args.dl_fname != args.prefix + dl_suffix) {
            std::ofstream ofs(args.prefix + dl_suffix);
            ofs << ifs.rdbuf();
            ofs.close();
        }
        ifs.close();
    }
    if (args.ft) {
        rbwt::RowBowt rb(bwt, {}, {}, {}, {});
        FTab ftab = rb.build_ftab(args.k);
        std::ofstream ft_ofs(args.prefix + ft_suffix);
        ftab.serialize(ft_ofs);
        ft_ofs.close();
    }
}

enum class LoadRbwtFlag {
      NONE // no extra aux structures
    , SA // load (toehold) suffix array
    , MA // load marker array
    , DL // load document list
    , FT // load ftab
};
inline constexpr LoadRbwtFlag operator|(LoadRbwtFlag a, LoadRbwtFlag b) {
    return static_cast<LoadRbwtFlag>(static_cast<int>(a) | static_cast<int>(b));
}
inline constexpr LoadRbwtFlag operator&(LoadRbwtFlag a, LoadRbwtFlag b) {
    return static_cast<LoadRbwtFlag>(static_cast<int>(a) & static_cast<int>(b));
}

RowBowt load_rowbowt(std::string prefix, LoadRbwtFlag flag) {
    rle_string_t bwt;
    MarkerArray<> ma;
    ToeholdSA tsa;
    DocList dl;
    FTab ft;
    std::ifstream bwt_ifs(prefix + rbwt_suffix);
    bwt.load(bwt_ifs);
    bwt_ifs.close();
    if (static_cast<bool>(flag & LoadRbwtFlag::SA)) {
        std::cerr << "loading: " << prefix + tsa_suffix << std::endl;
        uint64_t x;
        std::ifstream tsa_ifs(prefix + tsa_suffix, std::ios_base::in | std::ios_base::binary);
        if (!tsa_ifs.good()) {
            std::cerr << "bad tsa file" << std::endl;
            exit(1);
        }
        tsa.load(tsa_ifs);
        tsa_ifs.close();
    }
    if (static_cast<bool>(flag & LoadRbwtFlag::MA)) {
        std::cerr << "loading: " << prefix + ma_suffix << std::endl;
        std::ifstream ma_ifs(prefix + ma_suffix);
        if (!ma_ifs.good()) {
            std::cerr << "bad mab file" << std::endl;
            exit(1);
        }
        ma.load(ma_ifs);
        ma_ifs.close();
    }
    if (static_cast<bool>(flag & LoadRbwtFlag::DL)) {
        std::cerr << "loading: " << prefix + dl_suffix << std::endl;
        std::ifstream dl_ifs(prefix + dl_suffix);
        if (!dl_ifs.good()) {
            std::cerr << "bad docs file" << std::endl;
            exit(1);
        }
        dl.load(dl_ifs);
        dl_ifs.close();
    }
    if (static_cast<bool>(flag & LoadRbwtFlag::FT)) {
        std::cerr << "loading: " << prefix + ft_suffix << std::endl;
        std::ifstream ft_ifs(prefix + ft_suffix);
        if (!ft_ifs.good()) {
            std::cerr << "bad ftab file" << std::endl;
            exit(1);
        }
        ft.load(ft_ifs);
        ft_ifs.close();
    }
    return RowBowt(bwt, std::make_optional(ma), std::make_optional(tsa), std::make_optional(dl), std::make_optional(ft));
}
}

#endif

