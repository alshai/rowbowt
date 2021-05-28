#ifndef ROWBOWT_IO_HPP
#define ROWBOWT_IO_HPP

#include <fstream>
#include <iostream>
#include <string>
#include <optional>
#include "rowbowt.hpp"
#include "rle_string.hpp"
#include "fbb_string.hpp"
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
    int ft_only = 0;
    size_t k = 10;
    int fbb = 0;
};

template<typename StringT=rle_string_t>
typename std::enable_if<std::is_same<StringT, ri::rle_string_sd>::value, void>::type
construct_and_serialize_rowbowt(RowBowtConstructArgs args) {
    // std::ifstream bwt_ifs(args.bwt_fname);
    StringT bwt(args.bwt_fname);
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
        rbwt::RowBowt<StringT> rb(bwt, {}, {}, {}, {});
        FTab ftab = rb.build_ftab(args.k);
        std::ofstream ft_ofs(args.prefix + ft_suffix);
        ftab.serialize(ft_ofs);
        ft_ofs.close();
    }
}

template<typename StringT=rle_string_t>
typename std::enable_if<std::is_same<StringT, ri::fbb_string>::value, void>::type
construct_and_serialize_rowbowt(RowBowtConstructArgs args) {
    StringT bwt(args.bwt_fname);
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
    if (args.tsa) {
        std::cerr << "Warning: fbb_string does not support loading toehold suffix array\n" << std::endl;
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
        rbwt::RowBowt<StringT> rb(bwt, {}, {}, {}, {});
        FTab ftab = rb.build_ftab(args.k);
        std::ofstream ft_ofs(args.prefix + ft_suffix);
        ftab.serialize(ft_ofs);
        ft_ofs.close();
    }
}

template<typename StringT=rle_string_t>
void construct_and_serialize_ftab(RowBowtConstructArgs args) {
    std::ifstream rbwt_ifs(args.prefix + rbwt_suffix);
    rbwt::RowBowt<StringT> rb;
    StringT bwt;
    if (!rbwt_ifs.good())  {
        // std::ifstream bwt_ifs(args.bwt_fname);
        bwt = StringT(args.bwt_fname);
    } else {
        std::cerr << "loading rbwt file" << std::endl;
        bwt.load(rbwt_ifs);
    }
    rb = rbwt::RowBowt<StringT>(bwt, {}, {}, {}, {});
    FTab ftab = rb.build_ftab(args.k);
    std::ofstream ft_ofs(args.prefix + ft_suffix);
    ftab.serialize(ft_ofs);
    ft_ofs.close();
}

enum class LoadRbwtFlag {
      NONE=0 // no extra aux structures
    , SA=1 // load (toehold) suffix array
    , MA=2 // load marker array
    , DL=4 // load document list
    , FT=8 // load ftab
};
inline constexpr LoadRbwtFlag operator|(LoadRbwtFlag a, LoadRbwtFlag b) {
    return static_cast<LoadRbwtFlag>(static_cast<int>(a) | static_cast<int>(b));
}
inline constexpr LoadRbwtFlag operator&(LoadRbwtFlag a, LoadRbwtFlag b) {
    return static_cast<LoadRbwtFlag>(static_cast<int>(a) & static_cast<int>(b));
}

template<typename T>
T load_obj(std::string fname) {
    T t;
    std::cerr << "loading: " << fname << std::endl;
    uint64_t x;
    std::ifstream ifs(fname, std::ios_base::in | std::ios_base::binary);
    if (!ifs.good()) {
        std::cerr << "bad file" << std::endl;
        exit(1);
    }
    t.load(ifs);
    ifs.close();
    return t;

}

template<typename StringT=rle_string_t>
RowBowt<StringT> load_rowbowt(std::string prefix, LoadRbwtFlag flag) {
    StringT bwt;
    std::ifstream bwt_ifs(prefix + rbwt_suffix);
    bwt.load(bwt_ifs);
    bwt_ifs.close();
    // this would be so much cleaner in Rust... or am I just doing it wrong?
    std::optional<ToeholdSA> tsa    = static_cast<bool>(flag & LoadRbwtFlag::SA) ? std::make_optional(load_obj<ToeholdSA>(prefix+tsa_suffix))    : std::nullopt;
    std::optional<MarkerArray<>> ma = static_cast<bool>(flag & LoadRbwtFlag::MA) ? std::make_optional(load_obj<MarkerArray<>>(prefix+ma_suffix)) : std::nullopt;
    std::optional<DocList> dl       = static_cast<bool>(flag & LoadRbwtFlag::DL) ? std::make_optional(load_obj<DocList>(prefix+dl_suffix))       : std::nullopt;
    std::optional<FTab> ft          = static_cast<bool>(flag & LoadRbwtFlag::FT) ? std::make_optional(load_obj<FTab>(prefix+ft_suffix))          : std::nullopt;
    return RowBowt<StringT>(std::move(bwt), std::move(ma), std::move(tsa), std::move(dl), std::move(ft));
}
}

#endif

