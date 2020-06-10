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

namespace rbwt {

std::string rbwt_suffix = ".rbwt";
std::string tsa_suffix = ".tsa";
std::string ma_suffix = ".mab";

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
    std::string prefix;
    int ma = false;
    int tsa = false;
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
        std::cout << "constructing tsa with args: " << bwt.size() << " " << bwt.number_of_runs()  << " " << args.ssa_fname << " " << args.esa_fname << std::endl;
        ToeholdSA tsa(bwt.size(), bwt.number_of_runs(), args.ssa_fname, args.esa_fname);
        std::cout << "saving tsa to " << args.prefix + tsa_suffix << std::endl;
        std::ofstream tsa_ofs(args.prefix + tsa_suffix);
        tsa.serialize(tsa_ofs);
        tsa_ofs.close();
    }
}

enum class LoadRbwtFlag {NONE=0, SA=1, MA=2};
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
    std::ifstream bwt_ifs(prefix + rbwt_suffix);
    bwt.load(bwt_ifs);
    bwt_ifs.close();
    if (static_cast<bool>(flag & LoadRbwtFlag::SA)) {
        std::ifstream tsa_ifs(prefix + tsa_suffix);
        tsa.load(tsa_ifs);
        tsa_ifs.close();
    } if (static_cast<bool>(flag & LoadRbwtFlag::MA)) {
        std::ifstream ma_ifs(prefix + ma_suffix);
        ma.load(ma_ifs);
        ma_ifs.close();
    }
    return RowBowt(bwt, std::make_optional(ma), std::make_optional(tsa));
}

}
#endif
