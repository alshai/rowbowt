#ifndef ROWBOWT_HPP
#define ROWBOWT_HPP

#include <cstdlib>
#include <cinttypes>
#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <optional>
#include "marker_array.hpp"
#include "toehold_sa.hpp"
#include "rle_string.hpp"
#include "sdsl_bv_wrappers.hpp"
#include "sdsl/bit_vectors.hpp"

namespace rbwt {

using rle_string_t = ri::rle_string_sd;

class RowBowt {

    public:

    using range_t = std::pair<uint64_t, uint64_t>;

    RowBowt() { }
    RowBowt(rle_string_t& bwt
            ,const std::optional<MarkerArray<>>& ma
            ,const std::optional<ToeholdSA>& tsa
            )
        : bwt_(bwt)
        , f_(build_f(bwt))
        , ma_(ma)
        , tsa_(tsa)
    {
        r_ = bwt_.number_of_runs();
    }


    //backward navigation of BWT
    uint64_t LF(uint64_t  i) const {
        auto c = bwt_[i];
        return f_[c] + bwt_.rank(i,c);
    }

    // \param r inclusive range of a string w
    // \param c character
    // \return inclusive range of cw
    //
    range_t LF(range_t rn, uint8_t c) const {
        //if character does not appear in the text, return empty pair
        if((c==255 and f_[c]==bwt_.size()) || f_[c]>=f_[c+1])
            return {1,0};
        //number of c before the interval
        uint64_t c_before = bwt_.rank(rn.first,c);
        // number of c inside the interval rn
        // if rn.1 == rn.2, then this redudant rank query _should_ be avoided,
        // but this case won't be encountered too often if n/r is high
        uint64_t c_inside = bwt_.rank(rn.second+1,c) - c_before;
        //if there are no c in the interval, return empty range
        if(c_inside==0) return {1,0};
        uint64_t l = f_[c] + c_before;
        return {l,l+c_inside-1};
    }

    range_t LF(uint64_t s, uint64_t e, uint8_t c) const {
        return LF(range_t(s, e), c);
    }

    /*

    //forward navigation of the BWT
    uint64_t FL(uint64_t  i) {
        //i-th character in first BWT column
        auto c = f_at(i);
        //this c is the j-th (counting from 0)
        uint64_t j = i - f_[c];
        return bwt_.select(j,uint8_t(c));
    }

    //forward navigation of the BWT, where for efficiency we give c=F[i] as input
    uint64_t FL(uint64_t  i, uint8_t c) {
        //i-th character in first BWT column
        assert(c == f_at(i));
        //this c is the j-th (counting from 0)
        uint64_t j = i - f_[c];
        return bwt_.select(j,uint8_t(c));
    }
    */

    range_t full_range() const {
        //inclusive range
        return {0,bwt_.size()-1};
    }

    // Return BWT range of query string
    range_t find_range(const std::string& query) const {
        range_t range = full_range();
        uint64_t m = query.size();
        for(uint64_t i=0;i<m and range.second>=range.first;++i)
            range = LF(range,query[m-i-1]);
        return range;
    }

    // Return BWT range along with a 'toehold' sample of the suffix array. This toehold can be used find neighboring
    // values in the suffix array
    std::pair<range_t, uint64_t> find_range_w_toehold(const std::string& query) const {
        if (!tsa_) return {range_t(), 0};
        uint64_t m = query.size();
        range_t range = full_range();
        uint64_t k = tsa_->get_last_run_sample();
        for (uint32_t i = 0; i < query.size(); ++i) {
            std::tie(range, k) = LF_w_loc(range, query[m-i-1], k);
            if (range.second < range.first) {
                return std::make_pair(range_t(), 0);
            }
        }
        // range and k can be used to find locations
        return std::make_pair(range, k);
    }

    //
    // Return number of occurrences of query in the text
    ///
    uint64_t count(const std::string& query) const {
        auto rn = find_range(query);
        return rn.second>=rn.first ? (rn.second-rn.first)+1 : 0;
    }

    std::vector<MarkerT>& markers_at(uint64_t i, std::vector<MarkerT>& markers) const {
        if (!ma_) return markers;
        return ma_->at(i, markers);
    }

    std::vector<MarkerT> markers_at(uint64_t i) const {
        std::vector<MarkerT> markers;
        return markers_at(i, markers);
    }

    std::vector<MarkerT>& markers_at(range_t r, std::vector<MarkerT>& markers) const {
        if (!ma_) return markers;
        return ma_->at_range(r.first, r.second, markers);
    }

    std::vector<MarkerT> markers_at(range_t r) const {
        std::vector<MarkerT> markers;
        return markers_at(r, markers);
    }

    std::pair<range_t, uint64_t> LF_w_loc(const range_t range, uint8_t c, uint64_t k) const {
        uint64_t nk;
        range_t nrange = LF(range, c);
        if (nrange.first <= nrange.second) {
            if (bwt_[range.second] == c) { // trivial case
                nk = k-1;
            } else { // need to sample again
                uint64_t rnk = bwt_.rank(range.second, c) - 1;
                uint64_t j = bwt_.select(rnk, c);
                uint64_t run_of_j = bwt_.run_of_position(j);
                nk = tsa_->samples_last_at(run_of_j); // TODO: rename this function?
            }
        } else {
            return {{1,0},0};
        }
        return {nrange, nk};
    }

    std::vector<uint64_t>& locs_at(range_t range, uint64_t k, uint64_t max_hits, vector<uint64_t>& locs) const {
        return tsa_->locate_range(range.first, range.second, k, max_hits, locs);
    }

    std::vector<uint64_t>& locs_at(range_t range, uint64_t k, uint64_t max_hits) const {
        std::vector<uint64_t> locs;
        return tsa_->locate_range(range.first, range.second, k, max_hits, locs);
    }


    /*

    std::vector<uint64_t> locate_range(range_t range, uint64_t k, uint64_t max_hits = -1) {
    }

    std::vector<uint64_t> locate_pattern(std::string s, uint64_t max_hits = -1) {
    }

    */


    /* we should expect the user to serialize during construction */
    /*
    size_t serialize_bwt(std::string out_fname) {
        return bwt_.serialize(std::ostream(out_fname));
    }

    size_t serialize_toehold_sa(std::string out_fname) {
        return tsa_.serialize(std::ostream(out_fname));
    }

    size_t serialize_marker_array(std::string out_fname) {
        return ma_.serialize(std::ostream(out_fname));
    }
    */

    void load_bwt(std::ifstream& ifs) {
        bwt_.load(ifs);
        build_f(bwt_);
    }

    void load_toehold_sa(std::ifstream& ifs) {
        tsa_->load(ifs);
    }

    void load_marker_array(std::ifstream& ifs) {
        ma_->load(ifs);
    }

    /*
    void load_doc_list(std::ifstream& ifs) {
        docs_->load(ifs);
    }
    */

    const std::vector<uint64_t>&        get_f()              const {return f_;}
    const std::optional<ToeholdSA>&     get_toehold_sa()     const {return tsa_;}
    const std::optional<MarkerArray<>>& get_marker_array()   const {return ma_;}
    // const std::optional<DocList>&       get_doc_list()       const {return docs_;}

    private:

    std::vector<uint64_t> build_f(ri::rle_string_sd& bwt) const {
        std::vector<uint64_t> f_(256,0);
        uint64_t p = 0;
        for (size_t i = 0; i < 255; ++i) {
            p += bwt.rank(bwt.size(), i);
            f_[i+1] = p;
        }
        return f_;
    }

    uint8_t f_at(uint64_t i) const {
        uint64_t c = (std::upper_bound(f_.begin(), f_.end(),i) - f_.begin()) - 1;
        assert(c<256);
        assert(i>=f_[c]);
        return uint8_t(c);
    }


    rle_string_t bwt_;
    std::vector<uint64_t> f_;
    std::optional<ToeholdSA> tsa_;
    std::optional<MarkerArray<>> ma_;
    // std::optional<DocList> docs_;
    uint64_t r_;
};

/*
constexpr uint8_t TERMINATOR = 1;
std::vector<uint64_t> build_f(std::ifstream& ifs) {
    std::vector<uint64_t> f_(256,0);
    uint8_t c;
    uint64_t i = 0;
    while (ifs >> c) {
        if (c>TERMINATOR) f_[c]++;
        else f_[TERMINATOR]++;
        i++;
    }
    for(uint64_t i=255;i>0;--i)
        f_[i] = f_[i-1];
    f_[0] = 0;
    for(uint64_t i=1;i<256;++i)
        f_[i] += f_[i-1];
    return f_;
}
*/
}
#endif
