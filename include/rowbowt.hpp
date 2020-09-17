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
#include "doclist.hpp"

namespace rbwt {

using rle_string_t = ri::rle_string_sd;

class RowBowt {

    public:

    using range_t = std::pair<uint64_t, uint64_t>;

    RowBowt() { }

    RowBowt(rle_string_t& bwt
            ,const std::optional<MarkerArray<>>& ma
            ,const std::optional<ToeholdSA>& tsa
            ,const std::optional<DocList>& dl
            )
        : bwt_(bwt)
        , f_(build_f(bwt))
        , ma_(ma)
        , tsa_(tsa)
        , dl_(dl)
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
                return std::make_pair(range, 0);
            }
        }
        // range and k can be used to find locations
        return std::make_pair(range, k);
    }


    struct LFData {
        LFData() {};
        LFData(range_t r, uint64_t s, uint64_t e, uint64_t ss)
            : rn(r)
            , qstart(s)
            , qend(e)
            , ssamp(ss)
        {}
        range_t rn;
        uint64_t qstart;
        uint64_t qend;
        uint64_t ssamp;
    };

    // Return BWT range along with a 'toehold' sample of the suffix array. This toehold can be used find neighboring
    // values in the suffix array
    std::vector<LFData>& find_range_w_toehold_greedy(const std::string& query, uint64_t min_length, std::vector<LFData>& lfdata) const {
        lfdata.clear();
        if (!tsa_) return lfdata;
        uint64_t m = query.size();
        range_t range = full_range();
        range_t prev_range = full_range();
        const uint64_t first_k = tsa_->get_last_run_sample();
        uint64_t k = first_k;
        uint64_t pk = -1;
        uint64_t ei = m;
        for (uint64_t i = 0; i < query.size(); ++i) {
            std::tie(range, k) = LF_w_loc(range, query[m-i-1], k);
            if (range.second < range.first) {
                if (ei - (m-i) >= min_length) {
                    // m - i is start because we shouldn't count the current position in query
                    lfdata.push_back(LFData(prev_range, m-i, ei, pk));
                }
                // reset everything, skip to next i
                k = first_k;
                range = full_range();
                prev_range = full_range();
                ei = m-i-1; // this makes sure we skip current position in query for the next range
            } else {
                prev_range = range;
                pk = k;
            }
        }
        lfdata.push_back(LFData(prev_range, 0, ei, pk));
        // range and k can be used to find locations
        return lfdata;
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

    std::vector<uint64_t>& locs_at(range_t range, uint64_t k, uint64_t max_hits, std::vector<uint64_t>& locs) const {
        return tsa_->locate_range(range.first, range.second, k, max_hits, locs);
    }

    std::vector<uint64_t>& locs_at(range_t range, uint64_t k, uint64_t max_hits) const {
        std::vector<uint64_t> locs;
        return tsa_->locate_range(range.first, range.second, k, max_hits, locs);
    }

    std::pair<std::string, uint64_t> resolve_offset(uint64_t i) const {
        return dl_->doc_and_offset_at(i);
    }

    std::vector<uint64_t> locate_pattern(std::string s, uint64_t max_hits, std::vector<uint64_t>& locs) {
        auto range_k = find_range_w_toehold(s);
        return locs_at(range_k.first, range_k.second, max_hits, locs);
    }

    std::vector<uint64_t>& locate_pattern_greedy_seeding(std::string s, uint64_t min_length, uint64_t max_hits, std::vector<uint64_t>& locs) {
        locs.clear();
        std::vector<LFData> lfs;
        lfs = find_range_w_toehold_greedy(s, min_length, lfs);
        if (!lfs.size()) {
            return locs;
        }
        LFData best_range;
        uint64_t max_length = 0;
        for (const auto& lfd: lfs) {
            auto length = lfd.qend - lfd.qstart;
            if (length > max_length) {
                max_length = length;
                best_range = lfd;
                // std::cout << "best seed [" << lfd.rn.first << " " << lfd.rn.second << "), query [" << lfd.qstart << " " << lfd.qend << ") \n";
            }
        }
        // locate for best lf
        locs = locs_at(best_range.rn, best_range.ssamp, max_hits, locs);
        // correction based on where the seed is in the read
        for (uint64_t i = 0; i < locs.size(); ++i) {
            locs[i] = locs[i] - best_range.qstart;
        }
        return locs;
    }

    std::vector<uint64_t> locate_pattern_greedy_seeding(std::string s, uint64_t min_length, uint64_t max_hits) {
        std::vector<uint64_t> locs;
        return locate_pattern_greedy_seeding(s, min_length, max_hits, locs);
    }


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
    std::optional<DocList> dl_;
    uint64_t r_;
};
}

#endif
