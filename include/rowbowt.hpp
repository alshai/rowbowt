#ifndef ROWBOWT_HPP
#define ROWBOWT_HPP

#include <cstdlib>
#include <ctime>
#include <cinttypes>
#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <optional>
#include <random>
#include "marker_array.hpp"
#include "toehold_sa.hpp"
#include "rle_string.hpp"
#include "doclist.hpp"
#include "ftab.hpp"

namespace rbwt {

using rle_string_t = ri::rle_string_sd;

class RowBowt {

    public:

    using range_t = std::pair<uint64_t, uint64_t>;

    RowBowt() {
    }

    RowBowt(rle_string_t& bwt
            ,const std::optional<MarkerArray<>>& ma
            ,const std::optional<ToeholdSA>& tsa
            ,const std::optional<DocList>& dl
            ,const std::optional<FTab>& ft
            )
        : bwt_(bwt)
        , f_(build_f(bwt))
        , ma_(ma)
        , tsa_(tsa)
        , dl_(dl)
        , ft_(ft)
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
        size_t i = 0;
        if (ft_) {
            range = ft_->at(query.substr(query.size()-ft_->get_k(), ft_->get_k()));
            if (range.first <= range.second) {
                i = ft_->get_k();
            } else {
                range = full_range();
            }
        }
        size_t m = query.size();
        for(i; i < m && range.second>=range.first; ++i) {
            range = LF(range,query[m-i-1]);
        }
        return range;
    }

    /*
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
    */

    struct LFData {
        LFData() {}
        LFData(range_t r, uint64_t s, uint64_t e, uint64_t ss, std::vector<MarkerT>& m)
            : rn(r)
            , qstart(s)
            , qend(e)
            , ssamp(ss)
            , markers(m)
        {}
        LFData(range_t r, uint64_t s, uint64_t e, uint64_t ss)
            : rn(r)
            , qstart(s)
            , qend(e)
            , ssamp(ss)
        {}
        LFData(range_t r, uint64_t s, uint64_t e)
            : rn(r)
            , qstart(s)
            , qend(e)
        {}
        void clear() {
            rn = {1,0};
            qstart = 0;
            qend = 0;
            ssamp = 0;
            markers.clear();
        }
        range_t rn = {1,0};
        uint64_t qstart;
        uint64_t qend;
        uint64_t ssamp;
        std::vector<MarkerT> markers;
    };

    LFData find_range_w_toehold(const std::string& query) const {
        LFData lf;
        if (!tsa_) return lf;
        uint64_t m = query.size();
        lf.rn = full_range();
        lf.ssamp = tsa_->get_last_run_sample();
        for (uint32_t i = 0; i < query.size(); ++i) {
            std::tie(lf.rn, lf.ssamp) = LF_w_loc(lf.rn, query[m-i-1], lf.ssamp);
            if (lf.rn.second < lf.rn.first) {
                lf.clear();
                return lf;
            }
        }
        // lf.rn and lf.ssamp can be used to find locations
        return lf;
    }

    // use a greedy scheme to get seeds from a query.
    // start from rightmost position of seed, then LF-map until failure. If the resulting seed length >= min_length, add
    // its range, query positions, & SA sample to results.
    // Then, skip the mismatched base start procedure again from the next base over.
    // TODO: do we want to start from the mismatched base?
    std::vector<LFData>& get_seeds_greedy(const std::string& query, uint64_t min_length, std::vector<LFData>& lfdata) const {
        lfdata.clear();
        uint64_t m = query.size();
        range_t range = full_range();
        range_t prev_range = full_range();
        uint64_t ei = m;
        for (uint64_t i = 0; i < query.size(); ++i) {
            range = LF(range, query[m-i-1]);
            if (range.second < range.first) {
                if (ei - (m-i) >= min_length) {
                    // m - i is start because we shouldn't count the current position in query
                    lfdata.push_back(LFData(prev_range, m-i, ei));
                }
                // reset everything, skip to next i
                range = full_range();
                prev_range = full_range();
                ei = m-i-1; // this makes sure we skip current position in query for the next range
            } else {
                prev_range = range;
            }
        }
        lfdata.push_back(LFData(prev_range, 0, ei));
        // range and k can be used to find locations
        return lfdata;
    }

    // use a greedy scheme to get seeds from a query.
    // start from rightmost position of seed, then LF-map until failure. If the resulting seed length >= min_length, add
    // its range, query positions, & SA sample to results.
    // Then, skip the mismatched base start procedure again from the next base over.
    // TODO: do we want to start from // mismatched base?
    std::vector<LFData>& get_seeds_greedy_w_sample(const std::string& query, uint64_t min_length, std::vector<LFData>& lfdata) const {
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
                    // std::cerr << "seed size: " << ei - (m-i) << " " << m-i << "-" << ei << std::endl;
                    lfdata.push_back(LFData(prev_range, m-i, ei, pk));
                } else {
                    // std::cerr << "(discarded) seed size: " << ei - (m-i) << " " << m-i << "-" << ei << std::endl;
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
        if (ei >= min_length) {
            // std::cerr << "seed size: " << ei << " " << 0 << "-" << ei << std::endl;
            lfdata.push_back(LFData(prev_range, 0, ei, pk));
        }
        // else {
        //     std::cerr << "(discarded) seed size: " << ei << " " << 0 << "-" << ei << std::endl;
        // }
        // range and k can be used to find locations
        // result can be empty
        return lfdata;
    }

    std::vector<LFData> get_seeds_greedy_w_sample(const std::string& query, uint64_t min_length) const {
        std::vector<LFData> lfs;
        return get_seeds_greedy_w_sample(query, min_length, lfs);
    }

    /*
     * in progress.
     * find seeds by partitioning read into min_length sized regions and picking random seeds from each region
    std::vector<LFData> get_seeds_partitioned(const std::string& query, uint64_t min_length, std::vector<LFData>& lfs) const {
        // parition the read
        uint64_t partition_length = query.size() / min_length;
        for (int i = 0; i < query.size(); i += partition_length) {
            // pick a random position between i and i + parition_length
            // apparently you're not supposed to use mod with rand()
            // https://en.cppreference.com/w/cpp/numeric/random/rand
            uint64_t start_pos = i + (std::rand()/((RAND_MAX + 1u)/ parition_length));
            // do LF starting from start_pos
        }
    }
    */

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

    LFData find_range_w_markers(const std::string query, uint64_t wsize, uint64_t max_range) const {
        LFData lf;
        if (!ma_) {
            std::cerr << "warning: no marker array found!\n";
            return lf;
        }
        lf.markers.clear();
        if (query.size() < wsize) {
            std::cerr << "warning: query (size=" << query.size() << ") is less than wsize (" << wsize << ")\n";
            return lf;
        }
        lf.rn = full_range();
        uint64_t m = query.size();
        uint64_t window_ei = m;
        std::vector<MarkerT> mbuf;
        uint64_t i = 0;
        size_t nqueries = 0;
        for (i = 0; i < m; ++i) {
            lf.rn = LF(lf.rn, query[m-i-1]);
            if (lf.rn.second < lf.rn.first) {
                lf.clear();
                // std::cerr << "bad range at " << i << ": (" << lf.rn.first << ", " << lf.rn.second << ")\n";
                return lf;
            } else { // deal with windows here
                if (window_ei-(m-i) >= wsize) {
                    nqueries += 1;
                    mbuf.clear();
                    if (lf.rn.second - lf.rn.first + 1 <= max_range) {
                        mbuf = markers_at(lf.rn, mbuf);
                        lf.markers.insert(lf.markers.begin(), mbuf.begin(), mbuf.end());
                    }
                    window_ei = m-i; // current position is now window end
                }
            }
        }
        // deal with last bit (regardless of window size)
        // TODO: avoid redundancy here is query length % wsize == 0
        if (lf.rn.second >= lf.rn.first && (m - 1) % wsize != 0) {
            nqueries += 1;
            mbuf.clear();
            if (lf.rn.second - lf.rn.first + 1 <= max_range) {
                mbuf = markers_at(lf.rn, mbuf);
                lf.markers.insert(lf.markers.begin(), mbuf.begin(), mbuf.end());
            }
        }
        lf.qstart = 0;
        lf.qend = m;
        return lf;
    }

    // TODO: include counts/ranges in return value as well
    template<typename F>
    void get_markers_lmems(const std::string query, uint64_t wsize, uint64_t max_range, F fn) const {
        uint64_t m = query.size();
        for (int j = 0; j < query.size(); ++j) {
            range_t prev_range = full_range(), range = full_range();
            uint64_t window_ei = m-j, seed_ei = m-j;
            std::vector<MarkerT> mbuf;
            uint64_t i = 0;
            auto update_mbuf = [&](range_t r) {
                mbuf.clear();
                if (r.second - r.first + 1 <= max_range) {
                    mbuf = markers_at(r, mbuf);
                }
                fn(r, std::make_pair(m-i,seed_ei-1), mbuf);
            };
            for (i = j; i < query.size(); ++i) {
                range = LF(range, query[m-i-1]);
                if (range.second < range.first) { // this is when the seed fails
                    if (seed_ei-(m-i) >= wsize) { // check markers here if seed is large enough, regardless of window length
                        update_mbuf(prev_range);
                    }
                    break;
                } else { // this is for each window
                    if (window_ei-(m-i) >= wsize) {
                        update_mbuf(range);
                        window_ei = m-i; // current position is now window end
                    }
                    prev_range = range;
                }
            }
            std::cerr << m-j << ": " << prev_range.second-prev_range.first+1 << std::endl;
            // this is when the whole read is done and a seed hasn't finished yet
            if (seed_ei-(m-i) >= wsize) {
                update_mbuf(range);
            }
        }
    }

    // TODO: include counts/ranges in return value as well
    template<typename F>
    void get_markers_greedy_seeding(const std::string query, uint64_t wsize, uint64_t max_range, F fn) const {
        uint64_t m = query.size();
        range_t prev_range = full_range(), range = full_range();
        uint64_t window_ei = m, seed_ei = m;
        std::vector<MarkerT> mbuf;
        uint64_t i = 0;
        auto update_mbuf = [&](range_t r) {
            mbuf.clear();
            if (r.second - r.first + 1 <= max_range) {
                mbuf = markers_at(r, mbuf);
            }
            fn(r, std::make_pair(m-i,seed_ei-1), mbuf);
        };
        for (i = 0; i < query.size(); ++i) {
            range = LF(range, query[m-i-1]);
            if (range.second < range.first) { // this is when the seed fails
                if (seed_ei-(m-i) >= wsize) { // check markers here if seed is large enough, regardless of window length
                    update_mbuf(prev_range);
                } // then reset the seed
                range = full_range();
                prev_range = full_range();
                seed_ei = m-i-1; // this makes sure we skip current position in query for the next range
                window_ei = m-i-1;
            } else { // this is for each window
                if (window_ei-(m-i) >= wsize) {
                    update_mbuf(range);
                    window_ei = m-i; // current position is now window end
                }
                prev_range = range;
            }
        }

        // this is when the whole read is done and a seed hasn't finished yet
        if (seed_ei-(m-i) >= wsize) {
            update_mbuf(range);
        }
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

    std::vector<LFData>& find_range_w_toehold_chkpnts(const std::string& query, uint64_t wsize, std::vector<LFData>& lfs) const {
        lfs.clear();
        if (!tsa_) return lfs;
        uint64_t m = query.size();
        uint64_t window_ei = m;
        uint64_t i = 0;
        LFData lf;
        lf.rn = full_range();
        lf.ssamp = tsa_->get_last_run_sample();
        for (i = 0; i < query.size(); ++i) {
            std::tie(lf.rn, lf.ssamp) = LF_w_loc(lf.rn, query[m-i-1], lf.ssamp);
            if (lf.rn.second < lf.rn.first) {
                lfs.clear();
                return lfs;
            } else {
                if (window_ei-(m-i) >= wsize) {
                    lf.qstart = m-i;
                    lf.qend = window_ei;
                    lfs.push_back(LFData(lf));
                    window_ei = m-i; // current position is now window end
                }
            }
        }
        if (lf.rn.second >= lf.rn.first && (m - 1) % wsize != 0) {
            lf.qstart = 0;
            lf.qend = m;
            lfs.push_back(LFData(lf));
        }
        return lfs;
    }

    std::vector<LFData> find_range_w_toehold_chkpnts(const std::string& query, uint64_t wsize) const {
        std::vector<LFData> lfs;
        return find_range_w_toehold_chkpnts(query, wsize, lfs) ;
    }

    std::vector<uint64_t>& locs_at(range_t range, uint64_t k, uint64_t max_hits, std::vector<uint64_t>& locs) const {
        return tsa_->locate_range(range.first, range.second, k, max_hits, locs);
    }

    std::vector<uint64_t> locs_at(range_t range, uint64_t k, uint64_t max_hits) const {
        std::vector<uint64_t> locs;
        tsa_->locate_range(range.first, range.second, k, max_hits, locs);
        return locs;
    }

    std::pair<std::string, uint64_t> resolve_offset(uint64_t i) const {
        return dl_->doc_and_offset_at(i);
    }

    std::vector<uint64_t>& find_locs(std::string s, uint64_t max_hits, std::vector<uint64_t>& locs) const {
        auto lf = find_range_w_toehold(s);
        locs_at(lf.rn, lf.ssamp, max_hits, locs);
        return locs;
    }

    std::vector<uint64_t>& find_locs_greedy_seeding(std::string s, uint64_t min_length, uint64_t max_hits, std::vector<uint64_t>& locs) const {
        locs.clear();
        std::vector<LFData> lfs;
        lfs = get_seeds_greedy_w_sample(s, min_length, lfs);
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

    std::vector<uint64_t> find_locs_greedy_seeding(std::string s, uint64_t min_length, uint64_t max_hits) const {
        std::vector<uint64_t> locs;
        return find_locs_greedy_seeding(s, min_length, max_hits, locs);
    }

    std::vector<uint64_t>& locate_from_longest_seed(uint64_t max_hits, const std::vector<LFData>& lfs, std::vector<uint64_t>& locs) const {
        locs.clear();
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
            }
        }
        // std::cerr << "best seed [" << best_range.rn.first << " " << best_range.rn.second << "), query [" << best_range.qend-best_range.qstart <<  "] \n";
        // locate for best lf
        locs = locs_at(best_range.rn, best_range.ssamp, max_hits, locs);
        // correction based on where the seed is in the read
        for (uint64_t i = 0; i < locs.size(); ++i) {
            locs[i] = locs[i] - best_range.qstart;
        }
        return locs;
    }

    std::vector<uint64_t> locate_from_longest_seed(uint64_t max_hits, const std::vector<LFData>& lfs) {
        std::vector<uint64_t> locs;
        return locate_from_longest_seed(max_hits, lfs, locs);
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

    // kmer generation code from from heng li:
    // https://www.biostars.org/p/18096/#18107
    FTab build_ftab(size_t k = 10) {
        int i;
        unsigned long long x, y;
        std::string kmer(k, ' ');
        kmer.clear();
        FTab ftab;
        for (x = 0; x < 1ULL<<(2*k); ++x) {
            for (i = 0, y = x; i < k; ++i, y >>= 2) {
                kmer += "ACGT"[y&3];
            }
            auto range = find_range(kmer);
            if (range.first <= range.second) {
                ftab[kmer] = range;
            }
            kmer.clear();
        }
        return ftab;
    }

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
    std::optional<FTab> ft_;
    uint64_t r_;
};
}

#endif
