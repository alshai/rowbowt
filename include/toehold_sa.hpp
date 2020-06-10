#ifndef TOEHOLD_HPP
#define TOEHOLD_HPP

/* Header file describing a run-length compressed SA that requires a "toehold" for access
 * Prezza et. al Optimal-Time Text Indexing in BWT-runs Bounded Space. 2017. https://arxiv.org/abs/1705.10382
 * Code heavily derived from Nicola Prezza's r-index repo: https://github.com/nicolaprezza/r-index
 * NOTE: currently only supports the Phi function (citation here)
 *
 * TODO: use memory mapping when doing construction
 */

#include <string>
#include <cinttypes>
#include <fstream>
#include <vector>
#include <algorithm>
#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>
#include "sparse_sd_vector.hpp"
#include "utils.hpp"

class ToeholdSA {

    public:

    ToeholdSA() {}
    ToeholdSA(uint64_t n, uint64_t r, std::string ssa_fname, std::string esa_fname)
        : n_(n)
        , r_(r) {
        std::vector<std::pair<uint64_t, uint64_t>> samples_first_vec;
        std::vector<uint64_t> samples_last_vec;
        read_run_starts(ssa_fname, n, samples_first_vec);
        read_run_ends(esa_fname, n, samples_last_vec);
        build_phi(samples_first_vec, samples_last_vec);
    }

    std::vector<uint64_t>& locate_range(uint64_t l, uint64_t r, uint64_t k, uint64_t max_hits, vector<uint64_t>& locs) const {
        uint64_t n_occ = r>=l ? (r-l)+1 : 0;
        if (n_occ > max_hits) n_occ = max_hits;
        uint64_t k1 = k;
        if (n_occ > 0) {
            locs.push_back(k1);
            for(uint64_t i = 1; i < n_occ; ++i) {
                k1 = (this->phi)(k1);
                locs.push_back(k1);
            }
        }
        return locs;
    }

    std::vector<uint64_t>& locate_range(uint64_t l, uint64_t r, uint64_t k, uint64_t max_hits) const {
        std::vector<uint64_t> locs;
        return locate_range(l, r, k, max_hits, locs);
    }

    uint64_t phi(uint64_t i) const {
        assert(i != n_-1);
        //jr is the rank of the predecessor of i (circular)
        uint64_t jr = pred_.predecessor_rank_circular(i);
        assert(jr<=r_-1);
        //the actual predecessor
        uint64_t j = pred_.select(jr);
        assert(jr<r_-1 or j == n_-1);
        //distance from predecessor
        uint64_t delta = j<i ? i-j : i+1;
        //cannot fall on first run: this can happen only if I call Phi(SA[0])
        assert(pred_to_run_[jr]>0);
        //sample at the end of previous run
        assert(pred_to_run_[jr]-1 < samples_last_.size());
        uint64_t prev_sample = samples_last_[ pred_to_run_[jr]-1 ];
        return (prev_sample + delta) % n_;
    }

    size_t serialize(std::ostream& out) {
        size_t nbytes  = 0;
        nbytes = 2 * sizeof(uint64_t);
        out.write(reinterpret_cast<char*>(&r_), sizeof(uint64_t));
        out.write(reinterpret_cast<char*>(&n_), sizeof(uint64_t));
        nbytes += pred_.serialize(out);
        nbytes += samples_last_.serialize(out);
        nbytes += pred_to_run_.serialize(out);
        return nbytes;
    }

    void load(std::istream& in) {
        in.read(reinterpret_cast<char*>(&r_), sizeof(uint64_t));
        in.read(reinterpret_cast<char*>(&n_), sizeof(uint64_t));
        pred_.load(in);
        samples_last_.load(in);
        pred_to_run_.load(in);
    }

    uint64_t samples_last_at(uint64_t i) const {
        return samples_last_[i];
    }

    uint64_t get_last_run_sample() const {
        return (samples_last_[r_-1]+1) % n_;
    }

    private:

    using sparse_bv_type = ri::sparse_sd_vector;

    void build_phi(std::vector<std::pair<uint64_t,uint64_t>>& samples_first_vec, std::vector<uint64_t>& samples_last_vec) {
        int log_r = ri::utils::bitsize(uint64_t(r_));
        int log_n = ri::utils::bitsize(uint64_t(n_));
        std::sort(samples_first_vec.begin(), samples_first_vec.end());
        //build Elias-Fano predecessor
        {
            auto pred_bv = std::vector<bool>(n_);
            for(auto p : samples_first_vec) {
                assert(p.first < pred_bv.size());
                pred_bv[p.first] = true;
            }
            pred_ = sparse_bv_type(pred_bv);
        }
        assert(pred_.rank(pred_.size()) == r_);
        //last text position must be sampled
        assert(pred_[pred_.size()-1]);
        samples_last_ = sdsl::int_vector<>(r_,0,log_n); //text positions corresponding to last characters in BWT runs, in BWT order
        pred_to_run_ = sdsl::int_vector<>(r_,0,log_r); //stores the BWT run (0...R-1) corresponding to each position in pred, in text order
        for(uint64_t i=0;i<samples_last_vec.size();++i) {
            assert(ri::utils::bitsize(uint64_t(samples_last_vec[i])) <= log_n);
            samples_last_[i] = samples_last_vec[i];
        }
        for(uint64_t i=0;i<samples_first_vec.size();++i) {
            assert(ri::utils::bitsize(uint64_t(samples_first_vec[i].second)) <= log_r);
            pred_to_run_[i] = samples_first_vec[i].second;
        }
    }

    std::vector<std::pair<uint64_t,uint64_t>>& read_run_starts(std::string fname, uint64_t n, std::vector<std::pair<uint64_t,uint64_t>>& ssa) {
        ssa.clear();
        std::ifstream ifs(fname);
        uint64_t x = 0;
        uint64_t y = 0;
        uint64_t i = 0;
        while (ifs.read(reinterpret_cast<char*>(&x), 8) && ifs.read(reinterpret_cast<char*>(&y), 8)) {
            ssa.push_back(std::pair<uint64_t,uint64_t>(y ? y-1 : n-1, i));
            i++;
        }
        return ssa;
    }

    std::vector<uint64_t>& read_run_ends(std::string fname, uint64_t n, std::vector<uint64_t>& esa) {
        esa.clear();
        std::ifstream ifs(fname);
        uint64_t x = 0;
        uint64_t y = 0;
        while (ifs.read(reinterpret_cast<char*>(&x), 8) && ifs.read(reinterpret_cast<char*>(&y), 8)) {
            esa.push_back(y ? y-1 : n-1);
        }
        return esa;
    }

    sparse_bv_type pred_;
    sdsl::int_vector<> samples_last_; //text positions corresponding to last characters in BWT runs, in BWT order
    sdsl::int_vector<> pred_to_run_; //stores the BWT run (0...R-1) corresponding to each position in pred, in text order
    uint64_t r_;
    uint64_t n_;

};

#endif
