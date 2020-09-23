#ifndef DOCLIST_HPP
#define DOCLIST_HPP

#include <string>
#include <cinttypes>
#include "sdsl/bit_vectors.hpp"
#include "sdsl_bv_wrappers.hpp"

class DocList {

    public:

    DocList() {}

    DocList(const DocList& rhs)
        : doc_names_(rhs.doc_names_)
        , doc_bounds_(rhs.doc_bounds_)
    { }

    DocList(DocList&& rhs)
        : doc_names_(std::move(rhs.doc_names_))
        , doc_bounds_(std::move(rhs.doc_bounds_))
    { }

    DocList(std::istream& in) {
        load(in);
    }

    DocList& operator=(const DocList& rhs) {
        doc_names_ = rhs.doc_names_;
        doc_bounds_ = rhs.doc_bounds_;
        return *this;
    }

    DocList& operator=(DocList&& rhs) {
        doc_names_ = std::move(rhs.doc_names_);
        doc_bounds_ = std::move(rhs.doc_bounds_);
        return *this;
    }

    std::string doc_at(uint64_t i) const {
        uint64_t rank = doc_bounds_rank(i);
        return doc_names_[rank-1];
    }

    std::pair<std::string, uint64_t> doc_and_offset_at(uint64_t i) const {
        uint64_t rank = doc_bounds_rank(i);
        size_t offset = i - doc_bounds_.select(rank);
        return std::make_pair(doc_names_[rank-1], offset);
    }

    uint64_t offset_at(uint64_t i) const {
        uint64_t rank = doc_bounds_rank(i);
        return i - doc_bounds_.select(rank);
    }

    void load(std::istream& ifs) {
        std::string name;
        uint64_t pos = 0;
        std::vector<std::string> namev;
        std::vector<uint64_t> posv;
        while (ifs >> name >> pos) { // each doc name associated w/ start pos
            namev.push_back(name);
            posv.push_back(pos);
        }
        sdsl::bit_vector bv(pos+1);
        for (auto p: posv) {
            bv[p] = 1;
        }
        doc_bounds_ = bv_rs<sdsl::sd_vector<>>(bv);
        doc_bounds_.init_rs();
        doc_names_ = namev;
    }

    private:

    inline uint64_t doc_bounds_rank(uint64_t i) const {
        return i+1 > doc_bounds_.size() ? doc_bounds_.rank(doc_bounds_.size()) : doc_bounds_.rank(i+1);;
    }

    std::vector<std::string> doc_names_;
    bv_rs<sdsl::sd_vector<>> doc_bounds_;
};

#endif
