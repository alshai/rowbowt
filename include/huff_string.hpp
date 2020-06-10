/*
 * huff_string.hpp
 *
 *  Created on: May 18, 2015
 *      Author: nicola
 *
 *  Huffman-compressed string with access/rank/select. The class is a wrapper on sdsl::wt_huff, with a simpler constructor
 */

#ifndef HUFF_STRING_HPP_
#define HUFF_STRING_HPP_

#include <string>
#include <sdsl/wavelet_trees.hpp>

namespace ri {

class huff_string {

    public:

    huff_string() {}

    huff_string(std::string &s) {
        s.push_back(0);
        sdsl::construct_im(wt, s.c_str(), 1);
        assert(wt.size()==s.size()-1);
    }

    uint8_t operator[](uint64_t i) const {
        assert(i<wt.size());
        return wt[i];
    }

    uint64_t size() const{
        return wt.size();
    }

    uint64_t rank(uint64_t i, uint8_t c) const {
        assert(i<=wt.size());
        return wt.rank(i,c);
    }

    /*
     * position of i-th character c. i starts from 0!
     */
    uint64_t select(uint64_t i, uint8_t c) const {
        return wt.select(i+1,c);
    }

    /* serialize the structure to the ostream
     * \param out     the ostream
     */
    uint64_t serialize(std::ostream& out) const {
        return wt.serialize(out);
    }

    /* load the structure from the istream
     * \param in the istream
     */
    void load(std::istream& in) {
        wt.load(in);
    }

    private:
    sdsl::wt_huff<> wt;
};

}

#endif /* HUFF_STRING_HPP_ */
