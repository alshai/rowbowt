/*
 * fbb_string.hpp
 *
 *  Created on: May 18, 2015
 *      Author: nicola
 *
 *  Huffman-compressed string with access/rank/select. The class is a wrapper on sdsl::wt_fbb, with a simpler constructor
 */

#ifndef FBB_STRING_HPP_
#define FBB_STRING_HPP_

#include <string>
#include "sdsl/int_vector_buffer.hpp"
#include "wt_fbb.hpp"

namespace ri {

class fbb_string {

    public:

    fbb_string() {}

    fbb_string(std::string fname) {
        std::ifstream ifs(fname, std::ios_base::ate);
        if (!ifs.good()) {
            std::cerr << "cannot read " << fname << std::endl;
            exit(1);
        }
        size_t len = ifs.tellg();
        ifs.close();
        auto sbuf = sdsl::int_vector_buffer<8>(fname, std::ios::in, 1024*1024, 8, true);
        wt = wt_fbb(sbuf, len);
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
    /* TODO: SELECT NOT CURRENTLY SUPPORTED
     * uint64_t select(uint64_t i, uint8_t c) const {
     *    return wt.select(i+1,c);
     * }
     */

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
    wt_fbb<> wt;
};

}

#endif /* FBB_STRING_HPP_ */
