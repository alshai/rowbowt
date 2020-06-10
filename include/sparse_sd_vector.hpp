 /*
  * sparse_sd_vector: a wrapper on sdsl::sd_vector of the sdsl library, with support for rank/select1
  */

//============================================================================


#ifndef INTERNAL_SPARSE_SD_VECTOR_HPP_
#define INTERNAL_SPARSE_SD_VECTOR_HPP_

#include <vector>
#include "sdsl/bit_vectors.hpp"


namespace ri{

class sparse_sd_vector{

    public:

    sparse_sd_vector() {}

    sparse_sd_vector(std::vector<bool> &b) {
        if(b.size()==0) return;
        u = b.size();
        sdsl::bit_vector bv(b.size());
        for(uint64_t i=0;i<b.size();++i)
            bv[i] = b[i];
        sdv = sdsl::sd_vector<>(bv);
        rank1 = sdsl::sd_vector<>::rank_1_type(&sdv);
        select1 = sdsl::sd_vector<>::select_1_type(&sdv);
    }

    sparse_sd_vector(const sdsl::bit_vector& bv) {
        sdv = sdsl::sd_vector<>(std::move(bv));
        rank1 = sdsl::sd_vector<>::rank_1_type(&sdv);
        select1 = sdsl::sd_vector<>::select_1_type(&sdv);
    }

    sparse_sd_vector(const sparse_sd_vector& other) : sdv(other.sdv) {
        u = other.sdv.size();
        rank1 = sdsl::sd_vector<>::rank_1_type(&sdv);
        select1 = sdsl::sd_vector<>::select_1_type(&sdv);
    }

    sparse_sd_vector(sparse_sd_vector&& other) : sdv(std::move(other.sdv)) {
        u = other.sdv.size();
        rank1 = sdsl::sd_vector<>::rank_1_type(&sdv);
        select1 = sdsl::sd_vector<>::select_1_type(&sdv);
    }

    sparse_sd_vector& operator=(const sparse_sd_vector& other) {
        u = other.sdv.size();
        sdv = sdsl::sd_vector<>(other.sdv);
        rank1 = sdsl::sd_vector<>::rank_1_type(&sdv);
        select1 = sdsl::sd_vector<>::select_1_type(&sdv);
        return *this;
    }

    sparse_sd_vector& operator=(sparse_sd_vector&& other) {
        u = other.sdv.size();
        sdv = std::move(other.sdv);
        rank1 = sdsl::sd_vector<>::rank_1_type(&sdv);
        select1 = sdsl::sd_vector<>::select_1_type(&sdv);
        return *this;
    }

    sparse_sd_vector& operator=(const std::vector<bool>& b) {
        if(b.size()==0) return *this;
        u = b.size();
        sdsl::bit_vector bv(b.size());
        for(uint64_t i=0;i<b.size();++i)
            bv[i] = b[i];
        sdv = sdsl::sd_vector<>(bv);
        rank1 = sdsl::sd_vector<>::rank_1_type(&sdv);
        select1 = sdsl::sd_vector<>::select_1_type(&sdv);
        return *this;
    }

    sparse_sd_vector& operator=(std::vector<bool>&& b) {
        if(b.size()==0) return *this;
        u = b.size();
        sdsl::bit_vector bv(b.size());
        for(uint64_t i=0;i<b.size();++i)
            bv[i] = b[i];
        sdv = sdsl::sd_vector<>(bv);
        rank1 = sdsl::sd_vector<>::rank_1_type(&sdv);
        select1 = sdsl::sd_vector<>::select_1_type(&sdv);
        return *this;
    }

    /*
     * argument: position i in the bitvector
     * returns: bit in position i
     * only access! the bitvector is static.
     */
    bool operator[](uint64_t i) const {
        assert(i<size());
        return sdv[i];
    }

    bool at(uint64_t i) const {
        return operator[](i);
    }

    /*
     * argument: position i in the bitvector
     * returns: number of bits equal to 1 before position i excluded
     */
    uint64_t rank(uint64_t i) const {
        assert(i<=size());
        return rank1(i);
    }

    /*
     * input: position 0<=i<=n
     * output: predecessor of i (position i excluded) in
     * rank space (apply select to get bitvector space)
     */
    uint64_t predecessor_rank(uint64_t i) const {
        assert(rank(i)>0);
        return rank(i)-1;
    }

    /*
     * input: position 0<=i<=n
     * output: predecessor of i (i excluded) in
     * bitvector space
     */
    uint64_t predecessor(uint64_t i) const {
        assert(rank(i)>0);
        return select(rank(i)-1);
    }

    /*
     * input: position 0<=i<=n
     * output: rank of predecessor of i (i excluded) in
     * bitvector space. If i does not have a predecessor,
     * return rank of the last bit set in the bitvector
     */
    uint64_t predecessor_rank_circular(uint64_t i) const {
        return rank(i)==0 ? number_of_1()-1 : rank(i)-1;
    }

    /*
     * retrieve length of the i-th gap (i>=0). gap length includes the leading 1
     * \param i<number_of_1()
     *
     */
    uint64_t gapAt(uint64_t i) const {
        assert(i<number_of_1());
        if(i==0) return select(0)+1;
        return select(i)-select(i-1);
    }

    /*
     * argument: integer i>0
     * returns: position of the i-th one in the bitvector. i starts from 0!
     */
    uint64_t select(uint64_t i) const {
        assert(i<number_of_1());
        return select1(i+1);//in sdsl::sd_vector, i starts from 1
    }

    /*
    * returns: size of the bitvector
    */
    uint64_t size() const {
        return u;
    }

    /*
     * returns: number of 1s in the bitvector
     */
    uint64_t number_of_1() const {
        return rank1(size());
    }

    /* serialize the structure to the ostream
     * \param out     the ostream
     */
    uint64_t serialize(std::ostream& out) const {
        uint64_t w_bytes = 0;
        out.write((char*)&u, sizeof(u));
        w_bytes += sizeof(u);
        if(u==0) return w_bytes;
        w_bytes += sdv.serialize(out);
        return w_bytes;
    }

    /* load the structure from the istream
     * \param in the istream
     */
    void load(std::istream& in) {
        in.read((char*)&u, sizeof(u));
        if(u==0) return;
        sdv.load(in);
        rank1 = sdsl::sd_vector<>::rank_1_type(&sdv);
        select1 = sdsl::sd_vector<>::select_1_type(&sdv);
    }

    private:

    uint64_t u = 0;
    sdsl::sd_vector<> sdv;
    sdsl::sd_vector<>::rank_1_type rank1;
    sdsl::sd_vector<>::select_1_type select1;

};

}


#endif /* INTERNAL_SPARSE_SD_VECTOR_HPP_ */
