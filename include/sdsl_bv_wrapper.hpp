 /*
  * sdsl_bv_wrapper: an wrapper on biy_vector of the sdsl library, with support for rank and select
  */


#ifndef SDSL_BV_HPP
#define SDSL_BV_HPP

#include <vector>

namespace ri{

template<typename BvType=sdsl::bit_vector>
class sdsl_bv_wrapper : public BvType {

public:

	/*
	 * empty constructor. Initialize bitvector with length 0.
	 */
	sdsl_bv_wrapper() : BvType() {}

    sdsl_bv_wrapper(const sdsl_bv_wrapper& rhs) : BvType(rhs) {
        init_rs();
    }

    sdsl_bv_wrapper(sdsl_bv_wrapper&& rhs) : BvType(rhs) {
        init_rs();
    }

    sdsl_bv_wrapper(const sdsl::bit_vector& rhs) : BvType(rhs) {
        init_rs();
    }

    sdsl_bv_wrapper(sdsl::bit_vector&& rhs) : BvType(rhs) {
        init_rs();
    }

	BvType& operator=(const sdsl_bv_wrapper& rhs) {
        BvType::operator=(rhs);
        init_rs();
	    return *this;
	}

	BvType& operator=(sdsl_bv_wrapper&& rhs) {
        BvType::operator=(rhs);
        init_rs();
	    return *this;
	}

	/*
	 * argument: position i in the bitvector
	 * returns: number of bits equal to 1 before position i excluded
	 */
	uint64_t rank(uint64_t i){
		return rank1(i);
	}

	/*
	 * argument: integer i>=1
	 * returns: position of the i-th one in the bitvector. i starts from 1!
	 */
	uint64_t select(uint64_t i){
		return select1(i);
	}

	/*
	 * returns: number of 1s in the bitvector
	 */
	uint64_t number_of_1(){
        return rank1(BvType::size()); 
    }

	/* load the structure from the istream
	 * \param in the istream
	 */
	void load(std::istream& in) {
        BvType::load(in);
        init_rs();
	}

    void init_rs() {
        sdsl::util::init_support(rank1, this);
        sdsl::util::init_support(select1, this);
    }

private:

	typename BvType::rank_1_type rank1;
	typename BvType::select_1_type select1;

};

}


#endif /* SDSL_BV_HPP */
