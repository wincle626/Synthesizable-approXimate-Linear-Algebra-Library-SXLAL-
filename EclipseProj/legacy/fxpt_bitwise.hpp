/*
 *	This is light fixed point arithmetic
 *	Author: Yun Wu
 *	Created by: 2019-06-13
 *	Copyright @ Yun Wu
 *
 */

#ifndef SRC_FXPT_BITWISE_HPP_
#define SRC_FXPT_BITWISE_HPP_

#include "common.hpp"
#include "data.hpp"

#define BITWIDTH 8 // fixed point total bit length
#define FRACWIDTH 4 // fixed point fraction bit length

typedef struct{
	const size_t bitwidth;
	boost::dynamic_bitset<> data;
}BIN_DYNBITSET;

typedef struct{
	const size_t bitwidth;
	const size_t fracwidth;
	boost::dynamic_bitset<> data;
}FXPT_DYNBITSET;

typedef struct{
	const size_t bitwidth;
	std::bitset<BITWIDTH> data;
}BIN_BITSET;

typedef struct{
	const size_t bitwidth;
	const size_t fracwidth;
	std::bitset<BITWIDTH> data;
}FXPT_BITSET;

typedef bool BINARRAY[BITWIDTH] ;
typedef bool FXPT_INTE_ARRAY[BITWIDTH-FRACWIDTH];
typedef bool FXPT_FRAC_ARRAY[FRACWIDTH];

#endif /* SRC_FXPT_BITWISE_HPP_ */
