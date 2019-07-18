/*
 *	Author: Yun Wu
 *	Created by: 2019-06-13
 *	Copyright @ Yun Wu
 *
 */

#ifndef SRC_DATA_HPP_
#define SRC_DATA_HPP_

// Data structure
#define ROW 16 // matrix row number
#define COL 16 // matrix column number
#define DIAG (COL<ROW ? COL : ROW) // diagonal matrix size
#define DIAG_VALUE 50 // diagonal matrix value scale
#define DIAG_RATIO 0.5 // diagonal matrix sparse ratio

#define SIZE (COL>ROW ? COL : ROW) // square matrix size

#define ROW1 230 // left multiply matrix row number
#define COL1 220 // left multiply matrix column number
#define ROW2 COL1 // right multiply matrix row number
#define COL2 240 // right multiply matrix column number

// Data type
#define FLOAT_SIZE 1000 // floating point fraction scale
#define INTEGER_SCALE 1 // floating point integer scale

#endif /* SRC_DATA_HPP_ */
