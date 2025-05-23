
/*============================================================================

This C source file is part of the SoftPosit Posit Arithmetic Package
by S. H. Leong (Cerlane).

Copyright 2017 2018 A*STAR.  All rights reserved.

This C source file was based on SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3d, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014, 2015, 2016 The Regents of the University of
California.  All Rights Reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice,
    this list of conditions, and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions, and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

 3. Neither the name of the University nor the names of its contributors may
    be used to endorse or promote products derived from this software without
    specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS "AS IS", AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, ARE
DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=============================================================================*/

#include <stdint.h>

#include "platform.hpp"
#include "internals.hpp"

posit_2_t ui64_to_pX2 ( uint64_t a, int x ) {
	int_fast8_t k, log2 = 63;//length of bit (e.g. 18445618173802708991) in int (64 but because we have only 64 bits, so one bit off to accommodate that fact)
	union ui32_pX2 uZ;
	uint_fast64_t uiA=0;
	uint_fast64_t mask = 0x8000000000000000, frac64A;
	uint_fast32_t expA;

	//NaR
	if (a == 0x8000000000000000 || x<2 || x>32)
		uiA = 0x80000000;
	else if (x==2){
		if (a>0) uiA=0x40000000;
	}
	else if ( a > 0xFFFBFFFFFFFFFFFF){//18445618173802708991
		uiA = 0x7FFFC000; // 18446744073709552000
		if (x<18)  uiA&=((int32_t)0x80000000>>(x-1));
	}
	else if ( a < 0x2 )
		uiA = (a << 30);
	else {
		frac64A = a;
		while ( !(frac64A & mask) ) {
			log2--;
			frac64A <<= 1;
		}

		k = (log2 >> 2);

		expA = (log2 & 0x3) ;
		frac64A = (frac64A ^ mask);

		if(k>=(x-2)){//maxpos
			uiA = 0x7FFFFFFF & ((int32_t)0x80000000>>(x-1));
		}
		else if (k==(x-3)){//bitNPlusOne-> first exp bit //bitLast is zero
			uiA = (0x7FFFFFFF ^ (0x3FFFFFFF >> k));
			if( (expA & 0x2) && ((expA&0x1) | frac64A) ) //bitNPlusOne //bitsMore
				 uiA |= ((uint32_t)0x80000000>>(x-1));
		}
		else if (k==(x-4)){
			uiA = (0x7FFFFFFF ^ (0x3FFFFFFF >> k)) | ((expA &0x2)<< (27 - k));
			if(expA&0x1){
				if( (((uint32_t)0x80000000>>(x-1)) & uiA)|| frac64A)
						uiA  += ((uint32_t)0x80000000>>(x-1));
			}

		}
		else if (k==(x-5)){
			uiA = (0x7FFFFFFF ^ (0x3FFFFFFF >> k)) | (expA<< (27 - k));
			mask = (uint64_t)0x800000000 << (k + 32-x);
			if (mask & frac64A){ //bitNPlusOne
				if (((mask - 1) & frac64A) | (expA&0x1)) {
					uiA+= ((uint32_t)0x80000000>>(x-1));
				}
			}
		}
		else{
			uiA = (0x7FFFFFFF ^ (0x3FFFFFFF >> k)) | (expA<< (27 - k)) | ((frac64A>>(k+36)) & ((int32_t)0x80000000>>(x-1)));
			mask = (uint64_t)0x800000000 << (k + 32-x);  //bitNPlusOne position
			if (mask & frac64A) {
				if (((mask - 1) & frac64A) | ((mask << 1) & frac64A)) {
					uiA+= ((uint32_t)0x80000000>>(x-1));
				}
			}
		}
	}
	uZ.ui = uiA;
	return uZ.p;
}
