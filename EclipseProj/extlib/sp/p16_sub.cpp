
/*============================================================================

This C source file is part of the SoftPosit Posit Arithmetic Package
by S. H. Leong (Cerlane).

Copyright 2017 A*STAR.  All rights reserved.

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3d, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014, 2015, 2016 The Regents of the University of
California.  All rights reserved.

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


#include "../sp/include/internals.hpp"
#include "../sp/include/platform.hpp"


posit16_t p16_sub( posit16_t a, posit16_t b ){

    union ui16_p16 uA, uB;
    uint_fast16_t uiA, uiB;
    union ui16_p16 uZ;

    uA.p = a;
	uiA = uA.ui;
	uB.p = b;
	uiB = uB.ui;

#ifdef SOFTPOSIT_EXACT
		uZ.ui.exact = (uiA.ui.exact & uiB.ui.exact);
#endif

    //infinity
	if ( uiA==0x8000 || uiB==0x8000 ){
#ifdef SOFTPOSIT_EXACT
		uZ.ui.v = 0x8000;
		uZ.ui.exact = 0;
#else
		uZ.ui = 0x8000;
#endif
		return uZ.p;
	}
    //Zero
	else if ( uiA==0 || uiB==0 ){
#ifdef SOFTPOSIT_EXACT
		uZ.ui.v = (uiA | -uiB);
		uZ.ui.exact = 0;
#else
		uZ.ui = (uiA | -uiB);
#endif
		return uZ.p;
	}

	//different signs
	if ((uiA^uiB)>>15)
		return softposit_addMagsP16(uiA, (-uiB & 0xFFFF));
	else
		return softposit_subMagsP16(uiA, (-uiB & 0xFFFF));




}

