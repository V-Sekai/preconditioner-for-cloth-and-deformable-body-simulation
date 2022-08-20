/////////////////////////////////////////////////////////////////////////////////////////////
//  Copyright (C) 2002 - 2022,
//  All rights reserved.
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions
//  are met:
//     1. Redistributions of source code must retain the above copyright
//        notice, this list of conditions and the following disclaimer.
//     2. Redistributions in binary form must reproduce the above copyright
//        notice, this list of conditions and the following disclaimer in the
//        documentation and/or other materials provided with the distribution.
//     3. The names of its contributors may not be used to endorse or promote
//        products derived from this software without specific prior written
//        permission.
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//  A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
//  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//	NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//	SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#pragma once

#include <memory>
#include <vector>
#include <assert.h>

/*************************************************************************
***************************    Min_Max_Swap    ***************************
*************************************************************************/

#define SE_MIN(a, b)		(((a) < (b)) ? (a) : (b))
#define SE_MAX(a, b)		(((a) > (b)) ? (a) : (b))
#define SE_SWAP(a, b)		{ auto c = a; a = b; b = c;}

/*************************************************************************
****************************    Namespace    *****************************
*************************************************************************/

#define SE_NAMESPACE SE
#define SE_NAMESPACE_BEGIN namespace SE_NAMESPACE {
#define SE_USING_NAMESPACE using namespace SE_NAMESPACE; 
#define SE_NAMESPACE_END }

/*************************************************************************
***************************    Noncopyable    ****************************
*************************************************************************/

#define SE_NONCOPYABLE(className)										\
																		\
	className(const className&) = delete;								\
																		\
	void operator=(const className&) = delete;							\

/*************************************************************************
****************************    DLL_Export    ****************************
*************************************************************************/

	
#define SE_INLINE					__inline
#define SE_ALIGN(n)					__declspec(align(n))		//!	Only for windows platform.
#define SE_ASSERT(expression)		assert(expression)


/*************************************************************************
*****************************    Casting    ******************************
*************************************************************************/

#define SE_SCI(ClassObject)			static_cast<int>(ClassObject)
#define SE_SCF(ClassObject)			static_cast<float>(ClassObject)
#define SE_SCD(ClassObject)			static_cast<double>(ClassObject)
#define SE_SCU(ClassObject)			static_cast<unsigned int>(ClassObject)

#ifndef NOMINMAX
	#define NOMINMAX
#endif