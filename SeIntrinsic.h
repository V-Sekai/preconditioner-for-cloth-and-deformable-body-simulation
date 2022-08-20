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

#include <atomic>
#include <stdint.h>
#include <intrin.h>

#include "SeMatrix.h"

SE_NAMESPACE_BEGIN

namespace Intrinsic
{
	/*************************************************************************
	**************************    Count of bits    ***************************
	*************************************************************************/

#ifdef WIN32
    inline int Clz32(unsigned int Val) { return __lzcnt(Val); }
    inline int Clz64(unsigned __int64 Val) { return __lzcnt64(Val); }
    inline int Popcount32(unsigned int Val) { return __popcnt(Val); }
    inline int Popcount64(unsigned __int64 Val) { return __popcnt64(Val); }

    #pragma intrinsic(_BitScanForward)
    inline int Ffs(int s)
    {
        unsigned long index = -1;
        unsigned char isNonzero = _BitScanForward(&index, s);
        /*s should be greater than 0; index is zero-based*/
        return isNonzero ? (index + 1) : 0;
    }

#else
    inline int Clz32(unsigned int Val) { return std::countl_zero(Val); }
    inline int Clz64( uint64_t Val) { return std::countl_zero(Val); }
    inline int Popcount32(unsigned int Val) { return std::popcount(Val); }
    inline int Popcount64(uint64_t Val) { return std::popcount(Val); }
    inline int Ffs(int s) { return ffs(s); }
#endif

	/*************************************************************************
	************************    Atomic operations    *************************
	*************************************************************************/


	template<typename Type> inline Type AtomicOr(Type * Address, Type Val) { return reinterpret_cast<std::atomic<Type>*>(Address)->fetch_or(Val,std::memory_order_relaxed); }
	template<typename Type> inline Type AtomicXor(Type * Address, Type Val) { return reinterpret_cast<std::atomic<Type>*>(Address)->fetch_xor(Val,std::memory_order_relaxed); }
	template<typename Type> inline Type AtomicAnd(Type * Address, Type Val) { return reinterpret_cast<std::atomic<Type>*>(Address)->fetch_and(Val,std::memory_order_relaxed); }
	template<typename Type> inline Type AtomicAdd(Type * Address, Type Val) { return reinterpret_cast<std::atomic<Type>*>(Address)->fetch_add(Val,std::memory_order_relaxed); }
	template<typename Type> inline Type AtomicSub(Type * Address, Type Val) { return reinterpret_cast<std::atomic<Type>*>(Address)->fetch_sub(Val,std::memory_order_relaxed); }
	template<typename Type> inline Type AtomicExch(Type * Address, Type Val) { return reinterpret_cast<std::atomic<Type>*>(Address)->exchange(Val,std::memory_order_relaxed); }
	template<typename Type> inline Type AtomicCAS(Type* Address, Type Exp, Type Val) { reinterpret_cast<std::atomic<Type>*>(Address)->compare_exchange_strong(Exp, Val, std::memory_order_relaxed); return Exp; }

	template< typename Type, typename Op > Type AtomicOperationByCas(Type* addr, Type value, Op op)
	{
		Type oldMyValue = *addr;
		bool isChangedByOtherThread = true;
		while (isChangedByOtherThread)
		{
			oldMyValue = *addr;
			Type m = op(value, oldMyValue);
			float newMyValue = AtomicCAS<Type>(addr, oldMyValue, m);
			isChangedByOtherThread = oldMyValue != newMyValue;
		}
		return oldMyValue;
	}


	template<typename Type> inline Type AtomicMax(Type* addr, Type value) 
	{ 
		return AtomicOperationByCas<Type>(addr, value, [](const auto& a, const auto& b) { return a > b ? a : b; });
	}

	template<typename Type> inline Type AtomicMin(Type* addr, Type value) 
	{ 
		return AtomicOperationByCas<Type>(addr, value, [](const auto& a, const auto& b) { return a < b ? a : b; });
	}

#if defined(WIN32) && (_MSVC_LANG > 201703L) // only for C++ 20 on Windows
	inline void AtomicAdd(SeMatrix3f* address, SeMatrix3f value)
	{
		for (int i = 0; i < 3; i++)
        {
			for (int j = 0; j < 3; j++)
            {
                    AtomicAdd(&(*address)(i, j), value(i, j));
            }
        }
	}

	template<typename Type, int Desc> inline void AtomicAdd(SeVector<Type, Desc> * address, SeVector<Type, Desc> value)
	{
		for (int i = 0; i < SeVector<Type, Desc>::Elements; i++)
        {
            AtomicAdd(&(*address)[i], value[i]);
        }
	}
#else
    inline void AtomicAdd(SeMatrix3f* address, SeMatrix3f value)
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
            AtomicOperationByCas<float>( &(*address)(i, j), value(i, j), [](const auto& a, const auto& b) { return a + b; });
            }
        }
    }

    template<typename Type, int Desc> inline void AtomicAdd(SeVector<Type, Desc> * address, SeVector<Type, Desc> value)
    {
        for (int i = 0; i < SeVector<Type, Desc>::Elements; i++)
        {
            AtomicOperationByCas<float>(&(*address)[i], value[i], [](const auto& a, const auto& b) { return a + b; });
        }
    }

#endif

	/*************************************************************************
	*************************    Warp operations    **************************
	*************************************************************************/


	 inline void SyncWarp() {}
	 inline void SyncBlock() {}
	 inline unsigned int ActiveMask();
	 inline unsigned int LanemaskLt();
	 inline unsigned int LanemaskLt(unsigned int laneId) { return  (1U << laneId) - 1; }
}

SE_NAMESPACE_END