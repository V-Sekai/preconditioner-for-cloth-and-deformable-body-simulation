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

#include "SeVectorSimd.h"

#include "SeMath.h"

SE_NAMESPACE_BEGIN

namespace Math
{
	template<typename Type> struct HelperSimd
	{
		using Smv3 = SimdVector3<Type>;

		template<typename... Types> using Function = Type(*)(Types...);

		template<Function<Type> Fn1> static SE_INLINE				Smv3 Expand(const Smv3& a) { return Smv3(Fn1(a.x), Fn1(a.y), Fn1(a.z)); }
		
		template<Function<Type, Type> Fn2> static SE_INLINE			Smv3 Expand(const Smv3& a, Type s) { return Smv3(Fn2(a.x, s), Fn2(a.y, s), Fn2(a.z, s)); }
		
		template<Function<Type, Type> Fn2> static SE_INLINE			Smv3 Expand(const Smv3& a, const Smv3& b) { return Smv3(Fn2(a.x, b.x), Fn2(a.y, b.y), Fn2(a.z, b.z)); }
		
		template<Function<Type, Type, Type> Fn3> static SE_INLINE	Smv3 Expand(const Smv3& a, Type s0, Type s1) { return Smv3(Fn3(a.x, s0, s1), Fn3(a.y, s0, s1), Fn3(a.z, s0, s1)); }
		
		template<Function<Type, Type, Type> Fn3> static SE_INLINE	Smv3 Expand(const Smv3& a, const Smv3& b, Type scalar) { return Smv3(Fn3(a.x, b.x, scalar), Fn3(a.y, b.y, scalar), Fn3(a.z, b.z, scalar)); }
	};

	template<typename Type> SE_INLINE SimdVector3<Type> Abs(const SimdVector3<Type>& a) { return HelperSimd<Type>::Expand<Abs>(a); }
	template<typename Type> SE_INLINE SimdVector3<Type> Sign(const SimdVector3<Type>& a) { return HelperSimd<Type>::Expand<Sign>(a); }
	template<typename Type> SE_INLINE SimdVector3<Type> Recip(const SimdVector3<Type>& a) { return HelperSimd<Type>::Expand<Recip>(a); }
	template<typename Type> SE_INLINE SimdVector3<Type> Rsqrt(const SimdVector3<Type>& a) { return HelperSimd<Type>::Expand<Rsqrt>(a); }
	template<typename Type> SE_INLINE SimdVector3<Type> Cubic(const SimdVector3<Type>& a) { return HelperSimd<Type>::Expand<Cubic>(a); }
	template<typename Type> SE_INLINE SimdVector3<Type> Square(const SimdVector3<Type>& a) { return HelperSimd<Type>::Expand<Square>(a); }
	template<typename Type> SE_INLINE SimdVector3<Type> Min(const SimdVector3<Type>& a, const SimdVector3<Type>& b) { return HelperSimd<Type>::Expand<Min>(a, b); }
	template<typename Type> SE_INLINE SimdVector3<Type> Max(const SimdVector3<Type>& a, const SimdVector3<Type>& b) { return HelperSimd<Type>::Expand<Max>(a, b); }
	template<typename Type> SE_INLINE SimdVector3<Type> Clamp(const SimdVector3<Type>& a, Type minVal, Type maxVal) { return HelperSimd<Type>::Expand<Clamp>(a, minVal, maxVal); }

	
	template<> SE_INLINE SimdVector3<float> Sin(const SimdVector3<float>& a) { return HelperSimd<float>::Expand<Sin>(a); }
	template<> SE_INLINE SimdVector3<float> Cos(const SimdVector3<float>& a) { return HelperSimd<float>::Expand<Cos>(a); }
	template<> SE_INLINE SimdVector3<float> Tan(const SimdVector3<float>& a) { return HelperSimd<float>::Expand<Tan>(a); }
	template<> SE_INLINE SimdVector3<float> Asin(const SimdVector3<float>& a) { return HelperSimd<float>::Expand<Asin>(a); }
	template<> SE_INLINE SimdVector3<float> Acos(const SimdVector3<float>& a) { return HelperSimd<float>::Expand<Acos>(a); }
	template<> SE_INLINE SimdVector3<float> Atan(const SimdVector3<float>& a) { return HelperSimd<float>::Expand<Atan>(a); }
	template<> SE_INLINE SimdVector3<float> Ceil(const SimdVector3<float>& a) { return HelperSimd<float>::Expand<Ceil>(a); }
	template<> SE_INLINE SimdVector3<float> Floor(const SimdVector3<float>& a) { return HelperSimd<float>::Expand<Floor>(a); }
	template<> SE_INLINE SimdVector3<float> Round(const SimdVector3<float>& a) { return HelperSimd<float>::Expand<Round>(a); }
	template<> SE_INLINE SimdVector3<float> Radians(const SimdVector3<float>& a) { return HelperSimd<float>::Expand<Radians>(a); }

	/*********************************************************************
	*************************    Mathematics    **************************
	*********************************************************************/

	template<typename Type> SE_INLINE Type Sum(const SimdVector3<Type> & a) { return a.x + a.y + a.z; }
	template<typename Type> SE_INLINE Type Dot(const SimdVector3<Type> & a, const SimdVector3<Type> & b) { return Sum(a * b); }
	
	template<> SE_INLINE bool	IsNan(const Simd3f & a) { return isnan(a.x) || isnan(a.y) || isnan(a.z) || isnan(a.w); }
	template<> SE_INLINE bool	IsInf(const Simd3f & a) { return isinf(a.x) || isinf(a.y) || isinf(a.z) || isinf(a.w); }

	template<> SE_INLINE Simd3f Recip(const Simd3f & a) { return Simd3f(_mm_rcp_ps(a.pack)); }
	template<> SE_INLINE Simd3f Sqrt(const Simd3f & a) { return Simd3f(_mm_sqrt_ps(a.pack)); }
	template<> SE_INLINE Simd3f Rsqrt(const Simd3f & a) { return Simd3f(_mm_rsqrt_ps(a.pack)); }
	template<> SE_INLINE Simd3f Min(const Simd3f & a, const Simd3f & b) { return Simd3f(_mm_min_ps(a.pack, b.pack)); }
	template<> SE_INLINE Simd3f Max(const Simd3f & a, const Simd3f & b) { return Simd3f(_mm_max_ps(a.pack, b.pack)); }

	template<> SE_INLINE float  SqrLength(const Simd3f& a) { return Dot(a, a); }
	template<> SE_INLINE float  Length(const Simd3f& a) { return Sqrt(Dot(a, a)); }
	template<> SE_INLINE float  RelativeDistance(const Simd3f& a, const Simd3f& b) { float maxLength = Max(Length(b), Length(a)); if (maxLength == 0.f) return 0.f; else return Length(b - a) / maxLength; }

	template<> SE_INLINE Simd3f Normalize(const Simd3f & a)
	{
		float sqrLength = Math::SqrLength(a);
		return sqrLength > 0.f ? a * Rsqrt(Simd3f(sqrLength, sqrLength, sqrLength, 1.f)) : Simd3f(0.0f);
	}

	template<> SE_INLINE bool Normalized(Simd3f & a)
	{
		float sqrLength = Math::SqrLength(a);
		if (sqrLength > 0.f)
		{
			a *= Rsqrt(Simd3f(sqrLength, sqrLength, sqrLength, 1.f));
			return true;
		}
		a = Simd3f(0.f);
		return false;
	}

	SE_INLINE Simd3f Cross(const Simd3f & a, const Simd3f & b)
	{
		__m128 t0 = _mm_shuffle_ps(a.pack, a.pack, _MM_SHUFFLE(3, 0, 2, 1));	//!	(y0, z0, x0, 0)
		__m128 t1 = _mm_shuffle_ps(b.pack, b.pack, _MM_SHUFFLE(3, 1, 0, 2));	//!	(z1, x1, y1, 0)
		
		__m128 t2 = _mm_shuffle_ps(t0, t0, _MM_SHUFFLE(3, 0, 2, 1));			//!	(z0, x0, y0, 0)
		__m128 t3 = _mm_shuffle_ps(t1, t1, _MM_SHUFFLE(3, 1, 0, 2));			//!	(y1, z1, x1, 0)

		return _mm_sub_ps(_mm_mul_ps(t0, t1), _mm_mul_ps(t2, t3));
	}

	template<> SE_INLINE float TripleProduct(const SimdVector3<float>& a, const SimdVector3<float>& b, const SimdVector3<float>& c)
	{ 
		return Dot(a, Cross(b, c)); 
	}
};

SE_NAMESPACE_END

#undef min
#undef max
