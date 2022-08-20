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

#include "SeVector.h"
#include "SeMatrix.h"

SE_NAMESPACE_BEGIN

#define SE_E		2.71828182845904523536
#define SE_PI		3.14159265358979323846
#define SE_2_PI		6.28318530717958647692

namespace Math
{
	/*********************************************************************
	****************************    Helper    ****************************
	*********************************************************************/

	template<typename Type> struct Helper
	{
		using Vec2 = SeVector2<Type>;
		using Vec3 = SeVector3<Type>;
		using Vec4 = SeVector4<Type>;

		template<typename... Types> using Function = Type(*)(Types...);

		template<Function<Type> Fn1> static SE_INLINE Vec2 Expand(const Vec2 & a) { return Vec2(Fn1(a.x), Fn1(a.y)); }
		template<Function<Type> Fn1> static SE_INLINE Vec3 Expand(const Vec3 & a) { return Vec3(Fn1(a.x), Fn1(a.y), Fn1(a.z)); }
		template<Function<Type> Fn1> static SE_INLINE Vec4 Expand(const Vec4 & a) { return Vec4(Fn1(a.x), Fn1(a.y), Fn1(a.z), Fn1(a.w)); }

		template<Function<Type, Type> Fn2> static SE_INLINE Vec2 Expand(const Vec2 & a, Type s) { return Vec2(Fn2(a.x, s), Fn2(a.y, s)); }
		template<Function<Type, Type> Fn2> static SE_INLINE Vec3 Expand(const Vec3 & a, Type s) { return Vec3(Fn2(a.x, s), Fn2(a.y, s), Fn2(a.z, s)); }
		template<Function<Type, Type> Fn2> static SE_INLINE Vec4 Expand(const Vec4 & a, Type s) { return Vec4(Fn2(a.x, s), Fn2(a.y, s), Fn2(a.z, s), Fn2(a.w, s)); }

		template<Function<Type, Type> Fn2> static SE_INLINE Vec2 Expand(const Vec2 & a, const Vec2 & b) { return Vec2(Fn2(a.x, b.x), Fn2(a.y, b.y)); }
		template<Function<Type, Type> Fn2> static SE_INLINE Vec3 Expand(const Vec3 & a, const Vec3 & b) { return Vec3(Fn2(a.x, b.x), Fn2(a.y, b.y), Fn2(a.z, b.z)); }
		template<Function<Type, Type> Fn2> static SE_INLINE Vec4 Expand(const Vec4 & a, const Vec4 & b) { return Vec4(Fn2(a.x, b.x), Fn2(a.y, b.y), Fn2(a.z, b.z), Fn2(a.w, b.w)); }

		template<Function<Type, Type, Type> Fn3> static SE_INLINE Vec2 Expand(const Vec2 & a, Type s0, Type s1) { return Vec2(Fn3(a.x, s0, s1), Fn3(a.y, s0, s1)); }
		template<Function<Type, Type, Type> Fn3> static SE_INLINE Vec3 Expand(const Vec3 & a, Type s0, Type s1) { return Vec3(Fn3(a.x, s0, s1), Fn3(a.y, s0, s1), Fn3(a.z, s0, s1)); }
		template<Function<Type, Type, Type> Fn3> static SE_INLINE Vec4 Expand(const Vec4 & a, Type s0, Type s1) { return Vec4(Fn3(a.x, s0, s1), Fn3(a.y, s0, s1), Fn3(a.z, s0, s1), Fn3(a.w, s0, s1)); }

		template<Function<Type, Type, Type> Fn3> static SE_INLINE Vec2 Expand(const Vec2 & a, const Vec2 & b, Type scalar) { return Vec2(Fn3(a.x, b.x, scalar), Fn3(a.y, b.y, scalar)); }
		template<Function<Type, Type, Type> Fn3> static SE_INLINE Vec3 Expand(const Vec3 & a, const Vec3 & b, Type scalar) { return Vec3(Fn3(a.x, b.x, scalar), Fn3(a.y, b.y, scalar), Fn3(a.z, b.z, scalar)); }
		template<Function<Type, Type, Type> Fn3> static SE_INLINE Vec4 Expand(const Vec4 & a, const Vec4 & b, Type scalar) { return Vec4(Fn3(a.x, b.x, scalar), Fn3(a.y, b.y, scalar), Fn3(a.z, b.z, scalar), Fn3(a.w, b.w, scalar)); }
	};

	/*********************************************************************
	*************************    Mathematics    **************************
	*********************************************************************/
	
	SE_INLINE float Log(float x) { return logf(x);}
	SE_INLINE float Abs(float a) { return fabs(a); }
	SE_INLINE float Sin(float a) { return sinf(a); }
	SE_INLINE float Cos(float a) { return cosf(a); }
	SE_INLINE float Tan(float a) { return tanf(a); }
	SE_INLINE float Exp(float a) { return expf(a); }
	SE_INLINE float Asin(float a) { return asinf(a); }
	SE_INLINE float Acos(float a) { return acosf(a); }
	SE_INLINE float Atan(float a) { return atanf(a); }
	SE_INLINE float Sqrt(float a) { return sqrtf(a); }
	SE_INLINE float Cbrt(float x) { return cbrtf(x); }
	SE_INLINE float Ceil(float a) { return ceilf(a); }
	SE_INLINE float Recip(float a) { return 1.0f / a; }
	SE_INLINE float Floor(float a) { return floorf(a); }
	SE_INLINE float Round(float a) { return roundf(a); }
	SE_INLINE float Pow(float x, float y) { return powf(x, y); }
	SE_INLINE float Atan2(float y, float x) { return atan2f(y, x); }
	SE_INLINE constexpr float Radians(float a) { return a * SE_SCF(SE_PI / 180.0); }

	SE_INLINE float Rsqrt(float a) { return 1.0f / sqrtf(a); }

	template<typename Type> SE_INLINE constexpr Type Square(Type a) { return a * a; }
	template<typename Type> SE_INLINE constexpr Type Cubic(Type a) { return a * a * a; }
	template<typename Type> SE_INLINE constexpr Type Min(Type a, Type b) { return SE_MIN(a, b); }
	template<typename Type> SE_INLINE constexpr Type Max(Type a, Type b) { return SE_MAX(a, b); }
	template<typename Type> SE_INLINE constexpr Type Sign(Type a) { return Type(a > 0) - Type(a < 0); }
	template<typename Type> SE_INLINE constexpr Type Clamp(Type a, Type minVal, Type maxVal) { return Min(Max(minVal, a), maxVal); }

	template<typename Type> SE_INLINE bool IsNan(Type a) { return false; }
	template<typename Type> SE_INLINE bool IsInf(Type a) { return false; }
	template<typename Type> SE_INLINE Type Sum(const SeVector2<Type> & a) { return a.x + a.y; }
	template<typename Type> SE_INLINE Type Sum(const SeVector3<Type> & a) { return a.x + a.y + a.z; }
	template<typename Type> SE_INLINE Type Sum(const SeVector4<Type> & a) { return a.x + a.y + a.z + a.w; }
	template<typename Type> SE_INLINE Type Average(const Type & v0, const Type & v1) { return (v0 + v1) / 2; }
	template<typename Type> SE_INLINE Type Average(const Type & v0, const Type & v1, const Type & v2) { return (v0 + v1 + v2) / 3; }
	template<typename Type> SE_INLINE Type Lerp(const Type & a, const Type & b, float ratio) { return a * (1.0f - ratio) + b * ratio; }

	template<typename Type, int Desc> SE_INLINE Type Dot(const SeVector<Type, Desc> & a, const SeVector<Type, Desc> & b) { return Sum(a * b); }
	template<typename Type, int Desc> SE_INLINE SeVector<Type, Desc> Abs(const SeVector<Type, Desc> & a) { return Helper<Type>::Expand<Abs>(a); }
	template<typename Type, int Desc> SE_INLINE SeVector<Type, Desc> Sign(const SeVector<Type, Desc> & a) { return Helper<Type>::Expand<Sign>(a); }
	template<typename Type, int Desc> SE_INLINE SeVector<Type, Desc> Recip(const SeVector<Type, Desc> & a) { return Helper<Type>::Expand<Recip>(a); }
	template<typename Type, int Desc> SE_INLINE SeVector<Type, Desc> Rsqrt(const SeVector<Type, Desc> & a) { return Helper<Type>::Expand<Rsqrt>(a); }
	template<typename Type, int Desc> SE_INLINE SeVector<Type, Desc> Cubic(const SeVector<Type, Desc> & a) { return Helper<Type>::Expand<Cubic>(a); }
	template<typename Type, int Desc> SE_INLINE SeVector<Type, Desc> Square(const SeVector<Type, Desc> & a) { return Helper<Type>::Expand<Square>(a); }
	template<typename Type, int Desc> SE_INLINE SeVector<Type, Desc> Min(const SeVector<Type, Desc> & a, const SeVector<Type, Desc> & b) { return Helper<Type>::template Expand<Min<Type>>(a, b); }
	template<typename Type, int Desc> SE_INLINE SeVector<Type, Desc> Max(const SeVector<Type, Desc> & a, const SeVector<Type, Desc> & b) { return  Helper<Type>::template Expand<Max<Type>>(a, b); }
	template<typename Type, int Desc> SE_INLINE SeVector<Type, Desc> Clamp(const SeVector<Type, Desc> & a, Type minVal, Type maxVal)  { return Helper<Type>::Expand<Clamp>(a, minVal, maxVal); }
	
	template<typename Type, typename... Args> SE_INLINE Type Min(const Type & a, const Type & b, const Args &... args) { return Min(Min(a, b), args...); }
	template<typename Type, typename... Args> SE_INLINE Type Max(const Type & a, const Type & b, const Args &... args) { return Max(Max(a, b), args...); }

	SE_INLINE SeVector3<float> Cross(const SeVector2<float> & a, const SeVector2<float> & b) { return SeVector3<float>(0.0f, 0.0f, a.x * b.y - a.y * b.x); }
	SE_INLINE SeVector3<float> Cross(const SeVector3<float> & a, const SeVector3<float> & b) { return SeVector3<float>(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x); }
	SE_INLINE SeVector4<float> Cross(const SeVector4<float> & a, const SeVector4<float> & b) { return SeVector4<float>(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x, 0.0f); }

	template<int Desc> SE_INLINE float SqrLength(const SeVector<float, Desc> & a) { return Dot(a, a); }
	template<int Desc> SE_INLINE float Length(const SeVector<float, Desc> & a) { return Sqrt(Dot(a, a)); }
	template<int Desc> SE_INLINE float RelativeDistance(const SeVector<float, Desc>& a, const SeVector<float, Desc>& b) { float maxLength = Max(Length(b), Length(a)); if (maxLength == 0.f) return 0.f; else return Length(b - a) / maxLength; }
	template<int Desc> SE_INLINE SeVector<float, Desc> Sin(const SeVector<float, Desc> & a) { return Helper<float>::Expand<Sin>(a); }
	template<int Desc> SE_INLINE SeVector<float, Desc> Cos(const SeVector<float, Desc> & a) { return Helper<float>::Expand<Cos>(a); }
	template<int Desc> SE_INLINE SeVector<float, Desc> Tan(const SeVector<float, Desc> & a) { return Helper<float>::Expand<Tan>(a); }
	template<int Desc> SE_INLINE SeVector<float, Desc> Asin(const SeVector<float, Desc> & a) { return Helper<float>::Expand<Asin>(a); }
	template<int Desc> SE_INLINE SeVector<float, Desc> Acos(const SeVector<float, Desc> & a) { return Helper<float>::Expand<Acos>(a); }
	template<int Desc> SE_INLINE SeVector<float, Desc> Atan(const SeVector<float, Desc> & a) { return Helper<float>::Expand<Atan>(a); }
	template<int Desc> SE_INLINE SeVector<float, Desc> Sqrt(const SeVector<float, Desc> & a) { return Helper<float>::Expand<Sqrt>(a); }
	template<int Desc> SE_INLINE SeVector<float, Desc> Ceil(const SeVector<float, Desc> & a) { return Helper<float>::Expand<Ceil>(a); }
	template<int Desc> SE_INLINE SeVector<float, Desc> Floor(const SeVector<float, Desc> & a) { return Helper<float>::Expand<Floor>(a); }
	template<int Desc> SE_INLINE SeVector<float, Desc> Round(const SeVector<float, Desc> & a) { return Helper<float>::Expand<Round>(a); }
	template<int Desc> SE_INLINE SeVector<float, Desc> Radians(const SeVector<float, Desc> & a) { return Helper<float>::Expand<Radians>(a); }
	template<int Desc> SE_INLINE float TripleProduct(const SeVector<float, Desc> & a, const SeVector<float, Desc> & b, const SeVector<float, Desc> & c) { return Dot(a, Cross(b, c)); }

	//!	若输入a为NaN或0，返回0.0f。
	template<int Desc> SE_INLINE SeVector<float, Desc> Normalize(const SeVector<float, Desc> & a)
	{
		float sqrLen = SqrLength(a);		if (sqrLen > 0.0f) { return a * Rsqrt(sqrLen); }
		
	//	printf("Normalizing an invalid vector: (%f, %f, %f)!\n", a.x, a.y, a.z);

		return SeVector<float, Desc>(0.0f);
	}

	template<int Desc> SE_INLINE bool Normalized(SeVector<float, Desc> & a)
	{
		float sqrLen = SqrLength(a);
		if (sqrLen > 0.f) 
		{ 
			a *= Rsqrt(sqrLen); 
			return true;
		}
		return false;
	}
	
	template<> SE_INLINE bool IsNan(float a) { return isnan(a); }
	template<> SE_INLINE bool IsNan(double a) { return isnan(a); }
	template<> SE_INLINE bool IsNan(SeVector2<float> a) { return isnan(a.x) || isnan(a.y); }
	template<> SE_INLINE bool IsNan(SeVector3<float> a) { return isnan(a.x) || isnan(a.y) || isnan(a.z); }
	template<> SE_INLINE bool IsNan(SeVector4<float> a) { return isnan(a.x) || isnan(a.y) || isnan(a.z) || isnan(a.w); }
	template<> SE_INLINE bool IsNan(SeMatrix3f a) { return IsNan(a.Column(0)) || IsNan(a.Column(1)) || IsNan(a.Column(2)); }

	template<> SE_INLINE bool IsInf(float a) { return isinf(a); }
	template<> SE_INLINE bool IsInf(double a) { return isinf(a); }
	template<> SE_INLINE bool IsInf(SeVector2<float> a) { return isinf(a.x) || isinf(a.y); }
	template<> SE_INLINE bool IsInf(SeVector3<float> a) { return isinf(a.x) || isinf(a.y) || isinf(a.z); }
	template<> SE_INLINE bool IsInf(SeVector4<float> a) { return isinf(a.x) || isinf(a.y) || isinf(a.z) || isinf(a.w); }
	template<> SE_INLINE bool IsInf(SeMatrix3f a) { return IsInf(a.Column(0)) || IsInf(a.Column(1)) || IsInf(a.Column(2)); }
};

SE_NAMESPACE_END

#undef min
#undef max