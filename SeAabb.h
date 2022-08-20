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

#include "SeMath.h"
#include <cfloat>

SE_NAMESPACE_BEGIN

/*************************************************************************
******************************    SeAabb    ******************************
*************************************************************************/
 
/**
 *	@brief		Axis-align Bounding Box.
 */
template<typename VecType> struct SeAabb
{
	static_assert(std::is_floating_point<typename VecType::value_type>::value, "Only support floating type currently.");

public:

	//!	@brief	Default constructor.
	 SeAabb() : Lower(FLT_MAX), Upper(-FLT_MAX) {}

	//!	@brief	Constructed by a given point.
	 SeAabb(const VecType & Point) : Lower(Point), Upper(Point) {}

	//!	@brief	Constructed by a paired points explicitly.
	 explicit SeAabb(const VecType & P1, const VecType & P2) : Lower(Math::Min(P1, P2)), Upper(Math::Max(P1, P2)) {}

	//!	@brief	Constructed by a paired points explicitly.
	 explicit SeAabb(const VecType & P1, const VecType & P2, const VecType & P3) : Lower(Math::Min(P1, P2, P3)), Upper(Math::Max(P1, P2, P3)) {}

    
    

public:

	 SeAabb<VecType> operator+(const SeAabb<VecType> & Other) const
	{
		return SeAabb<VecType>(Math::Min(Lower, Other.Lower), Math::Max(Upper, Other.Upper));
	}

	 SeAabb<VecType> operator+(const VecType & Point) const
	{
		return SeAabb<VecType>(Math::Min(Lower, Point), Math::Max(Upper, Point));
	}

	 void operator+=(const SeAabb<VecType> & Other)
	{
		Lower = Math::Min(Lower, Other.Lower);		Upper = Math::Max(Upper, Other.Upper);
	}

	 void operator+=(const VecType & Point)
	{
		Lower = Math::Min(Lower, Point);		Upper = Math::Max(Upper, Point);
	}

	 void Enlarge(float Amount)
	{
		Lower -= Amount;		Upper += Amount;
	}

	 VecType Center() const
	{
		return Math::Average(Lower, Upper);
	}

	 VecType Extent() const
	{
		return Upper - Lower;
	}

public:

	VecType Lower, Upper;
};


/*************************************************************************
***********************    IsInside<AABB-AABB>    ************************
*************************************************************************/

template<typename VecT> 
static SE_INLINE bool IsContain(const SeAabb<VecT>& aabb, const VecT& P) 
{ 
	const int componentCount = VecT::Elements;
	for (int i = 0; i < componentCount; i++)
	{
		if (P[i]<aabb.Lower[i] || P[i]>aabb.Upper[i])
		{
			return false;
		}
	}
	return true; 

}
template<typename VecT> 
static SE_INLINE bool IsContain(const SeAabb<VecT>& aabb, const VecT& P,float radium) 
{ 
	auto enlargedAabb = aabb;
	enlargedAabb.Upper += radium;
	enlargedAabb.Lower -= radium;
	return IsContain(aabb, enlargedAabb);
}

template<> 
static SE_INLINE bool IsContain(const SeAabb<Float2>& aabb, const Float2& P) 
{ 
	return !((aabb.Upper.x < P.x) || (aabb.Upper.y < P.y) ||
			 (aabb.Lower.x > P.x) || (aabb.Lower.y > P.y));
}

template<> 
static SE_INLINE bool IsContain(const SeAabb<Float3>& aabb, const Float3& P) 
{ 
	return !((aabb.Upper.x < P.x) || (aabb.Upper.y < P.y) || (aabb.Upper.z < P.z) ||
			 (aabb.Lower.x > P.x) || (aabb.Lower.y > P.y) || (aabb.Lower.z > P.z));
}

/*************************************************************************
***********************    IsInside<AABB-Line>    ************************
*************************************************************************/

template<typename VecT>
static  bool IsIntersect(const SeAabb<VecT>& aabb, const VecT& Pa, const VecT& Pb, const VecT& iv);

template<typename VecT>
static  bool IsIntersect(const SeAabb<VecT>& aabb, const VecT& Pa, const VecT& Pb);


template<>
static SE_INLINE bool IsIntersect(const SeAabb<Float3>& aabb, const Float3& Pa, const Float3& Pb) 
{
	Float3 Dir = Float3(Pb - Pa);

	if (Dir.x == 0.0f)		Dir.x = 1e-6f;
	if (Dir.y == 0.0f)		Dir.y = 1e-6f;
	if (Dir.z == 0.0f)		Dir.z = 1e-6f;

	Float3 invDir = 1.0f / Dir;
	Float3 Left = Float3(aabb.Lower - Pa) * invDir;
	Float3 Right = Float3(aabb.Upper - Pa) * invDir;

	Float2 T1 = Float2(Math::Min(Left.x, Right.x), Math::Max(Left.x, Right.x));
	Float2 T2 = Float2(Math::Min(Left.y, Right.y), Math::Max(Left.y, Right.y));
	Float2 T3 = Float2(Math::Min(Left.z, Right.z), Math::Max(Left.z, Right.z));

	float range0 = Math::Max(T1.x, T2.x, T3.x, 0.0f);
	float range1 = Math::Min(T1.y, T2.y, T3.y, 1.0f);

	return range0 < range1;
}


template<typename VecT>
static SE_INLINE bool IsOverlap(const SeAabb<VecT>& aabb, const SeAabb<VecT>& P) { return false; }


SE_NAMESPACE_END
