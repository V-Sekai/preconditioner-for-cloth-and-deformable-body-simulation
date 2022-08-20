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

SE_NAMESPACE_BEGIN

struct EfSet
{
	int		m_eId;
	int		m_fId;
	float	stiff;
	Float3	m_bary;			// (x, 1-x) / (y, z, 1-y-z): barycentric weight of the intersection point
	SeVec3fSimd m_normal;	// repulsion direction
};

struct VfSet
{
	int		m_vId;
	int		m_fId;
	float	stiff;
	Float2	m_bary;			// (x, y, 1-x-y): barycentric weight of the vertex
	SeVec3fSimd m_normal;	// repulsion direction
};

struct EeSet
{
	int		m_eId0;
	int		m_eId1;
	float	stiff;
	Float2	m_bary;			// (x, 1-x) / (y, 1-y): barycentric weight of the two closest points
	SeVec3fSimd m_normal;	// repulsion direction
};

struct Stencil
{
	int verextNumPerStencil;
	int vertexNumOfFirstPrimitive;

	int		index[5];
	float	weight[5];
	float	stiff;
	SeVec3fSimd direction;
};

SE_NAMESPACE_END