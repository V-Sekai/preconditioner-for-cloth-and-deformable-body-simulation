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

#include "SePreDefine.h"

SE_NAMESPACE_BEGIN

namespace Utility
{
	template <typename Type>
	size_t ByteSize(const std::vector<Type> & data)
	{
		return sizeof(Type) * data.size();
	}

	template <typename Type>
	void Memcpy(Type * dstData, const Type * srcData, size_t size)
	{
		std::memcpy(dstData, srcData, size * sizeof(Type));
	}

	template <typename Type1, typename Type2>
	void Memcpy(std::vector<Type1> & data, const std::vector<Type2> & srcData, size_t size)
	{
		/*if(ByteSize(data) < ByteSize(srcData))
		{
			SE_WARNING_LOG("DesSize is less than Src Size!");
		}*/
		std::memcpy(data.data(), srcData.data(), SE_MIN(SE_MIN(ByteSize(data), ByteSize(srcData)), size * sizeof(Type1)));
	}

	template <typename Type>
	void Memset(std::vector<Type>& data, const Type & value)
	{
		std::fill(data.begin(), data.end(), value);
	}

	template <typename Type>
	void Memset(Type* data, const Type& value, size_t size)
	{
		std::fill(data, data + size, value);
	}

	template <typename Type>
	void MemsetZero(std::vector<Type> & data)
	{
		std::memset(data.data(), 0, ByteSize(data));
	}

	template <typename Type>
	void MemsetZero(Type * data, size_t size)
	{
		std::memset(data, 0, sizeof(Type) * size);
	}

	template <typename Type>
	void MemsetMinusOne(std::vector<Type> & data)
	{
		std::memset(data.data(), -1, ByteSize(data));
	}

	template <typename Type>
	void MemsetMinusOne(Type * data, size_t size)
	{
		std::memset(data, -1, sizeof(Type) * size);
	}

	template <typename Type>
	void ClearAndShrink(std::vector<Type> & data)
	{
		data.clear(); data.shrink_to_fit();
	}

	template <typename Type>
	void ResizeAndShrink(std::vector<Type> & data, size_t dim)
	{
		data.resize(dim); data.shrink_to_fit();
	}

	template <typename Type>
	void CopyAndShrink(std::vector<Type> & desData, const std::vector<Type> & srcData)
	{
		desData = srcData; desData.shrink_to_fit();
	}
}

SE_NAMESPACE_END
