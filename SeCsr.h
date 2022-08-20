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

#include "SeOmp.h"
#include "SeArray2D.h"
#include "SeUtility.h"

SE_NAMESPACE_BEGIN

template<typename Type> class SeCompressSparseData
{
public:

	SeCompressSparseData() {}
	SeCompressSparseData(const std::vector<int>& starts, const std::vector<int>& idxs, const std::vector<Type>& values)
		:
		m_starts(starts),
		m_idxs(idxs),
		m_values(values)
	{}

	virtual ~SeCompressSparseData() {}

	void InitIdxs(const std::vector<std::vector<int>> & array2D)
	{
		Clear();

		int size = array2D.size();

		m_starts.resize(size + 1); Utility::MemsetZero(m_starts);

		for (int i = 0; i < size; ++i)
		{
			int num = array2D[i].size();

			m_starts[i + 1] = m_starts[i] + num;
		}

		m_idxs.resize(m_starts.back());

		OMP_PARALLEL_FOR

		for (int i = 0; i < size; ++i)
		{
			int num = array2D[i].size();
			std::memcpy(m_idxs.data() + m_starts[i], array2D[i].data(), sizeof(int) * num);
		}
	}

	void InitIdxs(const SeArray2D<int> & array2D, bool isRowMajor)
	{
		Clear();

		int DIM = isRowMajor ? array2D.Rows() : array2D.Columns();

		m_starts.resize(DIM + 1); Utility::MemsetZero(m_starts);

		for (int i = 0; i < DIM; ++i)
		{
			int num = isRowMajor ? array2D[i][0] : array2D[0][i];

			m_starts[i + 1] = m_starts[i] + num;
		}

		m_idxs.resize(m_starts.back());

		OMP_PARALLEL_FOR

		for (int i = 0; i < DIM; ++i)
		{
			int * address = m_idxs.data() + m_starts[i];

			int num = Size(i);

			for (int k = 1; k <= num; ++k)
			{
				address[k - 1] = isRowMajor ? array2D[i][k] : array2D[k][i];
			}
		}
	}

	void Clear()
	{
		Utility::ClearAndShrink(m_starts);
		Utility::ClearAndShrink(m_idxs);
		Utility::ClearAndShrink(m_values);
	}

	int Size() const 
	{
		return m_starts.back();
	}

	int Size(int id) const
	{
		return -m_starts[id] + m_starts[id + 1];
	}

	int Start(int id) const 
	{
		return m_starts[id];
	}

	const int* StartPtr(int id) const 
	{
		return &(m_starts[id]);
	}

	const int Idx(int id) const
	{
		return m_idxs[id];
	}

	const int * IdxPtr(int id) const
	{
		return m_idxs.data() + m_starts[id];
	}

	const Type * ValuePtr(int id) const 
	{
		return m_values.data() + m_starts[id];
	}

	virtual const SeCompressSparseData * Ptr() const
	{
		return this;
	}

protected:

	std::vector<int>	m_starts;
	std::vector<int>	m_idxs;
	std::vector<Type>	m_values;
};

template<typename Type> class SeCsr : public SeCompressSparseData<Type>
{
public:
	SeCsr() {}
	SeCsr(const std::vector<int>& starts, const std::vector<int>& idxs, const std::vector<Type>& values) :SeCompressSparseData<Type>(starts, idxs, values) {}

	int Rows() const { return SeCompressSparseData<Type>::m_starts.size() - 1; }

	virtual const SeCsr * Ptr() const override
	{
		return this;
	}
};


template<typename Type> class SeCsc : public SeCompressSparseData<Type>
{
public:

	int Columns() const { return SeCompressSparseData<Type>::m_starts.size() - 1; }

	virtual const SeCsc * Ptr() const override
	{
		return this;
	}
};

SE_NAMESPACE_END