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


SE_NAMESPACE_BEGIN


template<typename Type, int N, int M> class SeMatrix
{

public:

	static constexpr int nRows = N;
	static constexpr int nCols = M;
	static constexpr int nElements = N * M;
	
private:

	union
	{
		Type				m_data[nElements];

		SeVector<Type, N>	m_columns[M];
	};

	 const SeMatrix & Self() const { return *this; }

public:
	
	 SeMatrix() {}

	 SeMatrix(const Type* d) 
	{
		for (int i = 0; i < nElements; ++i)
		{
			m_data[i] = d[i];
		}
	}

	 explicit SeMatrix(const Type & rhs)
	{
		for (int i = 0; i < nElements; ++i)
		{
			m_data[i] = rhs;
		}
	}

	 const	Type & operator() (int i, int j)	const	{ return m_data[j*nRows + i]; }
			Type & operator() (int i, int j)			{ return m_data[j*nRows + i]; }


	 SeMatrix operator - () const
	{
		SeMatrix rhs(Type(0));
		for (int i = 0; i < nElements; ++i)
		{
			rhs.m_data[i] = -m_data[i];
		}
		return rhs;
	}

	 SeMatrix operator + (const SeMatrix & rhs) const	
	{																	
		SeMatrix out(Type(0));											
		for (int i = 0; i < nElements; ++i)
		{
			out.m_data[i] = m_data[i] + rhs.m_data[i];
		}
		return out;
	}
	
	 SeMatrix operator - (const SeMatrix & rhs) const
	{
		SeMatrix out(Type(0));
		for (int i = 0; i < nElements; ++i)
		{
			out.m_data[i] = m_data[i] - rhs.m_data[i];
		}
		return out;
	}

	 SeMatrix & operator += (const SeMatrix & rhs)
	{															
		for (int i = 0; i < nElements; ++i)
		{
			m_data[i] += rhs.m_data[i];
		}
		return *this;							
	}

	 SeMatrix & operator -= (const SeMatrix & rhs)
	{
		for (int i = 0; i < nElements; ++i)
		{
			m_data[i] -= rhs.m_data[i];
		}
		return *this;
	}

	
	 SeMatrix operator * (const Type & rhs) const
	{															
		SeMatrix<Type, N, M> out(Type(0));
		for (int i = 0; i < nElements; ++i)
		{
			out.m_data[i] = m_data[i] * rhs;
		}
		return out;
	}

	 SeMatrix & operator *= (const Type & rhs)
	{													
		for (int i = 0; i < nElements; ++i)
		{
			m_data[i] *= rhs;
		}
		return *this;		
	}

	 SeVector<Type, N> operator * (const SeVector<Type, M> & rhs) const
	{
		SeVector<Type, N> ans(Type(0));
		for (int i = 0; i < N; ++i)
		{
			Type s(0);
			for (int j = 0; j < M; ++j)
			{
				s += Self()(i, j) * rhs[j];
			}
			ans[i] = s;
		}
		return ans;
	}

	template<int K>  SeMatrix<Type, N, K> operator * (const SeMatrix<Type, M, K>& rhs) const
	{
		SeMatrix<Type, N, K> ans(Type(0));
		for (int i = 0; i < N; ++i)
		{
			for (int k = 0; k < K; ++k)
			{
				Type s(0);
				for (int j = 0; j < M; ++j)
				{
					s += Self()(i, j) * rhs(j, k);
				}
				ans(i, k) = s;
			}
		}
		return ans;
	}


	 SeMatrix P2PProduct(const SeMatrix & rhs) const
	{
		SeMatrix out(Type(0));
		for (size_t i = 0; i < nElements; ++i)
		{
			out.m_data[i] = m_data[i] * rhs.m_data[i];
		}
		return out;
	}

	 SeMatrix P2PDivide(const SeMatrix & rhs) const
	{
		SeMatrix out(Type(0));
		for (int i = 0; i < nElements; ++i)
		{
			out.m_data[i] = m_data[i] / rhs.m_data[i];
		}
		return out;
	}

	template<int K, int L>  SeMatrix<Type, N*K, M*L> KroneckerProduct(const SeMatrix<Type, K, L>& rhs) const
	{
		SeMatrix<Type, N*K, M*L> ans(Type(0));
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < M; ++j)
			{
				int rStart = i * K;
				int cStart = j * L;

				Type ij = Self()(i, j);

				for (int p = 0; p < K; ++p)
				{
					for (int q = 0; q < L; ++q)
					{
						ans(rStart + p, cStart + q) = ij * rhs(p, q);
					}
				}
			}
		}
	}


	 SeVector<Type, M> Row(int rId) const
	{
		SeVector<Type, M> x;
		for (int k = 0; k < M; ++k)
		{
			x[k] = Self()(rId, k);
		}
		return x;
	}

	 SeVector<Type, N> Column(int cId) const
	{
		return m_columns[cId];
	}

	template<int K, int L>  SeMatrix<Type, K, L> Block(int row, int col) const
	{
		SeMatrix<Type, K, L> ans(Type(0));
		for (int i = 0; i < K; ++i)
		{
			for (int j = 0; j < L; ++j)
			{
				ans(i, j) = Self()(i + row, j + col);
			}
		}
		return ans;
	}

	 SeVector<Type, SE_MIN(N, M)> Diagonal()const
	{
		SeVector<Type, SE_MIN(N, M)> x;
		for (int i = 0; i < SE_MIN(N, M); ++i)
		{
			x[i] = Self()(i, i);
		}
		return x;
	}

	 void SetRow(const SeVector<Type, M> & row, int rId)
	{
		for (int k = 0; k < M; ++k)
		{
			m_data[rId + k * nRows] = row[k];
		}
	}

	 void SetColumn(const SeVector<Type, N> & column, int cId)
	{
		m_columns[cId] = column;
	}

	template<int K, int L>  void SetBlock(const SeMatrix<Type, K, L> & block, int rowStart, int colStart)
	{
		for (int i = 0; i < K; ++i)
		{
			for (int j = 0; j < L; ++j)
			{
				Self()(i + rowStart, j + colStart) = block(i, j);
			}
		}
	}

	 void SetDiagonal(const SeVector<Type, SE_MIN(N, M)> & diagonal)
	{
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < M; ++j)
			{
				Self()(i, j) = (i == j) ? diagonal[i] : Type(0);
			}
		}
	}

	 SeMatrix<Type, M, N> Transpose() const
	{
		SeMatrix<Type, M, N> trans(Type(0));
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < M; ++j)
			{
				trans(j, i) = Self()(i, j);
			}
		}
		return trans;
	}

	static  SeMatrix Identity(Type diagonal)
	{
		SeMatrix identity = SeMatrix(Type(0));
		for (int i = 0; i < SE_MIN(N, M); ++i)
		{
			identity(i, i) = diagonal;
		}
		return identity;
	}



	 Type FrobeniusNorm() const
	{
		Type sum(0);
		for (int i = 0; i < nElements; ++i)
		{
			sum += m_data[i] * m_data[i];
		}
		return std::sqrt(sum);
	}

	 Type Trace() const 
	{
		Type trace(0);
		for (int i = 0; i < SE_MIN(N, M); ++i)
		{
			trace += Self()(i, i);
		}
		return trace;
	}

};

template<typename Type, int N, int M> 
SE_INLINE SeMatrix<Type, N, M> operator * (const Type & scalar, const SeMatrix<Type, N, M> & mat)
{
	return mat * scalar;
}

template<typename Type, int M, int N>
SE_INLINE SeMatrix<Type, M, N> OuterProduct(const SeVector<Type, M>& vecM, const SeVector<Type, N>& vecN)
{
	SeMatrix<Type, M, N> ans(Type(0));
	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			ans(i, j) = vecM[i] * vecN[j];
		}
	}
	return ans;
}


//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////


//template<typename Type> class _declspec(align(16)) SeMatrix<Type, 2, 2>
template<typename Type> class SE_ALIGN(16) SeMatrix<Type,2,2>
{

public:

	static constexpr int DIM = 2;

	static constexpr int nElements = DIM * DIM;

private:

	union
	{
		Type			m_data[nElements]; 
		
		SeVector2<Type> m_columns[DIM]; 
	};

	 const SeMatrix & Self() const { return *this; }

public:
	
	 SeMatrix() {}

	 explicit SeMatrix(const Type & v)
	{
		m_data[0] = v; m_data[1] = v; m_data[2] = v; m_data[3] = v;
	}

	 SeMatrix(const Type & v0, const Type & v1, const Type & v2, const Type & v3)
	{
		m_data[0] = v0;	m_data[1] = v1; m_data[2] = v2;	m_data[3] = v3;
	}


	 const	Type & operator() (int i, int j)	const	{ return m_data[j*DIM + i]; }
			Type & operator() (int i, int j)			{ return m_data[j*DIM + i]; }

	 SeMatrix operator - () const
	{
		SeMatrix rhs(Type(0));
		for (int i = 0; i < nElements; ++i)
		{
			rhs.m_data[i] = -m_data[i];
		}
		return rhs;
	}

	 SeMatrix operator + (const SeMatrix & rhs) const
	{
		SeMatrix out(Type(0));
		for (int i = 0; i < nElements; ++i)
		{
			out.m_data[i] = m_data[i] + rhs.m_data[i];
		}
		return out;
	}

	 SeMatrix operator - (const SeMatrix & rhs) const
	{
		SeMatrix out(Type(0));
		for (int i = 0; i < nElements; ++i)
		{
			out.m_data[i] = m_data[i] - rhs.m_data[i];
		}
		return out;
	}

	 SeMatrix & operator += (const SeMatrix & rhs)
	{
		for (int i = 0; i < nElements; ++i)
		{
			m_data[i] += rhs.m_data[i];
		}
		return *this;
	}

	 SeMatrix & operator -= (const SeMatrix & rhs)
	{
		for (int i = 0; i < nElements; ++i)
		{
			m_data[i] -= rhs.m_data[i];
		}
		return *this;
	}


	 SeMatrix operator * (const Type & rhs) const
	{
		SeMatrix out(Type(0));
		for (int i = 0; i < nElements; ++i)
		{
			out.m_data[i] = m_data[i] * rhs;
		}
		return out;
	}

	 SeMatrix & operator *= (const Type & rhs)
	{
		for (int i = 0; i < nElements; ++i)
		{
			m_data[i] *= rhs;
		}
		return *this;
	}

	 SeVector<Type, DIM> operator * (const SeVector<Type, DIM> & rhs) const
	{
		SeVector<Type, DIM> ans(Type(0));
		for (int i = 0; i < DIM; ++i)
		{
			Type s(0);
			for (int j = 0; j < DIM; ++j)
			{
				s += Self()(i, j) * rhs[j];
			}
			ans[i] = s;
		}
		return ans;
	}

	template<int K>  SeMatrix<Type, DIM, K> operator * (const SeMatrix<Type, DIM, K>& rhs) const
	{
		SeMatrix<Type, DIM, K> ans(Type(0));
		for (int i = 0; i < DIM; ++i)
		{
			for (int k = 0; k < K; ++k)
			{
				Type s(0);
				for (int j = 0; j < DIM; ++j)
				{
					s += Self()(i, j) * rhs(j, k);
				}
				ans(i, k) = s;
			}
		}
		return ans;
	}


	 SeMatrix P2PProduct(const SeMatrix & rhs) const
	{
		SeMatrix out(Type(0));
		for (size_t i = 0; i < nElements; ++i)
		{
			out.m_data[i] = m_data[i] * rhs.m_data[i];
		}
		return out;
	}

	 SeMatrix P2PDivide(const SeMatrix & rhs) const
	{
		SeMatrix out(Type(0));
		for (int i = 0; i < nElements; ++i)
		{
			out.m_data[i] = m_data[i] / rhs.m_data[i];
		}
		return out;
	}

	template<int K, int L>  SeMatrix<Type, DIM*K, DIM*L> KroneckerProduct(const SeMatrix<Type, K, L>& rhs) const
	{
		SeMatrix<Type, DIM*K, DIM*L> ans(Type(0));
		for (int i = 0; i < DIM; ++i)
		{
			for (int j = 0; j < DIM; ++j)
			{
				int rStart = i * K;
				int cStart = j * L;

				Type ij = Self()(i, j);

				for (int p = 0; p < K; ++p)
				{
					for (int q = 0; q < L; ++q)
					{
						ans(rStart + p, cStart + q) = ij * rhs(p, q);
					}
				}
			}
		}
	}


	 SeVector<Type, DIM> Row(int rId) const
	{
		SeVector<Type, DIM> x;
		x[0] = m_data[rId];
		x[1] = m_data[rId + DIM];
		return x;
	}

	 SeVector<Type, DIM> Column(int cId) const
	{
		return m_columns[cId];
	}

	 SeVector<Type, DIM> Diagonal()const
	{
		SeVector<Type, DIM> x;
		x[0] = m_data[0];
		x[1] = m_data[3];
		return x;
	}

	 void SetRow(const SeVector<Type, DIM> & row, int rId)
	{
		m_data[rId]		= row[0];
		m_data[rId + 2] = row[1];
	}

	 void SetColumn(const SeVector<Type, DIM> & column, int cId)
	{
		m_columns[cId] = column;
	}

	 void SetDiagonal(const SeVector<Type, DIM> & diagonal)
	{
		m_data[0] = diagonal[0]; m_data[1] = 0;
		m_data[2] = 0; m_data[3] = diagonal[1];
	}

	 SeMatrix Transpose() const
	{
		SeMatrix trans = Self();
		trans(0, 1) = Self()(1, 0);
		trans(1, 0) = Self()(0, 1);
		return trans;
	}

	static  SeMatrix Identity(Type diagonal = Type(1))
	{
		SeMatrix identity(Type(0));
		identity(0, 0) = diagonal;
		identity(1, 1) = diagonal;
		return identity;
	}

	 Type FrobeniusNorm() const
	{
		Type sum(0);
		for (int i = 0; i < nElements; ++i)
		{
			sum += m_data[i] * m_data[i];
		}
		return std::sqrt(sum);
	}

	 Type Trace() const
	{	
		return m_data[0] + m_data[3];
	}

	 Type Det() const
	{
		return 	m_data[0] * m_data[3] - m_data[1] * m_data[2];
	}

	 void Inverse(SeMatrix & result) const
	{
		Type dt = 1 / Det();

		result(0, 0) = dt * Self()(1, 1);
		result(1, 0) = -dt * Self()(1, 0);
		result(0, 1) = -dt * Self()(0, 1);
		result(1, 1) = dt * Self()(0, 0);
	}

	 SeMatrix Inverse() const
	{
		SeMatrix r(Self());
		Inverse(r);
		return r;
	}
};



template<typename Type> class SeMatrix<Type, 3, 3>
{

public:

	static constexpr int DIM = 3;

	static constexpr int nElements = DIM * DIM;

private:

	union
	{
		Type			m_data[nElements];
		SeVector3<Type> m_columns[DIM];
	};

	 const SeMatrix & Self() const { return *this; }

public:

	 SeMatrix() {}

	 explicit SeMatrix(const Type & v)
	{
		for (int i = 0; i < nElements; ++i)
		{
			m_data[i] = v;
		}
	}

	 const	Type & operator() (int i, int j) const	{ return m_data[j*DIM + i]; }
			Type & operator() (int i, int j)		{ return m_data[j*DIM + i]; }

	 SeMatrix operator - () const
	{
		SeMatrix rhs(Type(0));
		for (int i = 0; i < nElements; ++i)
		{
			rhs.m_data[i] = -m_data[i];
		}
		return rhs;
	}

	 SeMatrix operator + (const SeMatrix & rhs) const
	{
		SeMatrix out(Type(0));
		for (int i = 0; i < nElements; ++i)
		{
			out.m_data[i] = m_data[i] + rhs.m_data[i];
		}
		return out;
	}

	 SeMatrix operator - (const SeMatrix & rhs) const
	{
		SeMatrix out(Type(0));
		for (int i = 0; i < nElements; ++i)
		{
			out.m_data[i] = m_data[i] - rhs.m_data[i];
		}
		return out;
	}

	 SeMatrix & operator += (const SeMatrix & rhs)
	{
		for (int i = 0; i < nElements; ++i)
		{
			m_data[i] += rhs.m_data[i];
		}
		return *this;
	}

	 SeMatrix & operator -= (const SeMatrix & rhs)
	{
		for (int i = 0; i < nElements; ++i)
		{
			m_data[i] -= rhs.m_data[i];
		}
		return *this;
	}

	 void operator = (const SeMatrix & rhs)
	{
		for (int i = 0; i < nElements; ++i)
		{
			m_data[i] = rhs.m_data[i];
		}
	}


	 SeMatrix operator * (const Type & rhs) const
	{
		SeMatrix out(Type(0));

		for (int i = 0; i < nElements; ++i)
		{
			out.m_data[i] = m_data[i] * rhs;
		}

		return out;
	}

	 SeMatrix operator / (const Type & rhs) const
	{
		SeMatrix out(Type(0));

		for (int i = 0; i < nElements; ++i)
		{
			out.m_data[i] = m_data[i] / rhs;
		}

		return out;
	}

	 SeMatrix & operator *= (const Type & rhs)
	{
		for (int i = 0; i < nElements; ++i)
		{
			m_data[i] *= rhs;
		}

		return *this;
	}

	 SeMatrix & operator /= (const Type & rhs)
	{
		for (int i = 0; i < nElements; ++i)
		{
			m_data[i] /= rhs;
		}

		return *this;
	}

	 SeVector<Type, DIM> operator * (const SeVector<Type, DIM> & rhs) const
	{
		SeVector<Type, DIM> ans(Type(0));
		for (int i = 0; i < DIM; ++i)
		{
			Type s(0);
			for (int j = 0; j < DIM; ++j)
			{
				s += Self()(i, j) * rhs[j];
			}
			ans[i] = s;
		}
		return ans;
	}

	template<int K>  SeMatrix<Type, DIM, K> operator * (const SeMatrix<Type, DIM, K>& rhs) const
	{
		SeMatrix<Type, DIM, K> ans(Type(0));
		for (int i = 0; i < DIM; ++i)
		{
			for (int k = 0; k < K; ++k)
			{
				Type s(0);
				for (int j = 0; j < DIM; ++j)
				{
					s += Self()(i, j) * rhs(j, k);
				}
				ans(i, k) = s;
			}
		}
		return ans;
	}


	 SeMatrix P2PProduct(const SeMatrix & rhs) const
	{
		SeMatrix out(Type(0));
		for (size_t i = 0; i < nElements; ++i)
		{
			out.m_data[i] = m_data[i] * rhs.m_data[i];
		}
		return out;
	}

	 SeMatrix P2PDivide(const SeMatrix & rhs) const
	{
		SeMatrix out(Type(0));
		for (int i = 0; i < nElements; ++i)
		{
			out.m_data[i] = m_data[i] / rhs.m_data[i];
		}
		return out;
	}

	template<int K, int L>  SeMatrix<Type, DIM*K, DIM*L> KroneckerProduct(const SeMatrix<Type, K, L>& rhs) const
	{
		SeMatrix<Type, DIM*K, DIM*L> ans(Type(0));
		for (int i = 0; i < DIM; ++i)
		{
			for (int j = 0; j < DIM; ++j)
			{
				int rStart = i * K;
				int cStart = j * L;

				Type ij = Self()(i, j);

				for (int p = 0; p < K; ++p)
				{
					for (int q = 0; q < L; ++q)
					{
						ans(rStart + p, cStart + q) = ij * rhs(p, q);
					}
				}
			}
		}
	}


	 SeVector<Type, DIM> Row(int rId) const
	{
		SeVector<Type, DIM> x;
		x[0] = m_data[rId];
		x[1] = m_data[rId + 3];
		x[2] = m_data[rId + 6];
		return x;
	}

	 SeVector<Type, DIM> Column(int cId) const
	{
		return m_columns[cId];
	}

	 SeVector<Type, DIM> Diagonal()const
	{
		SeVector<Type, DIM> x;
		x[0] = m_data[0];
		x[1] = m_data[4];
		x[2] = m_data[8];
		return x;
	}

	 void SetRow(const SeVector<Type, DIM> & row, int rId)
	{
		m_data[rId] = row[0];
		m_data[rId + 3] = row[1];
		m_data[rId + 6] = row[2];
	}

	 void SetColumn(const SeVector<Type, DIM> & column, int cId)
	{
		m_columns[cId] = column;
	}

	 void SetDiagonal(const SeVector<Type, DIM> & diagonal)
	{
		Type zero(0);

		m_data[0] = diagonal[0];	m_data[1] = zero;			m_data[2] = zero;
		m_data[3] = zero;			m_data[4] = diagonal[1];	m_data[5] = zero;
		m_data[6] = zero;			m_data[7] = zero;			m_data[8] = diagonal[2];
	}

	 SeMatrix Transpose() const
	{
		SeMatrix trans(Type(0));
		for (int i = 0; i < DIM; ++i)
		{
			for (int j = 0; j < DIM; ++j)
			{
				trans(i, j) = Self()(j, i);
			}
		}
		return trans;
	}

	static  SeMatrix Identity(Type diagonal = Type(1))
	{
		SeMatrix identity(Type(0));
		identity(0, 0) = diagonal;
		identity(1, 1) = diagonal;
		identity(2, 2) = diagonal;
		return identity;
	}

	 Type FrobeniusNorm() const
	{
		Type sum(0);
		for (int i = 0; i < nElements; ++i)
		{
			sum += m_data[i] * m_data[i];
		}
		return std::sqrt(sum);
	}

	 Type Trace() const
	{
		return m_data[0] + m_data[4] + m_data[8];
	}

	 Type Det() const
	{
		return 	m_data[0] * (m_data[4] * m_data[8] - m_data[7] * m_data[5])
			- m_data[3] * (m_data[1] * m_data[8] - m_data[7] * m_data[2])
			+ m_data[6] * (m_data[1] * m_data[5] - m_data[4] * m_data[2]);
	}

	 void Inverse(SeMatrix & result) const
	{
		Type dt = 1 / Det();

		result(0, 0) = dt * (Self()(1, 1) * Self()(2, 2) - Self()(1, 2) * Self()(2, 1));
		result(1, 0) = dt * (Self()(1, 2) * Self()(2, 0) - Self()(1, 0) * Self()(2, 2));
		result(2, 0) = dt * (Self()(1, 0) * Self()(2, 1) - Self()(1, 1) * Self()(2, 0));

		result(0, 1) = dt * (Self()(0, 2) * Self()(2, 1) - Self()(0, 1) * Self()(2, 2));
		result(1, 1) = dt * (Self()(0, 0) * Self()(2, 2) - Self()(0, 2) * Self()(2, 0));
		result(2, 1) = dt * (Self()(0, 1) * Self()(2, 0) - Self()(0, 0) * Self()(2, 1));

		result(0, 2) = dt * (Self()(0, 1) * Self()(1, 2) - Self()(0, 2) * Self()(1, 1));
		result(1, 2) = dt * (Self()(0, 2) * Self()(1, 0) - Self()(0, 0) * Self()(1, 2));
		result(2, 2) = dt * (Self()(0, 0) * Self()(1, 1) - Self()(0, 1) * Self()(1, 0));
	}
	
	 SeMatrix Inverse() const
	{
		SeMatrix r(Self());
		Inverse(r);
		return r;
	}
};

template<typename Type>
SE_INLINE SeMatrix<Type, 3, 3> operator * (const Type & scalar, const SeMatrix<Type, 3, 3> & mat)
{
	return mat * scalar;
}

/*inline Simd3f operator * (const SeMatrix<float, 3, 3> & mat, const Simd3f & vec)
{
	return Simd3f(mat * vec.xyz);
}*/


//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////


template<typename Type> class SE_ALIGN(16) SeMatrix<Type, 4, 4>
{

public:

	static constexpr int DIM = 4;
	static constexpr int nElements = DIM * DIM;

private:

	union
	{
		Type			m_data[nElements];
		
		SeVector4<Type> m_columns[DIM];
	};

	 const SeMatrix & Self() const { return *this; }

public:

	 SeMatrix() {}

	 explicit SeMatrix(const Type & scalar)
	{
		for (int i = 0; i < nElements; ++i)
		{
			m_data[i] = scalar;
		}
	}
	

	 const	Type & operator() (int i, int j)const	{ return m_data[j*DIM + i]; }
			Type & operator() (int i, int j)		{ return m_data[j*DIM + i]; }

	 SeMatrix operator - () const
	{
		SeMatrix rhs(Type(0));
		for (int i = 0; i < nElements; ++i)
		{
			rhs.m_data[i] = -m_data[i];
		}
		return rhs;
	}

	 SeMatrix operator + (const SeMatrix & rhs) const
	{
		SeMatrix out(Type(0));
		for (int i = 0; i < nElements; ++i)
		{
			out.m_data[i] = m_data[i] + rhs.m_data[i];
		}
		return out;
	}

	 SeMatrix operator - (const SeMatrix & rhs) const
	{
		SeMatrix out(Type(0));
		for (int i = 0; i < nElements; ++i)
		{
			out.m_data[i] = m_data[i] - rhs.m_data[i];
		}
		return out;
	}

	 SeMatrix & operator += (const SeMatrix & rhs)
	{
		for (int i = 0; i < nElements; ++i)
		{
			m_data[i] += rhs.m_data[i];
		}
		return *this;
	}

	 SeMatrix & operator -= (const SeMatrix & rhs)
	{
		for (int i = 0; i < nElements; ++i)
		{
			m_data[i] -= rhs.m_data[i];
		}
		return *this;
	}


	 SeMatrix operator * (const Type & rhs) const
	{
		SeMatrix out(Type(0));
		for (int i = 0; i < nElements; ++i)
		{
			out.m_data[i] = m_data[i] * rhs;
		}
		return out;
	}

	 SeMatrix & operator *= (const Type & rhs)
	{
		for (int i = 0; i < nElements; ++i)
		{
			m_data[i] *= rhs;
		}
		return *this;
	}

	 SeVector<Type, DIM> operator * (const SeVector<Type, DIM> & rhs) const
	{
		SeVector<Type, DIM> ans(Type(0));
		for (int i = 0; i < DIM; ++i)
		{
			Type s(0);
			for (int j = 0; j < DIM; ++j)
			{
				s += Self()(i, j) * rhs[j];
			}
			ans[i] = s;
		}
		return ans;
	}

	template<int K>  SeMatrix<Type, DIM, K> operator * (const SeMatrix<Type, DIM, K>& rhs) const
	{
		SeMatrix<Type, DIM, K> ans(Type(0));
		for (int i = 0; i < DIM; ++i)
		{
			for (int k = 0; k < K; ++k)
			{
				Type s(0);
				for (int j = 0; j < DIM; ++j)
				{
					s += Self()(i, j) * rhs(j, k);
				}
				ans(i, k) = s;
			}
		}
		return ans;
	}


	 SeMatrix P2PProduct(const SeMatrix & rhs) const
	{
		SeMatrix out(Type(0));
		for (size_t i = 0; i < nElements; ++i)
		{
			out.m_data[i] = m_data[i] * rhs.m_data[i];
		}
		return out;
	}

	 SeMatrix P2PDivide(const SeMatrix & rhs) const
	{
		SeMatrix out(Type(0));
		for (int i = 0; i < nElements; ++i)
		{
			out.m_data[i] = m_data[i] / rhs.m_data[i];
		}
		return out;
	}

	template<int K, int L>  SeMatrix<Type, DIM*K, DIM*L> KroneckerProduct(const SeMatrix<Type, K, L>& rhs) const
	{
		SeMatrix<Type, DIM*K, DIM*L> ans(Type(0));
		for (int i = 0; i < DIM; ++i)
		{
			for (int j = 0; j < DIM; ++j)
			{
				int rStart = i * K;
				int cStart = j * L;

				Type ij = Self()(i, j);

				for (int p = 0; p < K; ++p)
				{
					for (int q = 0; q < L; ++q)
					{
						ans(rStart + p, cStart + q) = ij * rhs(p, q);
					}
				}
			}
		}
	}


	 SeVector<Type, DIM> Row(int rId) const
	{
		SeVector<Type, DIM> x;
		x[0] = m_data[rId];
		x[1] = m_data[rId + DIM];
		x[2] = m_data[rId + DIM * 2];
		x[3] = m_data[rId + DIM * 3];
		return x;
	}

	 SeVector<Type, DIM> Column(int cId) const
	{
		return m_columns[cId];
	}

	 SeVector<Type, DIM> Diagonal()const
	{
		SeVector<Type, DIM> x;
		x[0] = m_data[0];
		x[1] = m_data[5];
		x[2] = m_data[10];
		x[3] = m_data[15];
		return x;
	}

	 void SetRow(const SeVector<Type, DIM> & row, int rId)
	{
		m_data[rId] = row[0];
		m_data[rId + DIM] = row[1];
		m_data[rId + DIM * 2] = row[2];
		m_data[rId + DIM * 3] = row[3];
	}

	 void SetColumn(const SeVector<Type, DIM> & column, int cId)
	{
		m_columns[cId] = column;
	}

	 void SetDiagonal(const SeVector<Type, DIM> & diagonal)
	{
		for (int i = 0; i < DIM; ++i)
		{
			for (int j = 0; j < DIM; ++j)
			{
				Self()(i, j) = (i == j) ? diagonal[i] : Type(0);
			}
		}
	}

	 SeMatrix Transpose() const
	{
		SeMatrix trans(Type(0));
		for (int i = 0; i < DIM; ++i)
		{
			for (int j = 0; j < DIM; ++j)
			{
				trans(j, i) = Self()(i, j);
			}
		}
		return trans;
	}

	static  SeMatrix Identity(Type diagonal = Type(1))
	{
		SeMatrix identity(Type(0));
		identity(0, 0) = diagonal;
		identity(1, 1) = diagonal;
		identity(2, 2) = diagonal;
		identity(3, 3) = diagonal;
		return identity;
	}

	 Type FrobeniusNorm() const
	{
		Type sum(0);
		for (int i = 0; i < nElements; ++i)
		{
			sum += m_data[i] * m_data[i];
		}
		return std::sqrt(sum);
	}

	 Type Trace() const
	{
		return m_data[0] + m_data[5] + m_data[10] + m_data[15];
	}

	

	 SeMatrix<Type, 3, 3> GetRotationPart()	 const
	{
		SeMatrix<Type, 3, 3> R;
		for (int y = 0; y < 3; ++y)
			for (int x = 0; x < 3; ++x)
				R(y, x) = (*this)(y, x);
		return R;
	}
	
	 SeVector3<Type> GetTranslationPart() const
	{
		SeVector3<Type> t;
		for (int y = 0; y < 3; ++y)
			t[y] = (*this)(y, 3);
		return t;
	}
	
	 void SetRotationPart(const SeMatrix<Type, 3, 3> & R)
	{
		for (int y = 0; y < 3; ++y)
			for (int x = 0; x < 3; ++x)
			{
				Self()(y, x) = R(y, x);
			}
	}
	
	 void SetTranslationPart(const SeVector3<Type> & t)
	{
		for (int y = 0; y < 3; ++y)
		{
			Self()(y, 3) = t[y];
		}
	}
};

//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////

template<typename Type, int N, int M, int K, int L> static  SeMatrix<Type, N* K, M* L> Kronecker(const SeMatrix<Type, N, M>& left, const SeMatrix<Type, K, L>& rhs)
{
	SeMatrix<Type, N* K, M* L> ans(Type(0));
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < M; ++j)
		{
			int rStart = i * K;
			int cStart = j * L;

			Type ij = left(i, j);

			for (int p = 0; p < K; ++p)
			{
				for (int q = 0; q < L; ++q)
				{
					ans(rStart + p, cStart + q) = ij * rhs(p, q);
				}
			}
		}
	}
	return ans;
}


template<typename Type, int N, int M> static  SeMatrix<Type, N, 1> GetColumn(const SeMatrix<Type, N, M>& a, int ci)
{
	SeMatrix<Type, N, 1> ret;
	for (int i = 0; i < N; i++)
	{
		ret(i, 0) = a(i, ci);
	}
	return ret;
}




//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////

template<typename Type, int N, int M> SE_INLINE bool HasNan(const SeMatrix<Type, N, M> & mat)
{
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < M; ++j)
		{
			if (std::isnan(mat(i, j)))
				return true;
		}
	}
	return false;
}

template<typename T,int N>
struct AnonymousArray { T data[N]; };

template<typename Type, int Row, int Col>
static  SeMatrix<Type, Row, Col> CreateWithColumnMajor(const AnonymousArray<Type, Row* Col>& a)
{
	SeMatrix<Type, Row, Col> m;

	for (int ri = 0; ri < Row; ri++)
	{
		for (int ci = 0; ci < Col; ci++)
		{
			m(ri, ci) = a.data[ri + ci * Row];
		}
	}

	return m;
}

//////////////////////////////////////////////////////////////////////////

template<typename Type> using SeMatrix2 = SeMatrix<Type, 2, 2>;
template<typename Type> using SeMatrix3 = SeMatrix<Type, 3, 3>;
template<typename Type> using SeMatrix4 = SeMatrix<Type, 4, 4>;
template<typename Type> using SeMatrix9 = SeMatrix<Type, 9, 9>;
template<typename Type> using SeMatrix12 = SeMatrix<Type, 12, 12>;

using SeMatrix2f = SeMatrix2<float>;
using SeMatrix3f = SeMatrix3<float>;
using SeMatrix4f = SeMatrix4<float>;
using SeMatrix9f = SeMatrix9<float>;
using SeMatrix12f = SeMatrix12<float>;

using SeMatrix2i = SeMatrix2<int>;
using SeMatrix3i = SeMatrix3<int>;
using SeMatrix4i = SeMatrix4<int>;

using SeMatrix23f = SeMatrix<float, 2, 3>;
using SeMatrix32f = SeMatrix<float, 3, 2>;
using SeMatrix34f = SeMatrix<float, 3, 4>;
using SeMatrix43f = SeMatrix<float, 4, 3>;


SE_NAMESPACE_END
