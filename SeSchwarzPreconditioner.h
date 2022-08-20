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

#include "SeCsr.h"
#include "SeMorton.h"
#include "SeAabbSimd.h"
#include "SeCollisionElements.h"


SE_NAMESPACE_BEGIN

class SeSchwarzPreconditioner
{

public:
	
	//==== input data

	const SeVec3fSimd* m_positions = nullptr;		// the mesh nodal positions

	//==== input mesh topology data for collision computation

	const Int4* m_edges = nullptr;				// indices of the two adjacent vertices and the two opposite vertices of edges
	const Int4* m_faces = nullptr;				// indices of the three adjacent vertices of faces, the fourth is useless
	
	const SeCsr<int>* m_neighbours = nullptr;	// the adjacent information of vertices stored in a csr format

public:

	//==== call before time integration once a frame
	void AllocatePrecoditioner(int numVerts, int numEdges, int numFaces);

	//==== call before PCG iteration loop
	void PreparePreconditioner(const SeMatrix3f* diagonal, const SeMatrix3f* csrOffDiagonals, const int* csrRanges,
		const EfSet* efSets, const EeSet* eeSets, const VfSet* vfSets, unsigned int* efCounts, unsigned int* eeCounts, unsigned int* vfCounts);

	//==== call during PCG iterations
	void Preconditioning(SeVec3fSimd* z, const SeVec3fSimd* residual, int dim);

private:

	int m_numVerts = 0;
	int m_numEdges = 0;
	int m_numFaces = 0;

	int m_frameIndex = 0;

	int m_totalSz = 0;
	int m_numLevel = 0;

	int m_totalNumberClusters;

	
	SeAabb<SeVec3fSimd>					m_aabb;


	SeArray2D<SeMatrix3f>				m_hessianMapped;
	std::vector<int>					m_mappedNeighborsNum;
	std::vector<int>					m_mappedNeighborsNumRemain;
	SeArray2D<int>						m_mappedNeighbors;
	SeArray2D<int>						m_mappedNeighborsRemain;

	SeArray2D<int>						m_CoarseSpaceTables;
	std::vector<int>                    m_prefixOrignal;
	std::vector<unsigned int>           m_fineConnectMask;
	std::vector<unsigned int>           m_nextConnectMsk;
	std::vector<unsigned int>           m_nextPrefix;

	std::vector<SeVec3fSimd>			m_mappedPos;

	std::vector<Int4>					m_coarseTables;
	std::vector<int>					m_goingNext;
	std::vector<int>					m_denseLevel;
	std::vector<Int2>					m_levelSize;  // .x = current level size    .y = (prefixed) current level begin index

	std::vector<SeVec3fSimd>			m_mappedR;
	std::vector<SeVec3fSimd>			m_mappedZ;
	std::vector<int>                    m_MapperSortedGetOriginal;          // sorted by morton
	std::vector<int>                    m_mapperOriginalGetSorted;
	std::vector<SeMorton64>             m_mortonCode;

	SeArray2D<SeMatrix3f>				m_hessian32;
	std::vector<SeMatrix3f>             m_additionalHessian32;
	std::vector<float>                  m_invSymR;

	
	int m_stencilNum;
	int m_maxStencilNum;
	std::vector<Stencil>				m_stencils;
	std::vector<Int5>					m_stencilIndexMapped;

private:

	void ComputeLevelNums(int bankSize);

	void DoAlllocation();

	void ComputeTotalAABB();

	void ComputeAABB();

	void SpaceSort();

	void FillSortingData();

	void DoingSort();

	void ComputeInverseMapper();

	void MapHessianTable();

	void PrepareCollisionStencils(const EfSet* efSets, const EeSet* eeSets, const VfSet* vfSets, unsigned int* efCounts, unsigned int* eeCounts, unsigned int* vfCounts);

	void MapCollisionStencilIndices();

	void ReorderRealtime();



	void BuildConnectMaskL0();

	void BuildCollisionConnection(unsigned int* pConnect, const int* pCoarseTable);

	void PreparePrefixSumL0();

	void BuildLevel1();

	void BuildConnectMaskLx(int level);

	void NextLevelCluster(int level);

	void PrefixSumLx(int level);

	void ComputeNextLevel(int level);

	void TotalNodes();

	void AggregationKernel();

	void AdditionalSchwarzHessian2(SeMatrix3f hessian, std::vector<SeMatrix3f>& pAdditionalHessian, SeArray2D<SeMatrix3f>& pDenseHessian, int v1, int v2, const std::vector<int>& pGoingNext, int nLevel, int vNum);

	void PrepareCollisionHessian();

	void PrepareHessian(const SeMatrix3f* diagonal, const SeMatrix3f* csrOffDiagonals, const int* csrRanges);
	
	void LDLtInverse512();

	void BuildResidualHierarchy(const SeVec3fSimd* m_cgResidual);

	void SchwarzLocalXSym();

	void CollectFinalZ(SeVec3fSimd* m_cgZ);
};

SE_NAMESPACE_END