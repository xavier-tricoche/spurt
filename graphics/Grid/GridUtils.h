#ifndef GRID_UTILS_H
#define GRID_UTILS_H

#include "Grid3D_Regular.h"
#include "Grid3D_Regular_Sparse_I1.h"
#include "Grid3D_Regular_CSR.h"
#include "Grid3D_Regular_Sparse_NoData.h"
#include "Grid3D_Regular_Block.h"
#include "Grid3D_Regular_Block_NoData.h"

#include <string>
using namespace std;

namespace GridCreation {
	template <class T> inline Grid3D_Regular_Base<T>* create(int type) {
		switch (type) {
			case GRID_REGULAR: return new Grid3D_Regular<T>();
			case GRID_REGULAR_SPARSE_I1: return new Grid3D_Regular_Sparse_I1<T>();
			//case GRID_REGULAR_SPARSE_I2: return new Grid3D_Regular_Sparse_I2<T>();
			//case GRID_REGULAR_SPARSE_I3: return new Grid3D_Regular_Sparse_I3<T>();
			//case GRID_REGULAR_SPARSE_VEC: return new Grid3D_Regular_Sparse_Vec<T>();
			case GRID_REGULAR_CSR: return new Grid3D_Regular_CSR<T>();
			case GRID_REGULAR_SPARSE_NODATA: return new Grid3D_Regular_Sparse_NoData<T>();
			case GRID_POW2: return new Grid3D_Pow2<T>();
			case GRID_POW2_CSR: return new Grid3D_Pow2_CSR<T>();
			case GRID_POW2_SPARSE_NODATA: return new Grid3D_Pow2_Sparse_NoData<T>();
			case GRID_REGULAR_BLOCK: return new Grid3D_Regular_Block<T>();
			case GRID_REGULAR_BLOCK_NODATA: return new Grid3D_Regular_Block_NoData<T>();
			default: printf("Unknown grid type passed to GridCreation::create\n");  return NULL;
		}
	}

	template <class T> inline Grid3D_Regular_Base<T>* createFrom(int outType, const Grid3D_Regular_Base<T> *grid, bool copyData) {
		Grid3D_Regular_Base<T> *newGrid = create<T>(outType);
		newGrid->rebuildGrid(grid->numCells);
		newGrid->setGridPosition(grid->startPos, grid->cellSize);
		newGrid->outsideVal = grid->outsideVal;

		if (!copyData) return newGrid;

		if (outType == GRID_POW2_CSR || outType == GRID_POW2_SPARSE_NODATA) {
			printf("to do: this will not copy quite right using space filling curve grid\n");
		}

		if (newGrid->getType() == GRID_REGULAR || newGrid->getType() == GRID_POW2) {
			if (grid->getType() == GRID_REGULAR || grid->getType() == GRID_POW2) {
				for (unsigned int i=0; i<grid->numCells.x; i++)
				for (unsigned int j=0; j<grid->numCells.y; j++)
				for (unsigned int k=0; k<grid->numCells.z; k++) {
					newGrid->set(i,j,k, grid->get(i,j,k));
				}
			}
			else if (grid->getType() == GRID_REGULAR_SPARSE_I1) {
				Grid3D_Regular_Sparse_I1<T> *unfixedRep = (Grid3D_Regular_Sparse_I1<T>*)grid;

				for (unsigned int x=0; x<unfixedRep->numCells.x; x++)
				for (unsigned int y=0; y<unfixedRep->numCells.y; y++)
				for (unsigned int z=0; z<unfixedRep->depthIndices[x][y].size(); z++) {
					newGrid->set(x, y, unfixedRep->depthIndices[x][y][z], unfixedRep->data[x][y][z]);
				}
			}
			else printf("Unknown copy from grid type given copy to type REGULAR\n");
		}

		else if (newGrid->getType() == GRID_REGULAR_CSR || newGrid->getType() == GRID_POW2_CSR) {
			if (grid->getType() == GRID_REGULAR || grid->getType() == GRID_POW2) {
				((Grid3D_Regular_CSR<T>*)newGrid)->copyDataFrom(*(Grid3D_Regular<T>*)grid, grid->outsideVal);
			}
			else if (grid->getType() == GRID_REGULAR_SPARSE_I1) {
				((Grid3D_Regular_CSR<T>*)newGrid)->copyDataFrom(*(Grid3D_Regular_Sparse_I1<T>*)grid);
			}
			else printf("Unknown copy from grid type given copy to type CSR\n");
		}

		else if (newGrid->getType() == GRID_REGULAR_SPARSE_NODATA || newGrid->getType() == GRID_POW2_SPARSE_NODATA) {
			if (grid->getType() == GRID_REGULAR || grid->getType() == GRID_POW2) {
				((Grid3D_Regular_Sparse_NoData<T>*)newGrid)->copyDataFrom(*(Grid3D_Regular<T>*)grid, grid->outsideVal);
			}
			else if (grid->getType() == GRID_REGULAR_SPARSE_I1) {
				((Grid3D_Regular_Sparse_NoData<T>*)newGrid)->copyDataFrom(*(Grid3D_Regular_Sparse_I1<T>*)grid);
			}
			else printf("Unknown copy from grid type given copy to type CSR (no data)\n");
		}

		else if (newGrid->getType() == GRID_REGULAR_BLOCK) {
			printf("Unknown copy from grid type given copy to type BLOCK\n");
		}
		else if (newGrid->getType() == GRID_REGULAR_BLOCK_NODATA) {
			if (grid->getType() == GRID_REGULAR_BLOCK) {
				((Grid3D_Regular_Block_NoData<T>*)newGrid)->copyDataFrom(*(Grid3D_Regular_Block<T>*)grid);
			}
			else printf("Unknown copy to grid type given copy to type BLOCK (no data)\n");
		}

		else printf("Unknown copy to grid type\n");

		return newGrid;
	}
}

namespace GridSize {
	template <class T> inline double calculateSparsityRatio(const Grid3D_Regular_Base<T> *grid) {
		if (grid->getType() == GRID_REGULAR_SPARSE_I1) {
			return (double)(((Grid3D_Regular_Sparse_I1<T>*)grid)->numberOfEntries)/(double)(grid->numCells.x*grid->numCells.y*grid->numCells.z);
		}
		else if (grid->getType() == GRID_REGULAR || grid->getType() == GRID_POW2) {
			Grid3D_Regular<T> *reg = (Grid3D_Regular<T>*)grid;
			double numberOfNonNull = 0;
			for (unsigned int n=0; n<reg->size(); n++) {
				if (reg->get(n) != reg->outsideVal) numberOfNonNull++;
			}
			return (double)numberOfNonNull/(double)reg->size();
		}
		else if (grid->getType() == GRID_REGULAR_BLOCK) return ((Grid3D_Regular_Block<T>*)grid)->getSparsityRatio();
		else if (grid->getType() == GRID_REGULAR_BLOCK_NODATA) return ((Grid3D_Regular_Block_NoData<T>*)grid)->getSparsityRatio();

		printf("Unknown grid type in GridSize::getSparsityRatio()\n");
		return 0;
	}

	template <class T> inline double getActualSparsityRatio(const Grid3D_Regular_Base<T> *grid, const double sparsityRatio) {
		int type = grid->getType();
		if (type == GRID_REGULAR_SPARSE_NODATA || type == GRID_POW2_SPARSE_NODATA) {
			if (!((Grid3D_Regular_Sparse_NoData<T>*)grid)->hasNormalMapping()) return 1.0 - sparsityRatio;
		}
		return sparsityRatio;
	}

	inline double getTheoreticalSize(int type, double sparsityRatio, Vector3ui numCells, double sizePerElement, double sizePerNullElement) {
		double N = numCells.x*numCells.y*numCells.z; // total possible points
		double M = N*sparsityRatio; // non-NULL points

		if (type == GRID_REGULAR || type == GRID_POW2) {
			return sizePerElement*M + sizePerNullElement*(N-M);
		}
		else if (type == GRID_REGULAR_CSR || type == GRID_POW2_CSR) {
			double idx1Size = numCells.x*numCells.y * sizeof(unsigned int);
			double idx2Size = M * sizeof(unsigned int);
			double dataSize = M * sizePerElement;
			return idx1Size + idx2Size + dataSize;
		}
		else if (type == GRID_REGULAR_SPARSE_NODATA || type == GRID_POW2_SPARSE_NODATA) {
			if (sparsityRatio >= 0.5) M = N * (sparsityRatio-0.5);

			double idx1Size = numCells.x*numCells.y * sizeof(unsigned int);
			double idx2Size = M * sizeof(unsigned int);
			return idx1Size + idx2Size;
		}
		/*else if (type == GRID_REGULAR_BLOCK) {
			
		}*/
		
		printf("Unknown grid type passed to GridSize::getTheoreticalSize()\n");

		return 0;
	}

	template <class T> inline void getSizes(const Grid3D_Regular_Base<T> *grid, double &totalSize, double &sizeOfData) {
		totalSize = sizeOfData = 0;
		
		if (grid->getType() == GRID_REGULAR) {
			Grid3D_Regular<T> *mat = (Grid3D_Regular<T>*)grid;
			double numBytes = mat->getNumberOfBytes();
			totalSize += numBytes;
			sizeOfData += numBytes;
		}
		else if (grid->getType() == GRID_REGULAR_CSR) {
			Grid3D_Regular_CSR<T> *mat = (Grid3D_Regular_CSR<T>*)grid;
			double numBytes = mat->getNumberOfBytes();
			totalSize += numBytes;

			vector<unsigned int> *rowCol = mat->getRowColIndexPtr();
			vector<unsigned int> *depth = mat->getDepthIndexPtr();
			vector<T> *dataPtr = mat->getDataPtr();

			sizeOfData += numBytes - (sizeof(*rowCol) + sizeof(*depth) + sizeof(*dataPtr) +
									  rowCol->size()*sizeof(unsigned int) + depth->size()*sizeof(unsigned int));
		}
		else if (grid->getType() == GRID_REGULAR_SPARSE_NODATA) {
			totalSize += grid->getNumberOfBytes();
		}
		else if (grid->getType() == GRID_REGULAR_BLOCK) {
			totalSize += grid->getNumberOfBytes();
		}
		else printf("Unknown grid type in GridSize::getSizes()\n");

		totalSize /= 1028.0*1028.0;
		sizeOfData /= 1028.0*1028.0;
	}
}

namespace GridConversion {
	inline bool canBeRepresentedAsFullMatrix(const Vector3ui &numCells) {
		Vector3<__int64> cells(numCells.x, numCells.y, numCells.z);
		__int64 maxIdx = (cells.x-1)*cells.y*cells.z + (cells.y-1)*cells.z + (cells.z-1);
		return maxIdx == (numCells.x-1)*numCells.y*numCells.z + (numCells.y-1)*numCells.z + (numCells.z-1);
	}

	template <class T> inline int calculateMostSpaceEfficientRep(Grid3D_Regular_Base<T> *grid, double sparsityRatio, bool useSpaceFillingCurve) {
		Vector3ui numCells = grid->numCells;
		if (useSpaceFillingCurve) {
			numCells.x = MathUtils::roundUpPowerOf2(numCells.x);
			numCells.y = MathUtils::roundUpPowerOf2(numCells.y);
			numCells.z = MathUtils::roundUpPowerOf2(numCells.z);
		}

		double sizePerNullElement = TemplateSize::getAverageNumberOfBytesPerElement(grid->outsideVal);

		double regSize = 0, sizePerElement = 0;
		if (grid->getType() == GRID_REGULAR_SPARSE_I1) {
			Grid3D_Regular_Sparse_I1<T> *unfixedRep = (Grid3D_Regular_Sparse_I1<T>*)grid;

			sizePerElement = unfixedRep->getAverageNumberOfBytesPerElement();
			
			regSize = GridSize::getTheoreticalSize(GRID_REGULAR, sparsityRatio, numCells, sizePerElement, sizePerNullElement);	
		}
		else if (grid->getType() == GRID_REGULAR || grid->getType() == GRID_POW2) {
			Grid3D_Regular<T> *reg = (Grid3D_Regular<T>*)grid;

			regSize = reg->getNumberOfBytes();
			sizePerElement = regSize / (double)reg->size(); // to do: not sure this is the best estimate
		}
		else {
			printf("Unknown grid type in GridConversion::calculateMostSpaceEfficientRep()\n");
			return UNKNOWN_GRID_TYPE;
		}

		double csrSize = GridSize::getTheoreticalSize(GRID_REGULAR_CSR, sparsityRatio, numCells, sizePerElement, sizePerNullElement);

		if (regSize <= csrSize) {
			if (useSpaceFillingCurve) return GRID_POW2;
			else return GRID_REGULAR;
		}
		else {
			if (useSpaceFillingCurve) return GRID_POW2_CSR;
			else return GRID_REGULAR_CSR;
		}
	}
}

namespace GridString {
	inline string get(const int type) {
		switch (type) {
			case GRID_REGULAR: return "Regular";
			case GRID_REGULAR_SPARSE_I1: return "Nested Sparse Vectors(1)";
			case GRID_REGULAR_SPARSE_I2: return "Nested Sparse Vectors(2)";
			case GRID_REGULAR_SPARSE_I3: return "Nested Sparse Vectors(3)";
			case GRID_REGULAR_SPARSE_VEC: return "Sparse Vector";
			case GRID_REGULAR_CSR: return "CSR";
			case GRID_REGULAR_SPARSE_NODATA: return "CSR (no data)";
			case GRID_POW2: return "Regular (pow 2)";
			case GRID_POW2_CSR: return "CSR (pow 2)";
			case GRID_POW2_SPARSE_NODATA: return "CSR (no data, pow 2)";
			case GRID_REGULAR_BLOCK: return "Sparse Block";
			case GRID_REGULAR_BLOCK_NODATA: return "Sparse Block (no data)";
			default: printf("Unknown grid type passed to GridString::get\n");  return NULL;
		}
	}
}

namespace GridType {
	template <class T> inline bool isVector(const Grid3D_Regular_Base<T> *grid) { return false; }
	template <class T> inline bool isVector(const Grid3D_Regular_Base< vector<T> > *grid) { return true; }
	template <class T> inline bool isVector(const Grid3D_Regular_Base< StlVector<T> > *grid) { return true; }	
}

namespace GridCuda {
	template <class T> inline void getMatrixTypeAndIndexSizes(const Grid3D_Regular_Base<T> *grid, int &matrixType,
															  unsigned int &rowColIndexSize, unsigned int &depthIndexSize,
															  unsigned int &blockIndexSize, bool forceSparseIfVectorType) {
		int type = grid->getType();

		matrixType = -1;
		rowColIndexSize = depthIndexSize = blockIndexSize = 0;

		if (type == GRID_REGULAR || type == GRID_POW2) {
			Grid3D_Regular<T> *full = (Grid3D_Regular<T>*)grid;

			if (forceSparseIfVectorType && TemplateType::isVector(full->outsideVal)) {
				Grid3D_Regular_CSR<T> *csr = new Grid3D_Regular_CSR<T>(*full, full->outsideVal);
				getMatrixTypeAndIndexSizes(csr, matrixType, rowColIndexSize, depthIndexSize, blockIndexSize, forceSparseIfVectorType);
				delete csr;  csr = NULL;
			}
			else matrixType = 0;
		}
		else if (type == GRID_REGULAR_CSR || type == GRID_POW2_CSR) {
			Grid3D_Regular_CSR<T> *csr = (Grid3D_Regular_CSR<T>*)grid;

			matrixType = 1;
			rowColIndexSize = csr->getRowColIndexSize();
			depthIndexSize = csr->getDepthIndexSize();
		}
		else if (type == GRID_REGULAR_SPARSE_NODATA || type == GRID_POW2_SPARSE_NODATA) {
			Grid3D_Regular_Sparse_NoData<T> *sparseNoData = (Grid3D_Regular_Sparse_NoData<T>*)grid;

			if (sparseNoData->hasNormalMapping()) matrixType = 1;
			else matrixType = 2;

			rowColIndexSize = sparseNoData->getRowColIndexSize();
			depthIndexSize = sparseNoData->getDepthIndexSize();
		}
		else if (type == GRID_REGULAR_BLOCK) {
			Grid3D_Regular_Block<T> *block = (Grid3D_Regular_Block<T>*)grid;

			matrixType = 3;
			getDataSize(&block->blockDataPointer, blockIndexSize, false);
		}
		else printf("unknown grid type in GridCuda::getMatrixTypesAndIndexOffsets()\n");
	}

	template <class T> inline void addToCSRIndexVectors(const Grid3D_Regular_Base<T> *grid, unsigned int *rowColIndices,
														unsigned int *depthIndices, bool forceSparseIfVectorType) {
		int type = grid->getType();

		if (type == GRID_REGULAR || type == GRID_POW2) {
			Grid3D_Regular<T> *full = (Grid3D_Regular<T>*)grid;

			if (forceSparseIfVectorType && TemplateType::isVector(full->outsideVal)) {
				Grid3D_Regular_CSR<T> *csr = new Grid3D_Regular_CSR<T>(*full, full->outsideVal);
				addToCSRIndexVectors(csr, rowColIndices, depthIndices, forceSparseIfVectorType);
				delete csr;  csr = NULL;
			}
			//else // no need to do anything since we are not copying the data in this function
		}
		else if (type == GRID_REGULAR_CSR || type == GRID_POW2_CSR) {
			Grid3D_Regular_CSR<T> *csr = (Grid3D_Regular_CSR<T>*)grid;

			vector<unsigned int> *rcIndices = csr->getRowColIndexPtr();
			vector<unsigned int> *dIndices = csr->getDepthIndexPtr();

			for (unsigned int j=0; j<rcIndices->size(); j++) rowColIndices[j] = (*rcIndices)[j];
			for (unsigned int j=0; j<dIndices->size(); j++) depthIndices[j] = (*dIndices)[j];
		}
		else if (type == GRID_REGULAR_SPARSE_NODATA || type == GRID_POW2_SPARSE_NODATA) {
			Grid3D_Regular_Sparse_NoData<T> *sparseNoData = (Grid3D_Regular_Sparse_NoData<T>*)grid;

			vector<unsigned int> *rcIndices = sparseNoData->getRowColIndexPtr();
			vector<unsigned int> *dIndices = sparseNoData->getDepthIndexPtr();

			for (unsigned int j=0; j<rcIndices->size(); j++) rowColIndices[j] = (*rcIndices)[j];
			for (unsigned int j=0; j<dIndices->size(); j++) depthIndices[j] = (*dIndices)[j];
		}
		else if (type == GRID_REGULAR_BLOCK || type == GRID_REGULAR_BLOCK_NODATA) { } // do nothing
		else printf("unknown grid type in GridCuda::addToCSRIndexVectors()\n");
	}

	template <class T> inline void addToBlockIndexVectors(const Grid3D_Regular_Base<T> *grid, unsigned int *blockIndices) {
		int type = grid->getType();

		if (type == GRID_REGULAR_BLOCK) {
			Grid3D_Regular_Block<T> *block = (Grid3D_Regular_Block<T>*)grid;
			addToDataVector(&block->blockDataPointer, blockIndices, false);
		}
		else if (type == GRID_REGULAR_BLOCK_NODATA) {
			Grid3D_Regular_Block_NoData<T> *block = (Grid3D_Regular_Block_NoData<T>*)grid;
			addToDataVector(&block->blockDataPointer, blockIndices, false);
		}
		else if (type == GRID_REGULAR || type == GRID_POW2 || type == GRID_REGULAR_CSR || type == GRID_POW2_CSR ||
				 type == GRID_REGULAR_SPARSE_NODATA || type == GRID_POW2_SPARSE_NODATA) { } // do nothing
		else printf("unknown grid type in GridCuda::addToBlockIndexVectors()\n");
	}

	// to do: merge with getMatrixTypeAndIndexSizes
	template <class T> inline void getDataSize(const Grid3D_Regular_Base<T> *grid, unsigned int &dataSize, bool forceSparseIfVectorType) {
		int type = grid->getType();

		dataSize = 0;

		if (type == GRID_REGULAR || type == GRID_POW2) {
			Grid3D_Regular<T> *full = (Grid3D_Regular<T>*)grid;

			if (forceSparseIfVectorType && TemplateType::isVector(full->outsideVal)) {
				Grid3D_Regular_CSR<T> *csr = new Grid3D_Regular_CSR<T>(*full, full->outsideVal);
				getDataSize(csr, dataSize, forceSparseIfVectorType);
				delete csr;  csr = NULL;
			}
			else dataSize = full->size();
		}
		else if (type == GRID_REGULAR_CSR || type == GRID_POW2_CSR) {
			Grid3D_Regular_CSR<T> *csr = (Grid3D_Regular_CSR<T>*)grid;
			dataSize = csr->getDepthIndexSize();
		}
		else if (type == GRID_REGULAR_SPARSE_NODATA || type == GRID_POW2_SPARSE_NODATA) {
			dataSize = 0;
		}
		else if (type == GRID_REGULAR_BLOCK) {
			Grid3D_Regular_Block<T> *block = (Grid3D_Regular_Block<T>*)grid;
			for (unsigned int i=0; i<block->blocks.size(); i++) {
				unsigned int size = 0;
				getDataSize(block->blocks[i], size, false);
				dataSize += size;
			}
		}
		else if (type == GRID_REGULAR_BLOCK_NODATA) {
			Grid3D_Regular_Block_NoData<T> *block = (Grid3D_Regular_Block_NoData<T>*)grid;
			for (unsigned int i=0; i<block->blocks.size(); i++) {
				unsigned int size = 0;
				getDataSize(block->blocks[i], size, false);
				dataSize += size;
			}
		}
		else printf("unknown grid type in GridCuda::getDataSize()\n");
	}

	template <class T> inline void addToDataVector(const Grid3D_Regular_Base<T> *grid, T *data, bool forceSparseIfVectorType) {
		int type = grid->getType();

		if (type == GRID_REGULAR || type == GRID_POW2) {
			Grid3D_Regular<T> *full = (Grid3D_Regular<T>*)grid;

			if (forceSparseIfVectorType && TemplateType::isVector(full->outsideVal)) {
				Grid3D_Regular_CSR<T> *csr = new Grid3D_Regular_CSR<T>(*full, full->outsideVal);
				addToDataVector(csr, data, forceSparseIfVectorType);
				delete csr;  csr = NULL;
			}
			else {
				unsigned int count = 0;
				/*for (unsigned int i=0; i<full->numCells.x; i++)
				for (unsigned int j=0; j<full->numCells.y; j++)
				for (unsigned int k=0; k<full->numCells.z; k++) {*/
				for (unsigned int n=0; n<full->size(); n++) {
					data[count] = full->get(n);
					count++;
				}
			}
		}
		else if (type == GRID_REGULAR_CSR || type == GRID_POW2_CSR) {
			Grid3D_Regular_CSR<T> *csr = (Grid3D_Regular_CSR<T>*)grid;
			vector<T> *dataPtr = csr->getDataPtr();

			for (unsigned int j=0; j<dataPtr->size(); j++) data[j] = (*dataPtr)[j];
		}
		else if (type == GRID_REGULAR_SPARSE_NODATA || type == GRID_POW2_SPARSE_NODATA) { }
		else if (type == GRID_REGULAR_BLOCK) {
			Grid3D_Regular_Block<T> *block = (Grid3D_Regular_Block<T>*)grid;
			unsigned int count = 0;
			for (unsigned int i=0; i<block->blocks.size(); i++) {
				unsigned int size = 0;
				getDataSize(block->blocks[i], size, false);
				
				T *blockData = new T[size];
				addToDataVector(block->blocks[i], blockData, false);

				for (unsigned int j=0; j<size; j++) data[count+j] = blockData[j];

				delete[] blockData;  blockData = NULL;
				count += size;
			}
		}
		else if (type == GRID_REGULAR_BLOCK_NODATA) {
			Grid3D_Regular_Block_NoData<T> *block = (Grid3D_Regular_Block_NoData<T>*)grid;
			printf("to do: need to figure out this guy\n");
		}
		else printf("unknown grid type in GridCuda::addToDataVector()\n");
	}

	// data[0] should already be initialized
	template <class T> inline void addToDataVector(const Grid3D_Regular_Base< vector<T> > *grid,
												   unsigned int *data, bool forceSparseIfVectorType) {
		int type = grid->getType();

		if (type == GRID_REGULAR || type == GRID_POW2) {
			Grid3D_Regular< vector<T> > *full = (Grid3D_Regular< vector<T> >*)grid;

			if (forceSparseIfVectorType && TemplateType::isVector(full->outsideVal)) {
				Grid3D_Regular_CSR< vector<T> > *csr = new Grid3D_Regular_CSR< vector<T> >(*full, full->outsideVal);
				addToDataVector(csr, /*offset,*/ data, forceSparseIfVectorType);
				delete csr;  csr = NULL;
			}
			else {
				unsigned int count = 0;
				/*for (unsigned int i=0; i<full->numCells.x; i++)
				for (unsigned int j=0; j<full->numCells.y; j++)
				for (unsigned int k=0; k<full->numCells.z; k++) {*/
				for (unsigned int i=0; i<full->size(); i++) {
					vector<T> val = full->get(i);
					data[count+1] = data[count] + (unsigned int)val.size();
					count++;
				}
			}
		}
		else if (type == GRID_REGULAR_CSR || type == GRID_POW2_CSR) {
			Grid3D_Regular_CSR< vector<T> > *csr = (Grid3D_Regular_CSR< vector<T> >*)grid;
			vector< vector<T> > *dataPtr = csr->getDataPtr();

			for (unsigned int i=0; i<dataPtr->size(); i++) {
				vector<T> val = (*dataPtr)[i];
				data[i+1] = data[i] + (unsigned int)val.size();
			}
		}
		else if (type == GRID_REGULAR_SPARSE_NODATA || type == GRID_POW2_SPARSE_NODATA) { }
		else if (type == GRID_REGULAR_BLOCK) {
			Grid3D_Regular_Block< vector<T> > *block = (Grid3D_Regular_Block< vector<T> >*)grid;
			unsigned int count = 0;

			for (unsigned int i=0; i<block->blocks.size(); i++) {
				unsigned int size = 0;
				getDataSize(block->blocks[i], size, false);
				
				if (size == 0) continue;

				for (unsigned int j=0; j<block->blocks[i]->size(); j++) {
					vector<T> val = block->blocks[i]->get(j);
					data[count+1] = data[count] + (unsigned int)val.size();
					count++;
				}
			}
		}
		else if (type == GRID_REGULAR_BLOCK_NODATA) {
			Grid3D_Regular_Block_NoData<T> *block = (Grid3D_Regular_Block_NoData<T>*)grid;
			printf("to do: need to figure out this guy\n");
		}
		else printf("unknown grid type in GridCuda::addToDataVector()\n");
	}

	template <class T> inline unsigned int addVectorsToDataVector(const Grid3D_Regular_Base<T> *grid,
																  T *vectorData, bool forceSparseIfVectorType) {
		return 0;
	}

	template <class T> inline unsigned int addVectorsToDataVector(const Grid3D_Regular_Base< vector<T> > *grid,
																  T *vectorData, bool forceSparseIfVectorType) {
		int type = grid->getType();
		unsigned int idx = 0;

		if (type == GRID_REGULAR || type == GRID_POW2) {
			Grid3D_Regular< vector<T> > *full = (Grid3D_Regular< vector<T> >*)grid;

			if (forceSparseIfVectorType && TemplateType::isVector(full->outsideVal)) {
				Grid3D_Regular_CSR< vector<T> > *csr = new Grid3D_Regular_CSR< vector<T> >(*full, full->outsideVal);
				return addVectorsToDataVector(csr, vectorData, forceSparseIfVectorType);
				delete csr;  csr = NULL;
			}
			else {
				/*for (unsigned int i=0; i<full->numCells.x; i++)
				for (unsigned int j=0; j<full->numCells.y; j++)
				for (unsigned int k=0; k<full->numCells.z; k++) {*/
				for (unsigned int i=0; i<full->size(); i++) {
					vector<T> val = full->get(i);

					for (unsigned int j=0; j<val.size(); j++) {
						vectorData[idx] = val[j];
						idx++;
					}
				}
			}
		}
		else if (type == GRID_REGULAR_CSR || type == GRID_POW2_CSR) {
			Grid3D_Regular_CSR< vector<T> > *csr = (Grid3D_Regular_CSR< vector<T> >*)grid;
			vector< vector<T> > *dataPtr = csr->getDataPtr();

			for (unsigned int i=0; i<dataPtr->size(); i++) {
				vector<T> val = (*dataPtr)[i];

				for (unsigned int j=0; j<val.size(); j++) {
					vectorData[idx] = val[j];
					idx++;
				}
			}
		}
		else if (type == GRID_REGULAR_SPARSE_NODATA || type == GRID_POW2_SPARSE_NODATA) { }
		else if (type == GRID_REGULAR_BLOCK) {
			Grid3D_Regular_Block< vector<T> > *block = (Grid3D_Regular_Block< vector<T> >*)grid;
			unsigned int count = 0;
			for (unsigned int i=0; i<block->blocks.size(); i++) {
				unsigned int size = 0;
				getDataSize(block->blocks[i], size, false);

				Grid3D_Regular< vector<T> > *full = (Grid3D_Regular< vector<T> >*)block->blocks[i];

				for (unsigned int i=0; i<full->size(); i++) {
					vector<T> *val = full->getPtr(i);

					for (unsigned int j=0; j<val->size(); j++) {
						vectorData[count] = (*val)[j];
						count++;
					}
				}
			}

			idx = count;
		}
		else if (type == GRID_REGULAR_BLOCK_NODATA) {
			Grid3D_Regular_Block_NoData<T> *block = (Grid3D_Regular_Block_NoData<T>*)grid;
			printf("to do: need to figure out this guy\n");
		}
		else printf("unknown grid type in GridCuda::addVectorsToDataVector()\n");

		return idx;
	}
}

#endif