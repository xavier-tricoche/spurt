#ifndef GRID_3D_COLLECTION_H
#define GRID_3D_COLLECTION_H

#include "GridUtils.h"

#include "StlVectorUtils.h"
#include "TemplateHelper.h"

template <class T> class Grid3D_Collection {
public:
	Grid3D_Collection() { }

	vector<Grid3D_Regular_Base<T>*> grids;

	~Grid3D_Collection() { clear(); }

	inline void clear() {
		for (unsigned int i=0; i<grids.size(); i++) {
			delete grids[i];  grids[i] = NULL;
		}
		grids.clear();
	}

	inline unsigned int size() const { return (unsigned int)grids.size(); }

	inline Grid3D_Regular_Base<T>*& operator[](int idx) { return grids[idx]; }

	inline Grid3D_Regular_Base<T>* operator[](int idx) const { return grids[idx]; }
	inline Grid3D_Regular_Base<T>* getFirst() const { return grids[0]; }
	inline Grid3D_Regular_Base<T>* getLast() const { return grids[grids.size()-1]; }

	inline void add(Grid3D_Regular_Base<T> *gridPtr) { grids.push_back(gridPtr); }
	inline void add(int type) { add(GridCreation::create<T>(type)); }
	template <class S> inline void add(Grid3D_Regular_Base<S> *gridPtr, const T outsideVal, int type) {
		add(type);
		Grid3D_Regular_Base<T>* grid = getLast();
		grid->rebuildGrid(gridPtr->numCells);
		grid->setGridPosition(gridPtr->startPos, gridPtr->cellSize);
	}

	inline void add_overlapPrevious(const int type, const Vector3ui numCellsFactor) {
		Grid3D_Regular_Base<T> *prev = getLast();

		add(type);

		Vector3ui numCells = prev->numCells * numCellsFactor;
		getLast()->rebuildGrid(numCells);

		Vector3d cellSize = prev->cellSize * Vector3d(prev->numCells) / Vector3d(numCells);
		Vector3d startPos = prev->startPos;
		getLast()->setGridPosition(startPos, cellSize);

		getLast()->outsideVal = prev->outsideVal;
	}

	inline void flatten() { // works when cell dimensions are some factor of each other
		for (unsigned int n=0; n<grids.size()-1; n++) {
			T nullVal = grids[n]->outsideVal;
			Vector3ui factor((grids[grids.size()-1]->numCells.x - 1) / (grids[n]->numCells.x - 1),
							 (grids[grids.size()-1]->numCells.y - 1) / (grids[n]->numCells.y - 1),
							 (grids[grids.size()-1]->numCells.z - 1) / (grids[n]->numCells.z - 1));

			// copy over the data
			for (unsigned int i=0; i<grids[n]->numCells.x; i++)
			for (unsigned int j=0; j<grids[n]->numCells.y; j++)
			for (unsigned int k=0; k<grids[n]->numCells.z; k++) {
				T val = grids[n]->get(i,j,k);

				if (val != nullVal) {
					Vector3ui idx = factor * Vector3ui(i,j,k);
					grids[grids.size()-1]->set(idx.x, idx.y, idx.z, val);
				}
			}

			// remove matrix
			delete grids[n];  grids[n] = NULL;
		}

		Grid3D_Regular_Base<T> *bottomPtr = grids[grids.size()-1];
		grids.clear();
		grids.push_back(bottomPtr);
	}

	inline double getSparsityRatio(const unsigned int index) const { return GridSize::calculateSparsityRatio(grids[index]); }

	inline void convertToNonDataRep(const unsigned int index) {
		int convertType = 0;
		if (grids[index]->getType() == GRID_POW2_CSR) convertType = GRID_POW2_SPARSE_NODATA;
		else if (grids[index]->getType() == GRID_REGULAR_CSR) convertType = GRID_REGULAR_SPARSE_NODATA;
		else if (grids[index]->getType() == GRID_REGULAR_BLOCK) convertType = GRID_REGULAR_BLOCK_NODATA;
		else return;

		Grid3D_Regular_Base<T> *newGrid = GridCreation::createFrom(convertType, grids[index], true);
		delete grids[index];
		grids[index] = newGrid;
	}

	inline void convertToMoreEfficientRep(const unsigned int index, bool useSpaceFillingCurve, bool removeData, bool forceFull) {
		if (index >= grids.size()) return;
		
		printf(" - %d converted to ", index);

		double sparsityRatio = GridSize::calculateSparsityRatio(grids[index]);

		int convertType;
		if (forceFull) {
			if (useSpaceFillingCurve) convertType = GRID_POW2;
			else convertType = GRID_REGULAR;
		}
		else if (removeData) {
			if (useSpaceFillingCurve) convertType = GRID_POW2_SPARSE_NODATA;
			else convertType = GRID_REGULAR_SPARSE_NODATA;
		}
		else if (!GridConversion::canBeRepresentedAsFullMatrix(grids[index]->numCells)) {
			if (useSpaceFillingCurve) convertType = GRID_POW2_CSR;
			else convertType = GRID_REGULAR_CSR;
		}
		else convertType = GridConversion::calculateMostSpaceEfficientRep(grids[index], sparsityRatio, useSpaceFillingCurve);

		if (convertType != grids[index]->getType()) {
			Grid3D_Regular_Base<T> *newGrid = GridCreation::createFrom(convertType, grids[index], true);
			delete grids[index];
			grids[index] = newGrid;
		}

		printf("%s (sparsity = %f)\n", GridString::get(convertType).c_str(), GridSize::getActualSparsityRatio(grids[index], sparsityRatio));
	}

	inline void convertToMoreEfficientReps(bool useSpaceFillingCurve, bool removeData, bool forceFull) {
		printf("Converting %d matrices...\n", grids.size());
		for (unsigned int i=0; i<grids.size(); i++) {
			convertToMoreEfficientRep(i, useSpaceFillingCurve, removeData, forceFull);
		}
	}

	inline void write(FILE *file) const {
		unsigned int size = (unsigned int)grids.size();
		fwrite(&size, sizeof(unsigned int), 1, file);

		for (unsigned int n=0; n<grids.size(); n++) {
			grids[n]->writeType(file);
			grids[n]->writeGridProperties(file);
			grids[n]->writeOutsideValue(file);
			grids[n]->writeData(file);
		}
	}

	inline void read(FILE *file) {
		clear();

		unsigned int size = 0;
		fread(&size, sizeof(unsigned int), 1, file);

		for (unsigned int n=0; n<size; n++) {
			int type;  fread(&type, sizeof(int), 1, file);
			
			add(type);

			grids[n]->readGridProperties(file);
			grids[n]->readOutsideValue(file);
			grids[n]->readData(file);
		}
	}

	inline void getCudaRep_BasicData(int **matrixTypes, unsigned int **rowColIndexOffsets, unsigned int **depthIndexOffsets,
									 unsigned int **blockIndexOffsets, unsigned int **numberOfCells, float **cellSizes,
									 bool forceSparseIfVectorType) {
		if (grids.size() == 0) return;

		*matrixTypes = new int[grids.size()];
		*rowColIndexOffsets = new unsigned int[grids.size()+1];
		*depthIndexOffsets = new unsigned int[grids.size()+1];
		*blockIndexOffsets = new unsigned int [grids.size()+1];
		*numberOfCells = new unsigned int[grids.size()*3];
		*cellSizes = new float[grids.size()*3];

		// figure out the offsets
		(*rowColIndexOffsets)[0] = (*depthIndexOffsets)[0] = (*blockIndexOffsets)[0] = 0;

		for (unsigned int n=0; n<grids.size(); n++) {
			(*numberOfCells)[ n*3 ] = grids[n]->numCells.x;
			(*numberOfCells)[n*3+1] = grids[n]->numCells.y;
			(*numberOfCells)[n*3+2] = grids[n]->numCells.z;

			(*cellSizes)[ n*3 ] = grids[n]->cellSize.x;
			(*cellSizes)[n*3+1] = grids[n]->cellSize.y;
			(*cellSizes)[n*3+2] = grids[n]->cellSize.z;

			int type = grids[n]->getType();

			unsigned int rowColIndexSize = 0, depthIndexSize  = 0, blockIndexSize = 0;
			GridCuda::getMatrixTypeAndIndexSizes(grids[n], (*matrixTypes)[n], rowColIndexSize, depthIndexSize, blockIndexSize, forceSparseIfVectorType);

			(*rowColIndexOffsets)[n+1] = (*rowColIndexOffsets)[n] + rowColIndexSize;
			(*depthIndexOffsets)[n+1] = (*depthIndexOffsets)[n] + depthIndexSize;
			(*blockIndexOffsets)[n+1] = (*blockIndexOffsets)[n] + blockIndexSize;
		}
	}

	inline void getCudaRep_CSRIndexOffsets(unsigned int **rowColIndices, unsigned int **depthIndices, unsigned int *rowColIndexOffsets,
										   unsigned int *depthIndexOffsets, bool forceSparseIfVectorType) {
		// copy over the indices
		*rowColIndices = *depthIndices = NULL;

		if (rowColIndexOffsets[grids.size()] == 0 && depthIndexOffsets[grids.size()] == 0) return;

		if (rowColIndexOffsets[grids.size()] != 0) *rowColIndices = new unsigned int[rowColIndexOffsets[grids.size()]];
		if (depthIndexOffsets[grids.size()] != 0) *depthIndices = new unsigned int[depthIndexOffsets[grids.size()]];

		for (unsigned int n=0; n<grids.size(); n++) {
			GridCuda::addToCSRIndexVectors(grids[n], (*rowColIndices)+rowColIndexOffsets[n],
										   (*depthIndices)+depthIndexOffsets[n], forceSparseIfVectorType);
		}
	}

	inline void getCudaRep_BlockIndexOffsets(unsigned int **blockIndices, unsigned int *blockIndexOffsets) {
		// copy over the indices
		*blockIndices = NULL;

		if (blockIndexOffsets[grids.size()] == 0) return;

		*blockIndices = new unsigned int[blockIndexOffsets[grids.size()]];

		for (unsigned int n=0; n<grids.size(); n++) {
			GridCuda::addToBlockIndexVectors(grids[n], (*blockIndices)+blockIndexOffsets[n]);
		}
	}

	inline void getCudaRep_DataOffsets(unsigned int **dataOffsets, bool forceSparseIfVectorType) {
		*dataOffsets = new unsigned int[grids.size()+1];
		(*dataOffsets)[0] = 0;

		for (unsigned int n=0; n<grids.size(); n++) {
			unsigned int dataIndexSize;
			GridCuda::getDataSize(grids[n], dataIndexSize, forceSparseIfVectorType);
			(*dataOffsets)[n+1] = (*dataOffsets)[n] + dataIndexSize;
		}
	}

	// S == T when T != vector<>, else S == unsigned int (i.e. it will get the offsets to the vectors)
	template <class S> inline void getCudaRep_Data(S **data, unsigned int *dataOffsets, bool forceSparseIfVectorType) {
		*data = NULL;
		if (dataOffsets[grids.size()] > 0) {
			dataOffsets[grids.size()]++;

			*data = new S[dataOffsets[grids.size()]];
			(*data)[0] = 0;

			for (unsigned int n=0; n<grids.size(); n++) {
				GridCuda::addToDataVector(grids[n], (*data)+dataOffsets[n], forceSparseIfVectorType);
			}
		}
	}

	// should only be called when T == vector<S>
	template <class S> inline void getCudaRep_VectorData(S **vectorData, unsigned int totalData, bool forceSparseIfVectorType) {
		*vectorData = NULL;
		if (totalData > 0) {
			*vectorData = new S[totalData];
			unsigned int count = 0;

			for (unsigned int n=0; n<grids.size(); n++) {
				unsigned int numberAdded = GridCuda::addVectorsToDataVector(grids[n], (*vectorData)+count, forceSparseIfVectorType);
				count += numberAdded;
			}

			if (count != totalData) {
				printf("problem in getCudaRep_VectorData\n");
			}
		}
	}



private:

};

#endif