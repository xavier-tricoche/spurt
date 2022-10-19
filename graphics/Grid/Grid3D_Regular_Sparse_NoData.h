#ifndef GRID_3D_REGULAR_SPARSE_NODATA_H
#define GRID_3D_REGULAR_SPARSE_NODATA_H

#include "Grid3D_Regular_Base.h"
#include "Grid3D_Regular_Sparse_I1.h"
#include "Grid3D_Regular.h"

#include "StlVectorUtils.h"

template <class T> class Grid3D_Regular_Sparse_NoData : public Grid3D_Regular_Base<T> {
public:
	Grid3D_Regular_Sparse_NoData() { }

	Grid3D_Regular_Sparse_NoData(Vector3ui numberOfCells, T outsideValue) {
		outsideVal = outsideValue;
		rebuildGrid(numberOfCells);
		setGridPosition(Vector3d(0,0,0), Vector3d(1,1,1));
	}

	Grid3D_Regular_Sparse_NoData(Vector3ui numberOfCells, T outsideValue, Vector3d startPosition, Vector3d cellSizes) {
		outsideVal = outsideValue;
		rebuildGrid(numberOfCells);
		setGridPosition(startPosition, cellSizes);
	}

	Grid3D_Regular_Sparse_NoData(const Grid3D_Regular<T> &reg, const T &nullValue) {
		outsideVal = nullValue;
		rebuildGrid(reg.numCells);
		setGridPosition(reg.startPos, reg.cellSize);
		copyDataFrom(reg, nullValue);
	}

	Grid3D_Regular_Sparse_NoData(const Grid3D_Regular_Sparse_I1<T> &unfixedRep) {
		outsideVal = unfixedRep.outsideVal;
		rebuildGrid(unfixedRep.numCells);
		setGridPosition(unfixedRep.startPos, unfixedRep.cellSize);
		copyDataFrom(unfixedRep);
	}

	~Grid3D_Regular_Sparse_NoData() { clear(); }

	inline void clear() {
		rowColIndices.clear();
		depthIndices.clear();
	}

	inline void deletePointers() { }

	inline int getType() const { return GRID_REGULAR_SPARSE_NODATA; }

	void rebuildGrid(Vector3ui numberOfCells) {
		Grid3D_Regular_Base::rebuildGrid(numberOfCells);

		rowColIndices.resize(numberOfCells.x*numberOfCells.y+1);
		for (unsigned int i=0; i<rowColIndices.size(); i++) rowColIndices.vec[i] = i;

		depthIndices.resize(numberOfCells.x*numberOfCells.y);
		for (unsigned int i=0; i<depthIndices.size(); i++) depthIndices.vec[i] = i;
	}

	inline void set(const unsigned int x, const unsigned int y, const unsigned int z, const T &val) {
		unsigned int depth, xyIndex = rowColIndexMap(x,y);
		if (depthIndices.binarySearch_sortedLowToHigh(z, rowColIndices[xyIndex], rowColIndices[xyIndex+1]-1, &depth)) { }
		else {
			depthIndices.insert(depth, z);

			for (unsigned int i=xyIndex+1; i<rowColIndices.size(); i++) rowColIndices.vec[i]++;
		}
	}

	inline T set_returnOld(const unsigned int x, const unsigned int y, const unsigned int z, const T &val) {
		unsigned int depth, xyIndex = rowColIndexMap(x,y);
		if (depthIndices.binarySearch_sortedLowToHigh(z, rowColIndices[xyIndex], rowColIndices[xyIndex+1]-1, &depth)) {
			return outsideVal;
		}
		else {
			depthIndices.insert(depth, z);

			for (unsigned int i=xyIndex+1; i<rowColIndices.size(); i++) rowColIndices.vec[i]++;

			return outsideVal;
		}
	}

	inline T get(const unsigned int x, const unsigned int y, const unsigned int z) const { return outsideVal; }
	inline T* getPtr(const unsigned int x, const unsigned int y, const unsigned int z) { return NULL; }

	inline bool entryExists(const unsigned int x, const unsigned int y, const unsigned int z) const {
		unsigned int depth, xyIndex = rowColIndexMap(x,y);
		if (depthIndices.binarySearch_sortedLowToHigh(z, rowColIndices[xyIndex], rowColIndices[xyIndex+1]-1, &depth)) {
			return normalMapping;
		}
		return !normalMapping;
	}

	inline unsigned int getRowColIndexSize() const { return rowColIndices.size(); }
	inline unsigned int getDepthIndexSize() const { return depthIndices.size(); }

	inline vector<unsigned int>* getRowColIndexPtr() { return &rowColIndices.vec; }
	inline vector<unsigned int>* getDepthIndexPtr() { return &depthIndices.vec; }

	inline bool hasNormalMapping() const { return normalMapping; }
	inline void setNormalMapping(const bool isNormalMapping) { normalMapping = isNormalMapping; }

	inline double getNumberOfBytes() const { return sizeof(rowColIndices) + sizeof(depthIndices) + rowColIndices.size()*sizeof(unsigned int) + depthIndices.size()*sizeof(unsigned int); }

	inline void writeData(FILE *file) const {
		TemplateFileIO::writeBinary(file, normalMapping);
		TemplateFileIO::writeBinary(file, rowColIndices);
		TemplateFileIO::writeBinary(file, depthIndices);
	}
	inline void readData(FILE *file) {
		TemplateFileIO::readBinary(file, &normalMapping);
		TemplateFileIO::readBinary(file, &rowColIndices);
		TemplateFileIO::readBinary(file, &depthIndices);
	}

	void copyDataFrom(const Grid3D_Regular<T> &reg, const T &nullValue) {
		double numberOfActualEntries = 0;
		for (unsigned int n=0; n<reg.size(); n++) {
			if (reg.get(n) != nullValue) numberOfActualEntries++;
		}

		if (2*numberOfActualEntries > numCells.x*numCells.y*numCells.z) normalMapping = false;
		else normalMapping = true;

		rowColIndices.vec[0] = 0;

		for (unsigned int n=0; n<rowColIndices.size()-1; n++) {
			unsigned int x, y;
			getReversedRowColMapping(n, x, y);

			unsigned int rowColSize = 0;

			if (x < reg.numCells.x && y < reg.numCells.y) {
				for (unsigned int z=0; z<reg.numCells.z; z++) {
					if (reg.get(x,y,z) != nullValue) rowColSize++;
				}
			}

			if (!normalMapping) rowColSize = reg.numCells.z-rowColSize;

			rowColIndices.vec[n+1] = rowColIndices.vec[n] + rowColSize;
		}

		depthIndices.resize(rowColIndices[rowColIndices.size()-1]);

		unsigned int count = 0;
		for (unsigned int n=0; n<rowColIndices.size()-1; n++) {
			unsigned int x, y;
			getReversedRowColMapping(n, x, y);

			if (x < reg.numCells.x && y < reg.numCells.y) {
				for (unsigned int z=0; z<reg.numCells.z; z++) {
					bool add = false;
					if (reg.get(x,y,z) != nullValue) {
						if (normalMapping) add = true;
					}
					else if (!normalMapping) add = true;

					if (add) {
						depthIndices.vec[count] = z;
						count++;
					}
				}
			}
		}
	}

	void copyDataFrom(const Grid3D_Regular_Sparse_I1<T> &unfixedRep) {
		if (2*unfixedRep.numberOfEntries > numCells.x*numCells.y*numCells.z) normalMapping = false;
		else normalMapping = true;

		unsigned int n = rowColIndices.size()-1;

		rowColIndices.vec[0] = 0;
		for (unsigned int i=0; i<n; i++) {
			unsigned int x, y;
			getReversedRowColMapping(i, x, y);

			unsigned int rowColSize = 0;
			if (x < unfixedRep.numCells.x && y < unfixedRep.numCells.y) {
				rowColSize = unfixedRep.depthIndices[x][y].size();
				if (!normalMapping) rowColSize = unfixedRep.numCells.z-rowColSize;
			}

			rowColIndices.vec[i+1] = rowColIndices[i] + rowColSize;
		}

		depthIndices.resize(rowColIndices[n]);

		unsigned int count = 0;
		for (unsigned int i=0; i<n; i++) {
			unsigned int x, y;
			getReversedRowColMapping(i, x, y);

			if (x < unfixedRep.numCells.x && y < unfixedRep.numCells.y) {
				if (normalMapping) {
					for (unsigned int k=0; k<unfixedRep.depthIndices[x][y].size(); k++) {
						depthIndices.vec[count] = unfixedRep.depthIndices[x][y][k];
						count++;
					}
				}
				else {
					// go through all the values and add in the indices to the "zeros"
					for (unsigned int k=0; k<unfixedRep.depthIndices[x][y].size(); k++) {
						unsigned int startVal = 0;
						if (k != 0) startVal = unfixedRep.depthIndices[x][y][k-1]+1;

						for (unsigned int n=startVal; n<unfixedRep.depthIndices[x][y][k]; n++) {
							depthIndices.vec[count] = n;
							count++;
						}

						if (k == unfixedRep.depthIndices[x][y].size()-1) {
							for (unsigned int n=unfixedRep.depthIndices[x][y][k]+1; n<numCells.z; n++) {
								depthIndices.vec[count] = n;
								count++;
							}
						}		
					}
				}
			}
		}
	}

	inline bool getDataIndex(const unsigned int x, const unsigned int y, const unsigned int z, unsigned int &index) const {
		unsigned int xyIndex = rowColIndexMap(x,y);
		if (depthIndices.binarySearch_sortedLowToHigh(z, rowColIndices[xyIndex], rowColIndices[xyIndex+1]-1, &index)) {
			return normalMapping;
		}
		else {
			// depth = index of first value that is less than z
			// depth indices in prior row/columns = numCells.z*xyIndex - rowColIndices[xyIndex]
			// depth indices in this row/column before z = z - depth + rowColIndices[xyIndex]
			index = numCells.z*xyIndex + z - index;
			return !normalMapping;
		}
	}

	inline unsigned int getDataSize() const {
		unsigned int dataSize = depthIndices.size();
		if (!normalMapping) dataSize = numCells.x*numCells.y*numCells.z - dataSize;
		return dataSize;
	}

protected:
	StlVector<unsigned int> rowColIndices;
	StlVector<unsigned int> depthIndices;
	bool normalMapping; // if false, keep track of "zeros" instead those places with data

	virtual inline unsigned int rowColIndexMap(const unsigned int x, const unsigned int y) const { return x*numCells.y + y; }
	virtual inline void getReversedRowColMapping(const unsigned int index, unsigned int &x, unsigned int &y) {
		x = index * numCellsInv.y;
		y = index % numCells.y;
	}
};

template <class T> class Grid3D_Pow2_Sparse_NoData : public Grid3D_Regular_Sparse_NoData<T> {
public:
	Grid3D_Pow2_Sparse_NoData() { }

	Grid3D_Pow2_Sparse_NoData(Vector3ui numberOfCells, T outsideValue)
		{ Grid3D_Regular_Sparse_NoData::Grid3D_Regular_Sparse_NoData(getPowerOf2Indices(numberOfCells), outsideValue); }

	Grid3D_Pow2_Sparse_NoData(Vector3ui numberOfCells, T outsideValue, Vector3d startPosition, Vector3d cellSizes)
		{ Grid3D_Regular_Sparse_NoData::Grid3D_Regular_Sparse_NoData(getPowerOf2Indices(numberOfCells), outsideValue, startPosition, cellSizes); }

	Grid3D_Pow2_Sparse_NoData(const Grid3D_Regular<T> &reg, const T &nullValue) {
		outsideVal = nullValue;
		rebuildGrid(reg.numCells);
		setGridPosition(reg.startPos, reg.cellSize);
		copyDataFrom(reg, nullValue);
	}

	Grid3D_Pow2_Sparse_NoData(const Grid3D_Regular_Sparse_I1<T> &unfixedRep) {
		outsideVal = unfixedRep.outsideVal;
		rebuildGrid(unfixedRep.numCells);
		setGridPosition(unfixedRep.startPos, unfixedRep.cellSize);
		copyDataFrom(unfixedRep);
	}


	void rebuildGrid(Vector3ui numberOfCells) { Grid3D_Regular_Sparse_NoData::rebuildGrid(getPowerOf2Indices(numberOfCells)); }

	inline int getType() { return GRID_POW2_SPARSE_NODATA; };

	Vector3ui initialNumCells;

	inline unsigned int rowColIndexMap(const unsigned int x, const unsigned int y) const { return MathUtils::interleaveBits(x,y); }
	inline void getReversedRowColMapping(const unsigned int index, unsigned int &x, unsigned int &y) { MathUtils::unInterleaveBits(index,x,y); }

protected:
	inline Vector3ui getPowerOf2Indices(const Vector3ui originalNumberOfCells) {
		initialNumCells = originalNumberOfCells;
		return Vector3ui(MathUtils::roundUpPowerOf2(originalNumberOfCells.x), MathUtils::roundUpPowerOf2(originalNumberOfCells.y), MathUtils::roundUpPowerOf2(originalNumberOfCells.z));
	}
};

#endif