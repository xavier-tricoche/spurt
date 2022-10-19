#ifndef GRID_3D_REGULAR_CSR_H
#define GRID_3D_REGULAR_CSR_H

#include "Grid3D_Regular_Base.h"
#include "Grid3D_Regular_Sparse_I1.h"
#include "Grid3D_Regular.h"

#include "StlVectorUtils.h"

template <class T> class Grid3D_Regular_CSR : public Grid3D_Regular_Base<T> {
public:
	Grid3D_Regular_CSR() { }

	Grid3D_Regular_CSR(Vector3ui numberOfCells, const T &outsideValue) {
		outsideVal = outsideValue;
		rebuildGrid(numberOfCells);
		setGridPosition(Vector3d(0,0,0), Vector3d(1,1,1));
	}

	Grid3D_Regular_CSR(Vector3ui numberOfCells, const T &outsideValue, Vector3d startPosition, Vector3d cellSizes) {
		outsideVal = outsideValue;
		rebuildGrid(numberOfCells);
		setGridPosition(startPosition, cellSizes);
	}

	Grid3D_Regular_CSR(const Grid3D_Regular<T> &reg, const T &nullValue) {
		outsideVal = nullValue;
		rebuildGrid(reg.numCells);
		setGridPosition(reg.startPos, reg.cellSize);
		copyDataFrom(reg, nullValue);
	}

	Grid3D_Regular_CSR(const Grid3D_Regular_Sparse_I1<T> &unfixedRep) {
		outsideVal = unfixedRep.outsideVal;
		rebuildGrid(unfixedRep.numCells);
		setGridPosition(unfixedRep.startPos, unfixedRep.cellSize);
		copyDataFrom(unfixedRep);
	}

	~Grid3D_Regular_CSR() { clear(); }

	inline void clear() {
		if (data.size() == 0) return;

		rowColIndices.clear();
		depthIndices.clear();
		data.deleteElements();
		data.clear();
	}

	inline void deletePointers() { data.deleteElements(); }

	inline int getType() const { return GRID_REGULAR_CSR; }

	virtual void rebuildGrid(Vector3ui numberOfCells) {
		Grid3D_Regular_Base::rebuildGrid(numberOfCells);

		rowColIndices.resize(numberOfCells.x*numberOfCells.y+1);
		for (unsigned int i=0; i<rowColIndices.size(); i++) rowColIndices.vec[i] = i; // to do: I don't like this, should be zero

		depthIndices.resize(numberOfCells.x*numberOfCells.y);
		for (unsigned int i=0; i<depthIndices.size(); i++) depthIndices.vec[i] = i;

		data.resize(numberOfCells.x*numberOfCells.y, outsideVal);
	}

	inline void set(const unsigned int x, const unsigned int y, const unsigned int z, const T &val) {
		unsigned int depth, xyIndex = rowColIndexMap(x,y);
		if (depthIndices.binarySearch_sortedLowToHigh(z, rowColIndices[xyIndex], rowColIndices[xyIndex+1]-1, &depth)) {
			data.vec[depth] = val;
		}
		else {
			depthIndices.insert(depth, z);
			data.insert(depth, val);

			for (unsigned int i=xyIndex+1; i<rowColIndices.size(); i++) rowColIndices.vec[i]++;
		}
	}

	inline T set_returnOld(const unsigned int x, const unsigned int y, const unsigned int z, const T &val) {
		unsigned int depth, xyIndex = rowColIndexMap(x,y);
		if (depthIndices.binarySearch_sortedLowToHigh(z, rowColIndices[xyIndex], rowColIndices[xyIndex+1]-1, &depth)) {
			T old = data[depth];  data.vec[depth] = val;  return old;
		}
		else {
			depthIndices.insert(depth, z);
			data.insert(depth, val);

			for (unsigned int i=xyIndex+1; i<rowColIndices.size(); i++) rowColIndices.vec[i]++;

			return outsideVal;
		}
	}

	inline T get(const unsigned int x, const unsigned int y, const unsigned int z) const {
		unsigned int depth, xyIndex = rowColIndexMap(x,y);
		if (depthIndices.binarySearch_sortedLowToHigh(z, rowColIndices[xyIndex], rowColIndices[xyIndex+1]-1, &depth)) {
			return data.vec[depth];
		}
		return outsideVal;
	}

	inline T* getPtr(const unsigned int x, const unsigned int y, const unsigned int z) {
		unsigned int depth, xyIndex = rowColIndexMap(x,y);
		if (depthIndices.binarySearch_sortedLowToHigh(z, rowColIndices[xyIndex], rowColIndices[xyIndex+1]-1, &depth)) {
			return &data.vec[depth];
		}
		return NULL;
	}

	inline void get2x2x2Block(const unsigned int lowX, const unsigned int lowY, const unsigned int lowZ,
							T &val000, T &val001, T &val010, T &val011, T &val100, T &val101, T &val110, T &val111) const {
		unsigned int depth, xyIndex = rowColIndexMap(lowX,lowY);
		if (depthIndices.binarySearch_sortedLowToHigh(lowZ, rowColIndices[xyIndex], rowColIndices[xyIndex+1]-1, &depth)) {
			val000 = data.vec[depth];
		}
		else val000 = outsideVal;

		// since we are fetching a 2x2 block, we just need to check the next depth index
		//   this will work even if the above value does not exist since the binary search returns the index of the
		//   first value that is less than what we searched for
		if (depth+1 < rowColIndices[xyIndex+1] && depthIndices[depth+1] == lowZ+1) val001 = data.vec[depth+1];
		else val001 = outsideVal;


		// chances are the position of what we just found will be close to the other rows, so use it to speed up the search
		//   worst case, this will basically increase the numbers of iterations of the binary search by 1
		unsigned int offset = depth-rowColIndices[xyIndex];

		xyIndex = getRowColIndexMapUp(xyIndex);
		unsigned int dim = rowColIndices[xyIndex+1]-rowColIndices[xyIndex];

		bool found;
		if (offset < dim) {
			depth = rowColIndices[xyIndex]+offset;
			if (lowZ == depthIndices[depth]) found = true;
			else if (lowZ < depthIndices[depth]) found = depthIndices.binarySearch_sortedLowToHigh(lowZ, rowColIndices[xyIndex], depth, &depth);
			else/*(lowZ > depthIndices[depth])*/ found = depthIndices.binarySearch_sortedLowToHigh(lowZ, depth, rowColIndices[xyIndex+1]-1, &depth);
		}
		else found = depthIndices.binarySearch_sortedLowToHigh(lowZ, rowColIndices[xyIndex], rowColIndices[xyIndex+1]-1, &depth);
		
		if (found) val010 = data.vec[depth];
		else val010 = outsideVal;

		if (depth+1 < rowColIndices[xyIndex+1] && depthIndices[depth+1] == lowZ+1) val011 = data.vec[depth+1];
		else val011 = outsideVal;



		offset = depth-rowColIndices[xyIndex];
		xyIndex = getRowColIndexMapRight(xyIndex);
		dim = rowColIndices[xyIndex+1]-rowColIndices[xyIndex];

		if (offset < dim) {
			depth = rowColIndices[xyIndex]+offset;
			if (lowZ == depthIndices[depth]) found = true;
			else if (lowZ < depthIndices[depth]) found = depthIndices.binarySearch_sortedLowToHigh(lowZ, rowColIndices[xyIndex], depth, &depth);
			else/*(lowZ > depthIndices[depth])*/ found = depthIndices.binarySearch_sortedLowToHigh(lowZ, depth, rowColIndices[xyIndex+1]-1, &depth);
		}
		else found = depthIndices.binarySearch_sortedLowToHigh(lowZ, rowColIndices[xyIndex], rowColIndices[xyIndex+1]-1, &depth);
		
		if (found) val110 = data.vec[depth];
		else val110 = outsideVal;

		if (depth+1 < rowColIndices[xyIndex+1] && depthIndices[depth+1] == lowZ+1) val111 = data.vec[depth+1];
		else val111 = outsideVal;



		offset = depth-rowColIndices[xyIndex];
		xyIndex = getRowColIndexMapDown(xyIndex);
		dim = rowColIndices[xyIndex+1]-rowColIndices[xyIndex];
		
		if (offset < dim) {
			depth = rowColIndices[xyIndex]+offset;
			if (lowZ == depthIndices[depth]) found = true;
			else if (lowZ < depthIndices[depth]) found = depthIndices.binarySearch_sortedLowToHigh(lowZ, rowColIndices[xyIndex], depth, &depth);
			else/*(lowZ > depthIndices[depth])*/ found = depthIndices.binarySearch_sortedLowToHigh(lowZ, depth, rowColIndices[xyIndex+1]-1, &depth);
		}
		else found = depthIndices.binarySearch_sortedLowToHigh(lowZ, rowColIndices[xyIndex], rowColIndices[xyIndex+1]-1, &depth);
		
		if (found) val100 = data.vec[depth];
		else val100 = outsideVal;

		if (depth+1 < rowColIndices[xyIndex+1] && depthIndices[depth+1] == lowZ+1) val101 = data.vec[depth+1];
		else val101 = outsideVal;
	}

	inline bool entryExists(const unsigned int x, const unsigned int y, const unsigned int z) const {
		unsigned int depth, xyIndex = rowColIndexMap(x,y);
		if (depthIndices.binarySearch_sortedLowToHigh(z, rowColIndices[xyIndex], rowColIndices[xyIndex+1]-1, &depth)) {
			return true;
		}
		return false;
	}

	inline bool entryExists(const unsigned int x, const unsigned int y, const unsigned int z, T **dataPtr) {
		unsigned int depth, xyIndex = rowColIndexMap(x,y);
		if (depthIndices.binarySearch_sortedLowToHigh(z, rowColIndices[xyIndex], rowColIndices[xyIndex+1]-1, &depth)) {
			*dataPtr = &data.vec[depth];
			return true;
		}
		*dataPtr = NULL;
		return false;
	}

	inline unsigned int getRowColIndexSize() const { return rowColIndices.size(); }
	inline unsigned int getDepthIndexSize() const { return depthIndices.size(); }

	inline vector<unsigned int>* getRowColIndexPtr() { return &rowColIndices.vec; }
	inline vector<unsigned int>* getDepthIndexPtr() { return &depthIndices.vec; }
	inline vector<T>* getDataPtr() { return &data.vec; }

	inline double getNumberOfBytes() const {
		return TemplateSize::getNumberOfBytes(rowColIndices) + TemplateSize::getNumberOfBytes(depthIndices) + 
			   TemplateSize::getNumberOfBytes(data);
	}

	inline void writeData(FILE *file) const {
		TemplateFileIO::writeBinary(file, rowColIndices);
		TemplateFileIO::writeBinary(file, depthIndices);
		TemplateFileIO::writeBinary(file, data);
	}
	inline void readData(FILE *file) {
		TemplateFileIO::readBinary(file, &rowColIndices);
		TemplateFileIO::readBinary(file, &depthIndices);
		TemplateFileIO::readBinary(file, &data);
	}

	void copyDataFrom(const Grid3D_Regular<T> &reg, const T &nullValue) {
		double numberOfActualEntries = 0;
		for (unsigned int n=0; n<reg.size(); n++) {
			if (reg.get(n) != nullValue) numberOfActualEntries++;
		}

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

			rowColIndices.vec[n+1] = rowColIndices.vec[n] + rowColSize;
		}

		depthIndices.resize(rowColIndices[rowColIndices.size()-1]);
		data.resize(rowColIndices[rowColIndices.size()-1], outsideVal);

		unsigned int count = 0;
		for (unsigned int n=0; n<rowColIndices.size()-1; n++) {
			unsigned int x, y;
			getReversedRowColMapping(n, x, y);

			if (x < reg.numCells.x && y < reg.numCells.y) {
				for (unsigned int z=0; z<reg.numCells.z; z++) {
					if (reg.get(x,y,z) != nullValue) {
						depthIndices.vec[count] = z;
						data.vec[count] = reg.get(x,y,z);
						count++;
					}
				}
			}
		}
	}

	void copyDataFrom(const Grid3D_Regular_Sparse_I1<T> &unfixedRep) {
		unsigned int n = rowColIndices.size()-1;

		rowColIndices.vec[0] = 0;
		for (unsigned int i=0; i<n; i++) {
			unsigned int x, y;
			getReversedRowColMapping(i, x, y);

			unsigned int rowColSize = 0;
			if (x < unfixedRep.numCells.x && y < unfixedRep.numCells.y) {
				rowColSize = unfixedRep.depthIndices[x][y].size();
			}

			rowColIndices.vec[i+1] = rowColIndices.vec[i] + rowColSize;
		}

		depthIndices.resize(rowColIndices[n]);
		data.resize(rowColIndices[n], outsideVal);

		unsigned int count = 0;
		for (unsigned int i=0; i<n; i++) {
			unsigned int x, y;
			getReversedRowColMapping(i, x, y);

			if (x < unfixedRep.numCells.x && y < unfixedRep.numCells.y) {
				for (unsigned int k=0; k<unfixedRep.depthIndices[x][y].size(); k++) {
					depthIndices.vec[count] = unfixedRep.depthIndices[x][y][k];
					data.vec[count] = unfixedRep.data[x][y][k];

					count++;
				}
			}
		}
	}

	inline bool getDataIndex(const unsigned int x, const unsigned int y, const unsigned int z, unsigned int &index) const {
		unsigned int depth, xyIndex = rowColIndexMap(x,y);
		if (depthIndices.binarySearch_sortedLowToHigh(z, rowColIndices[xyIndex], rowColIndices[xyIndex+1]-1, &depth)) {
			index = depth;
			return true;
		}
		else return false;
	}

protected:
	StlVector<unsigned int> rowColIndices;
	StlVector<unsigned int> depthIndices;
	StlVector<T> data;

	virtual inline unsigned int rowColIndexMap(const unsigned int x, const unsigned int y) const { return x*numCells.y + y; }
	virtual inline void getReversedRowColMapping(const unsigned int index, unsigned int &x, unsigned int &y) {
		x = index * numCellsInv.y;
		y = index % numCells.y;
	}

	virtual inline unsigned int getRowColIndexMapLeft(const unsigned int idx) const { return idx-numCells.y; }
	virtual inline unsigned int getRowColIndexMapRight(const unsigned int idx) const { return idx+numCells.y; }
	virtual inline unsigned int getRowColIndexMapDown(const unsigned int idx) const { return idx-1; }
	virtual inline unsigned int getRowColIndexMapUp(const unsigned int idx) const { return idx+1; }
};

template <class T> class Grid3D_Pow2_CSR : public Grid3D_Regular_CSR<T> {
public:
	Grid3D_Pow2_CSR() { }

	Grid3D_Pow2_CSR(Vector3ui numberOfCells, T outsideValue)
		{ Grid3D_Regular_CSR::Grid3D_Regular_CSR(getPowerOf2Indices(numberOfCells), outsideValue); }

	Grid3D_Pow2_CSR(Vector3ui numberOfCells, T outsideValue, Vector3d startPosition, Vector3d cellSizes)
		{ Grid3D_Regular_CSR::Grid3D_Regular_CSR(getPowerOf2Indices(numberOfCells), outsideValue, startPosition, cellSizes); }

	Grid3D_Pow2_CSR(const Grid3D_Regular<T> &reg, const T &nullValue) {
		outsideVal = nullValue;
		rebuildGrid(reg.numCells);
		setGridPosition(reg.startPos, reg.cellSize);
		copyDataFrom(reg, nullValue);
	}

	Grid3D_Pow2_CSR(const Grid3D_Regular_Sparse_I1<T> &unfixedRep) {
		outsideVal = unfixedRep.outsideVal;
		rebuildGrid(unfixedRep.numCells);
		setGridPosition(unfixedRep.startPos, unfixedRep.cellSize);
		copyDataFrom(unfixedRep);
	}

	void rebuildGrid(Vector3ui numberOfCells) { Grid3D_Regular_CSR::rebuildGrid(getPowerOf2Indices(numberOfCells)); }

	inline int getType() { return GRID_POW2_CSR; };

	Vector3ui initialNumCells;

	inline unsigned int rowColIndexMap(const unsigned int x, const unsigned int y) const { return MathUtils::interleaveBits(x,y); }
	inline void getReversedRowColMapping(const unsigned int index, unsigned int &x, unsigned int &y) { MathUtils::unInterleaveBits(index,x,y); }

protected:
	inline Vector3ui getPowerOf2Indices(const Vector3ui originalNumberOfCells) {
		initialNumCells = originalNumberOfCells;
		return Vector3ui(MathUtils::roundUpPowerOf2(originalNumberOfCells.x), MathUtils::roundUpPowerOf2(originalNumberOfCells.y), MathUtils::roundUpPowerOf2(originalNumberOfCells.z));
	}
};

template <> inline bool* Grid3D_Regular_CSR<bool>::getPtr(const unsigned int x, const unsigned int y, const unsigned int z) { return NULL; }
template <> inline bool Grid3D_Regular_CSR<bool>::entryExists(const unsigned int x, const unsigned int y, const unsigned int z, bool **dataPtr) { return entryExists(x,y,z); }



#endif