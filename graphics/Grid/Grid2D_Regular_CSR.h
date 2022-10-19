#ifndef GRID_2D_REGULAR_CSR_H
#define GRID_2D_REGULAR_CSR_H

#include "Grid2D_Regular_Base.h"
#include "Grid2D_Regular_Sparse_I1.h"
//#include "Grid2D_Regular.h"

#include "StlVectorUtils.h"

template <class T> class Grid2D_Regular_CSR : public Grid2D_Regular_Base<T> {
public:
	Grid2D_Regular_CSR() { }

	Grid2D_Regular_CSR(Vector2ui numberOfCells, const T &outsideValue) {
		rebuildGrid(numberOfCells, outsideValue);
		setGridPosition(Vector2d(0,0), Vector2d(1,1));
	}

	Grid2D_Regular_CSR(Vector2ui numberOfCells, const T &outsideValue, Vector2d startPosition, Vector2d cellSizes) {
		rebuildGrid(numberOfCells, outsideValue);
		setGridPosition(startPosition, cellSizes);
	}

	//Grid2D_Regular_CSR(const Grid3D_Regular<T> &reg, const T &nullValue) {
	//	rebuildGrid(reg.numCells, reg.outsideVal);
	//	convertFrom(reg, nullValue);
	//}

	Grid2D_Regular_CSR(const Grid2D_Regular_Sparse_I1<T> &unfixedRep) {
		rebuildGrid(unfixedRep.numCells, unfixedRep.outsideVal);
		convertFrom(unfixedRep);
	}

	~Grid2D_Regular_CSR() { clear(); }

	inline void clear() {
		if (data.size() == 0) return;

		rowIndices.clear();
		columnIndices.clear();
		data.deleteElements();
		data.clear();
	}

	inline void deletePointers() { data.deleteElements(); }

	inline int getType() { return GRID2D_REGULAR_CSR; }

	virtual void rebuildGrid(Vector2ui numberOfCells, T outsideValue) {
		Grid2D_Regular_Base::rebuildGrid(numberOfCells, outsideValue);

		rowIndices.resize(numberOfCells.x+1);
		for (unsigned int i=0; i<rowIndices.size(); i++) rowIndices.vec[i] = i; // to do: I don't like this, should be zero
																			//   but typically not used since this type is
																			//   usually built from another matrix

		columnIndices.resize(numberOfCells.x);
		for (unsigned int i=0; i<columnIndices.size(); i++) columnIndices.vec[i] = i;

		data.resize(numberOfCells.x, outsideValue);
	}

	inline void set(const unsigned int x, const unsigned int y, const T &val) {
		unsigned int column;
		if (columnIndices.binarySearch_sortedLowToHigh(y, rowIndices[x], rowIndices[x+1]-1, &column)) {
			data.vec[column] = val;
		}
		else {
			columnIndices.insert(column, y);
			data.insert(column, val);

			for (unsigned int i=x+1; i<rowIndices.size(); i++) rowIndices.vec[i]++;
		}
	}

	inline T set_returnOld(const unsigned int x, const unsigned int y, const T &val) {
		unsigned int column;
		if (columnIndices.binarySearch_sortedLowToHigh(y, rowIndices[x], rowIndices[x+1]-1, &column)) {
			T old = data[column];  data.vec[column] = val;  return old;
		}
		else {
			columnIndices.insert(column, y);
			data.insert(column, val);

			for (unsigned int i=x+1; i<rowIndices.size(); i++) rowIndices.vec[i]++;

			return outsideVal;
		}
	}

	inline T get(const unsigned int x, const unsigned int y) const {
		unsigned int column;
		if (columnIndices.binarySearch_sortedLowToHigh(y, rowIndices[x], rowIndices[x+1]-1, &column)) {
			return data.vec[column];
		}
		return outsideVal;
	}

	inline T* getPtr(const unsigned int x, const unsigned int y) {
		unsigned int column;
		if (columnIndices.binarySearch_sortedLowToHigh(y, rowIndices[x], rowIndices[x+1]-1, &column)) {
			return &data.vec[column];
		}
		return NULL;
	}

	inline void get2x2Block(const unsigned int lowX, const unsigned int lowY, T &val00, T &val01, T &val10, T &val11) const {
		unsigned int column;
		if (columnIndices.binarySearch_sortedLowToHigh(lowY, rowIndices[lowX], rowIndices[lowX+1]-1, &column)) {
			val00 = data.vec[column];
		}
		else val00 = outsideVal;

		// since we are fetching a 2x2 block, we just need to check the next depth index
		//   this will work even if the above value does not exist since the binary search returns the index of the
		//   first value that is less than what we searched for
		if (column+1 < rowIndices[lowX+1] && columnIndices[column+1] == lowY+1) val01 = data.vec[column+1];
		else val01 = outsideVal;


		// chances are the position of what we just found will be close to the other rows, so use it to speed up the search
		//   worst case, this will basically increase the numbers of iterations of the binary search by 1
		unsigned int offset = column-rowIndices[lowX];

		unsigned int x = lowX+1;
		unsigned int dim = rowIndices[x+1]-rowIndices[x];

		bool found;
		if (offset < dim) {
			column = rowIndices[x]+offset;
			if (lowY == columnIndices[column]) found = true;
			else if (lowY < columnIndices[column]) found = columnIndices.binarySearch_sortedLowToHigh(lowY, rowIndices[x], column, &column);
			else/*(lowY > columnIndices[column])*/ found = columnIndices.binarySearch_sortedLowToHigh(lowY, column, rowIndices[x+1]-1, &column);
		}
		else found = columnIndices.binarySearch_sortedLowToHigh(lowY, rowIndices[x], rowIndices[x+1]-1, &column);
		
		if (found) val10 = data.vec[column];
		else val10 = outsideVal;

		if (column+1 < rowIndices[x+1] && columnIndices[column+1] == lowY+1) val11 = data.vec[column+1];
		else val11 = outsideVal;
	}

	inline bool entryExists(const unsigned int x, const unsigned int y) const
		{ return columnIndices.binarySearch_sortedLowToHigh(y, rowIndices[x], rowIndices[x+1]-1); }

	inline unsigned int getRowColIndexSize() const { return rowIndices.size(); }
	inline unsigned int getDepthIndexSize() const { return columnIndices.size(); }

	inline vector<unsigned int>* getRowIndexPtr() { return &rowIndices.vec; }
	inline vector<unsigned int>* getColumnIndexPtr() { return &columnIndices.vec; }
	inline vector<T>* getDataPtr() { return &data.vec; }

	inline double getNumberOfBytes() const;

protected:
	StlVector<unsigned int> rowIndices;
	StlVector<unsigned int> columnIndices;
	StlVector<T> data;

	virtual inline unsigned int getRowColIndexMapLeft(const unsigned int idx) const { return idx-numCells.y; }
	virtual inline unsigned int getRowColIndexMapRight(const unsigned int idx) const { return idx+numCells.y; }
	virtual inline unsigned int getRowColIndexMapDown(const unsigned int idx) const { return idx-1; }
	virtual inline unsigned int getRowColIndexMapUp(const unsigned int idx) const { return idx+1; }

	/*void convertFrom(const Grid2D_Regular<T> &reg, const T &nullValue) {
		setGridPosition(reg.startPos, reg.cellSize);

		double numberOfActualEntries = 0;
		for (unsigned int n=0; n<reg.size(); n++) {
			if (reg.get(n) != nullValue) numberOfActualEntries++;
		}

		rowColIndices[0] = 0;

		for (unsigned int n=0; n<rowColIndices.size()-1; n++) {
			unsigned int x, y;
			getReversedRowColMapping(n, x, y);

			unsigned int rowColSize = 0;

			if (x < reg.numCells.x && y < reg.numCells.y) {
				for (unsigned int z=0; z<reg.numCells.z; z++) {
					if (reg.get(x,y,z) != nullValue) rowColSize++;
				}
			}

			rowColIndices[n+1] = rowColIndices[n] + rowColSize;
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
						depthIndices[count] = z;
						data.vec[count] = reg.get(x,y,z);
						count++;
					}
				}
			}
		}
	}*/

	void convertFrom(const Grid2D_Regular_Sparse_I1<T> &unfixedRep) {
		setGridPosition(unfixedRep.startPos, unfixedRep.cellSize);

		unsigned int n = rowIndices.size()-1;

		rowIndices[0] = 0;
		for (unsigned int i=0; i<n; i++) {
			rowIndices[i+1] = rowIndices[i] + unfixedRep.columnIndices[i].size();
		}

		columnIndices.resize(rowIndices[n]);
		data.resize(rowIndices[n], outsideVal);

		unsigned int count = 0;
		for (unsigned int i=0; i<n; i++) {
			for (unsigned int j=0; j<unfixedRep.columnIndices[i].size(); j++) {
				columnIndices[count] = unfixedRep.columnIndices[i][j];
				data.vec[count] = unfixedRep.data[i][j];
				count++;
			}
		}
	}
};

template <class T> inline double Grid2D_Regular_CSR<T>::getNumberOfBytes() const {
	return sizeof(rowIndices) + sizeof(columnIndices) + sizeof(data) + rowIndices.size()*sizeof(unsigned int) + columnIndices.size()*sizeof(unsigned int) + data.size()*sizeof(T);
}

// I don't know how to get the template pattern matching to work nicely with this stuff (without making a specialized template class)
//   so instead, when a new type is added, just add a specialization here
template <> inline double Grid2D_Regular_CSR< vector<int> >::getNumberOfBytes() const {
	double size = sizeof(rowIndices) + sizeof(columnIndices) + sizeof(data) + rowIndices.size()*sizeof(unsigned int) + columnIndices.size()*sizeof(unsigned int);
	for (unsigned int i=0; i<data.size(); i++) { size += (double)data[i].size()*sizeof(int); }  return size;
}
template <> inline double Grid2D_Regular_CSR< vector<unsigned int> >::getNumberOfBytes() const {
	double size = sizeof(rowIndices) + sizeof(columnIndices) + sizeof(data) + rowIndices.size()*sizeof(unsigned int) + columnIndices.size()*sizeof(unsigned int);
	for (unsigned int i=0; i<data.size(); i++) { size += (double)data[i].size()*sizeof(unsigned int); }  return size;
}


#endif