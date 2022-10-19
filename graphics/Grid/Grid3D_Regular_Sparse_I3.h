#ifndef GRID_3D_REGULAR_SPARSE_I3_H
#define GRID_3D_REGULAR_SPARSE_I3_H

#include "Grid3D_Regular_Base.h"

#include "StlVectorUtils.h"

template <class T> class Grid3D_Regular_Sparse_I3 : public Grid3D_Regular_Base<T> {
public:
	Grid3D_Regular_Sparse_I3() { }

	Grid3D_Regular_Sparse_I3(Vector3ui numberOfCells, T outsideValue) {
		rebuildGrid(numberOfCells, outsideValue);
		setGridPosition(Vector3d(0,0,0), Vector3d(1,1,1));
	}

	Grid3D_Regular_Sparse_I3(Vector3ui numberOfCells, T outsideValue, Vector3d startPosition, Vector3d cellSizes) {
		rebuildGrid(numberOfCells, outsideValue);
		setGridPosition(startPosition, cellSizes);
	}

	~Grid3D_Regular_Sparse_I3() { clear(false); }

	void clear(bool freePointers) {
		if (data.size() == 0) return;

		/*if (freePointers) {
			for (unsigned int i=0; i<data.size(); i++)
			for (unsigned int j=0; j<data[i].size(); j++)
			for (unsigned int k=0; k<data[i][j].size(); k++)
				if (data[i][j][k] != NULL) { delete data[i][j][k];  data[i][j][k] = NULL; }
		}*/

		for (unsigned int i=0; i<data.size(); i++) {
			for (unsigned int j=0; j<data[i].size(); j++) {
				data[i][j].clear();
				depthIndices[i][j].clear();
			}
			
			data[i].clear();
			colIndices[i].clear();
			depthIndices[i].clear();
		}
		data.clear();
		rowIndices.clear();
		colIndices.clear();
		depthIndices.clear();
	}

	inline int getType() { return GRID_REGULAR_SPARSE_I3; }

	inline void set(const unsigned int x, const unsigned int y, const unsigned int z, const T &val) {
		unsigned int row;
		if (rowIndices.find(x, &row)) {
			unsigned int col;
			if (colIndices[row].find(y, &col)) {
				unsigned int depth;
				if (depthIndices[row][col].find(z, &depth)) data[row][col][depth] = val; 
				else {
					depthIndices[row][col].insert(depth, z);
					data[row][col].insert(depth, val);
				}
			}
			else {
				colIndices[row].insert(col, y);

				depthIndices[row].insert(col, SortedStlVector<unsigned int>());
				depthIndices[row][col].add(z);

				data[row].insert(col, StlVector<T>());
				data[row][col].add(val);
			}
		}
		else {
			rowIndices.insert(row, x);

			colIndices.insert(row, SortedStlVector<unsigned int>());
			colIndices[row].add(y);

			depthIndices.insert(row, StlVector< SortedStlVector<unsigned int> >());
			depthIndices[row].add(SortedStlVector<unsigned int>());
			depthIndices[row][depthIndices[row].size()-1].add(z);

			data.insert(row, StlVector< StlVector<T> >());
			data[row].add(StlVector<T>());
			data[row][data[row].size()-1].add(val);
		}
	}

	inline T get(const unsigned int x, const unsigned int y, const unsigned int z) const {
		unsigned int row;
		if (rowIndices.find(x, &row)) {
			unsigned int col;
			if (colIndices[row].find(y, &col)) {
				unsigned int depth;
				if (depthIndices[row][col].find(z, &depth)) return data[row][col][depth]; 
			}
		}
		return outsideVal;

		// old iterative way of doing it, keeping this for future reference if needed
		/*for (unsigned int row=0; row<rowIndices.size(); row++) {
			if (rowIndices[row] == x) {
				for (unsigned int col=0; col<colIndices[row].size(); col++) {
					if (colIndices[row][col] == y) {
						for (unsigned int depth=0; depth<depthIndices[row][col].size(); depth++) {
							if (depthIndices[row][col][depth] == z) return data[row][col][depth];
							else if (depthIndices[row][col][depth] > z) return outsideVal;
						}
						return outsideVal;
					}
					else if (colIndices[row][col] > y) return outsideVal;
				}
				return outsideVal;
			}
			else if (rowIndices[row] > x) return outsideVal;
		}
		return outsideVal;*/
	}

	inline T* getPtr(const unsigned int x, const unsigned int y, const unsigned int z) {
		unsigned int row;
		if (rowIndices.find(x, &row)) {
			unsigned int col;
			if (colIndices[row].find(y, &col)) {
				unsigned int depth;
				if (depthIndices[row][col].find(z, &depth)) return &data[row][col][depth];
			}
		}
		return NULL;
	}

	inline bool entryExists(const unsigned int x, const unsigned int y, const unsigned int z) const {
		unsigned int row;
		if (rowIndices.find(x, &row)) {
			unsigned int col;
			if (colIndices[row].find(y, &col)) {
				unsigned int depth;
				if (depthIndices[row][col].find(z, &depth)) return true;
			}
		}
		return false;
	}

protected:
	SortedStlVector<unsigned int> rowIndices;
	StlVector< SortedStlVector<unsigned int> >colIndices;
	StlVector< StlVector< SortedStlVector<unsigned int> > >depthIndices;
	StlVector< StlVector< StlVector<T> > > data;
};

#endif