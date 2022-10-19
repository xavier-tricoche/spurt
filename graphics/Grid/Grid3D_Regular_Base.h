#ifndef GRID_3D_REGULAR_BASE_H
#define GRID_3D_REGULAR_BASE_H

#include "Math/VectorN.h"
#include "Math/DirectionConstants.h"
#include "TemplateHelper.h"

#include "Math/MathUtils.h"
using namespace MathUtils;

#include <vector>
using namespace std;

enum GridTypes { UNKNOWN_GRID_TYPE, GRID_REGULAR, GRID_REGULAR_SPARSE_I1, GRID_REGULAR_SPARSE_I2,
				 GRID_REGULAR_SPARSE_I3, GRID_REGULAR_SPARSE_VEC, GRID_REGULAR_CSR,
				 GRID_REGULAR_SPARSE_NODATA, GRID_POW2, GRID_POW2_CSR, GRID_POW2_SPARSE_NODATA,
				 GRID_REGULAR_BLOCK, GRID_REGULAR_BLOCK_NODATA };

enum Grid2DTypes { UNKNOWN_GRID2D_TYPE, GRID2D_REGULAR, GRID2D_REGULAR_SPARSE_I1, GRID2D_REGULAR_SPARSE_I2,
				   GRID2D_REGULAR_SPARSE_I3, GRID2D_REGULAR_SPARSE_VEC, GRID2D_REGULAR_CSR,
				   GRID2D_REGULAR_SPARSE_NODATA, GRID2D_POW2, GRID2D_POW2_CSR, GRID2D_POW2_SPARSE_NODATA };

template <class T> class Grid3D_Regular_Base {
public:
	Grid3D_Regular_Base() { }
	virtual ~Grid3D_Regular_Base() { }

	Vector3ui numCells;
	Vector3i numCellsMinus1;
	Vector3d startPos, endPos, cellSize, invCellSize, halfCellSize, invHalfCellSize, numCellsInv, doubleCellSize, invDoubleCellSize;

	T outsideVal;

	virtual inline int getType() const = 0;
	virtual inline void deletePointers() = 0;
	virtual inline void clear() = 0;

	void setGridPosition(Vector3d startPosition, Vector3d cellSizes) {
		startPos = startPosition;
		cellSize = cellSizes;
		endPos = getCellMaxPos(numCells.x-1, numCells.y-1, numCells.z-1);
		invCellSize = 1.0 / cellSize;
		halfCellSize = cellSize*0.5;
		invHalfCellSize = 1.0 / halfCellSize;
		doubleCellSize = 2.0 * cellSize;
		invDoubleCellSize = 1.0 / doubleCellSize;
	}

	void setGridPositionUsingBoundingBox(Vector3d startPosition, Vector3d endPosition)
		{ setGridPosition(startPosition, (endPosition-startPosition)/Vector3d(numCells.x, numCells.y, numCells.z)); }

	virtual void rebuildGrid(Vector3ui numberOfCells) {
		numCells = numberOfCells;

		if (numCells.x == 0) numCells.x++;
		if (numCells.y == 0) numCells.y++;
		if (numCells.z == 0) numCells.z++;

		numCellsInv = 1.0/Vector3d(numCells.x, numCells.y, numCells.z);
		numCellsYZ = numCells.y*numCells.z;
		numCellsYZInv = 1.0/(double)numCellsYZ;
		numCellsMinus1 = numCells-1;
	}

	virtual inline void resetAllToValue(const T &defaultValue) { }

	virtual inline bool getDataIndex(const unsigned int x, const unsigned int y, const unsigned int z, unsigned int &index) const {
		index = getCellIndex(x,y,z);
		return true;
	}

	virtual inline unsigned int getCellIndex(const unsigned int x, const unsigned int y, const unsigned int z) const { return x*numCellsYZ + y*numCells.z + z; }
	inline int getCellIndex(const Vector3d &worldSpacePos) const {
		Vector3d index = worldSpaceToCellSpaceTransform(worldSpacePos);
		Vector3i idx(round(index.x), round(index.y), round(index.z));
		return getCellIndex(idx.x, idx.y, idx.z);
	}
	inline int getCellIndex_checkBounds(const Vector3d &worldSpacePos, const bool clampToNearest = false) const {
		Vector3d index = worldSpaceToCellSpaceTransform(worldSpacePos);
		Vector3i idx(round(index.x), round(index.y), round(index.z));

		if (areIndicesInBounds(idx.x, idx.y, idx.z)) return getCellIndex(idx.x, idx.y, idx.z);
		else if (clampToNearest) return getCellIndex(clamp(idx.x, 0, numCellsMinus1.x), clamp(idx.y, 0, numCellsMinus1.y), clamp(idx.z, 0, numCellsMinus1.z));
		else return -1;
	}

	virtual inline unsigned int getCellIndexComponentX(const unsigned int index) const { return (unsigned int)(index*numCellsYZInv); }
	virtual inline unsigned int getCellIndexComponentY(const unsigned int index) const { return ((int)(index*numCellsInv.z))%numCells.y; }
	virtual inline unsigned int getCellIndexComponentZ(const unsigned int index) const { return index%numCells.z; }
	virtual inline Vector3ui getCellIndices(const unsigned int index) const { return Vector3ui(getCellIndexComponentX(index), getCellIndexComponentY(index), getCellIndexComponentZ(index)); }
	inline Vector3i getCellIndices(const Vector3d &worldSpacePos) const {
		Vector3d index = worldSpaceToCellSpaceTransform(worldSpacePos);
		return Vector3i(round(index.x), round(index.y), round(index.z));
	}
	inline Vector3ui getCellIndices_clampToBounds(const Vector3d &worldSpacePos) const {
		Vector3d index = worldSpaceToCellSpaceTransform(worldSpacePos);
		return Vector3ui(clamp(round(index.x), 0, numCellsMinus1.x), clamp(round(index.y), 0, numCellsMinus1.y), clamp(round(index.z), 0, numCellsMinus1.z));
	}

	inline int getNeighborCellIndex(const int index, const int direction) const {
		switch (direction) {
			case FACE_RIGHT: return getRightCellIndex(index);	case FACE_LEFT: return getLeftCellIndex(index);
			case FACE_UP:	 return getUpCellIndex(index);		case FACE_DOWN: return getDownCellIndex(index);
			case FACE_FRONT: return getFrontCellIndex(index);	case FACE_BACK: return getBackCellIndex(index);
			default: printf("Bad direction specified\n");  return -1;
		}
	}
	inline int getRightCellIndex(const int index) const { return index+numCellsYZ; }
	inline int getLeftCellIndex(const int index) const { return index-numCellsYZ; }
	inline int getUpCellIndex(const int index) const { return index+numCells.z; }
	inline int getDownCellIndex(const int index) const { return index-numCells.z; }
	inline int getFrontCellIndex(const int index) const { return index+1; }
	inline int getBackCellIndex(const int index) const { return index-1; }

	inline int getNeighborCellIndex_checkBounds(const int index, const int direction) const {
		switch (direction) {
			case FACE_RIGHT: return getRightCellIndex_checkBounds(index);	case FACE_LEFT: return getLeftCellIndex_checkBounds(index);
			case FACE_UP: return getUpCellIndex_checkBounds(index);			case FACE_DOWN: return getDownCellIndex_checkBounds(index);
			case FACE_FRONT: return getFrontCellIndex_checkBounds(index);	case FACE_BACK: return getBackCellIndex_checkBounds(index);
			default: printf("Bad direction specified\n");  return -1;
		}
	}
	inline int getRightCellIndex_checkBounds(const int index) const { if ((int)getCellIndexComponentX(index) < numCellsMinus1.x) return index+numCellsYZ;  else return -1; }
	inline int getLeftCellIndex_checkBounds(const int index) const { if ((int)getCellIndexComponentX(index) > 0) return index-numCellsYZ;  else return -1; }
	inline int getUpCellIndex_checkBounds(const int index) const { if ((int)getCellIndexComponentY(index) < numCellsMinus1.y) return index+numCells.z;  else return -1; }
	inline int getDownCellIndex_checkBounds(const int index) const { if ((int)getCellIndexComponentY(index)) return index-numCells.z;  else return -1; }
	inline int getFrontCellIndex_checkBounds(const int index) const { if ((int)getCellIndexComponentZ(index) < numCellsMinus1.z) return index+1;  else return -1; }
	inline int getBackCellIndex_checkBounds(const int index) const { if ((int)getCellIndexComponentZ(index) > 0) return index-1;  else return -1; }

	inline bool areIndicesInBounds(const int x, const int y, const int z) const { return x >= 0 && x <= numCellsMinus1.x && y >= 0 && y <= numCellsMinus1.y && z >= 0 && z <= numCellsMinus1.z; }
	inline bool areIndicesInBounds(Vector3i indices) const { return areIndicesInBounds(indices.x, indices.y, indices.z); }

	inline void clampIndices(int &x, int &y, int &z) const {
		x = clamp(x, 0, numCellsMinus1.x);
		y = clamp(y, 0, numCellsMinus1.y);
		z = clamp(z, 0, numCellsMinus1.z);
	}
	inline void clampIndices(Vector3i &indices) const { clampIndices(indices.x, indices.y, indices.z); }
	
	virtual inline void set(const unsigned int x, const unsigned int y, const unsigned int z, const T &val) = 0;
	virtual inline T set_returnOld(const unsigned int x, const unsigned int y, const unsigned int z, const T &val) = 0;

	virtual inline T get(const unsigned int x, const unsigned int y, const unsigned int z) const = 0;
	virtual inline T* getPtr(const unsigned int x, const unsigned int y, const unsigned int z) = 0;

	virtual inline void get2x2x2Block(const unsigned int lowX, const unsigned int lowY, const unsigned int lowZ,
									T &val000, T &val001, T &val010, T &val011, T &val100, T &val101, T &val110, T &val111) const { }

	virtual inline bool entryExists(const unsigned int x, const unsigned int y, const unsigned int z) const { return true; }
	virtual inline bool entryExists(const unsigned int x, const unsigned int y, const unsigned int z, T **dataPtr)
		{ *dataPtr = getPtr(x,y,z);  return true; }
	
	inline T getNearest(const Vector3d &worldSpacePos) const { Vector3i idx = getCellIndices_clampToBounds(worldSpacePos);  return get(idx.x, idx.y, idx.z); }
	
	inline T getLinear(const Vector3d &worldSpacePos) const {
		int i0, i1, j0, j1, k0, k1;
		Vector3d cellSpacePos = worldSpaceToCellSpaceTransform(worldSpacePos);
		getBounds(cellSpacePos, i0, i1, j0, j1, k0, k1);

		return triLinearInterpolation(cellSpacePos, i0, i1, j0, j1, k0, k1);
	}
	
	inline T getCubic(const Vector3d &worldSpacePos) const {
		int i0, i1, j0, j1, k0, k1;
		Vector3d cellSpacePos = worldSpaceToCellSpaceTransform(worldSpacePos);
		getBounds(cellSpacePos, i0, i1, j0, j1, k0, k1);

		return catmullRomInterpolation(cellSpacePos, i0, i1, j0, j1, k0, k1, true, true);
	}

	inline void findIndexBoundsOfValues(Vector3i &minIndex, Vector3i &maxIndex, const T &lowVal, const T &highVal) const {
		minIndex.set(numCells.x+1, numCells.y+1, numCells.z+1);
		maxIndex.set(-1,-1,-1);
		for (unsigned int i=0; i<numCells.x; i++)  for (unsigned int j=0; j<numCells.y; j++)  for (unsigned int k=0; k<numCells.z; k++) {
			T val = get(i,j,k);
			if (val >= lowVal && val <= highVal) { minIndex = minIndex.minVector(i,j,k);  maxIndex = maxIndex.maxVector(i,j,k); }
		}
	}

	inline Vector3d getBoundingBoxMin() const { return startPos; }
	inline Vector3d getBoundingBoxMax() const { return endPos; }
	inline Vector3d getBoundingBoxCenter() const { return (endPos+startPos)*0.5; }
	inline Vector3d getBoundingBoxRadius() const { return (endPos-startPos)*0.5; }
	inline Vector3d getBoundingBoxDimensions() const { return endPos-startPos; }

	inline Vector3d getPositionOfData(const unsigned int x, const unsigned int y, const unsigned int z) const { return getCellCenter(x,y,z); }
	
	inline Vector3d getCellCenter(const unsigned int x, const unsigned int y, const unsigned int z) const { return cellSpaceToWorldSpaceTransform(Vector3d(x,y,z)); }

	inline Vector3d getCellMaxPos(const unsigned int x, const unsigned int y, const unsigned int z) const { return getCellMaxPos(getCellMinPos(x,y,z)); }
	inline Vector3d getCellMaxPos(const Vector3d &cellMinPos) const { return cellMinPos+cellSize; }
	
	inline Vector3d getCellMinPos(const unsigned int x, const unsigned int y, const unsigned int z) const { return (Vector3d(x,y,z) * cellSize) + startPos; }
	inline Vector3d getCellMinPos(const Vector3d &cellMaxPos) const { return cellMaxPos-cellSize; }

	inline void getCellMinMaxPos(const unsigned int x, const unsigned int y, const unsigned int z, Vector3d &minPos, Vector3d &maxPos) const
		{ minPos = getCellMinPos(x,y,z);  maxPos = getCellMaxPos(minPos); }

	inline Vector3d getLeftFaceCenter(const unsigned int x, const unsigned int y, const unsigned int z) const { return cellSpaceToWorldSpaceTransform(Vector3d((double)x-0.5,y,z)); }
	inline Vector3d getRightFaceCenter(const unsigned int x, const unsigned int y, const unsigned int z) const { return cellSpaceToWorldSpaceTransform(Vector3d((double)x+0.5,y,z)); }
	inline Vector3d getDownFaceCenter(const unsigned int x, const unsigned int y, const unsigned int z) const { return cellSpaceToWorldSpaceTransform(Vector3d(x,(double)y-0.5,z)); }
	inline Vector3d getUpFaceCenter(const unsigned int x, const unsigned int y, const unsigned int z) const { return cellSpaceToWorldSpaceTransform(Vector3d(x,(double)y+0.5,z)); }
	inline Vector3d getBackFaceCenter(const unsigned int x, const unsigned int y, const unsigned int z) const { return cellSpaceToWorldSpaceTransform(Vector3d(x,y,(double)z-0.5)); }
	inline Vector3d getFrontFaceCenter(const unsigned int x, const unsigned int y, const unsigned int z) const { return cellSpaceToWorldSpaceTransform(Vector3d(x,y,(double)z+0.5)); }
	inline Vector3d getFaceCenter(const unsigned int x, const unsigned int y, const unsigned int z, const int direction) const {
		switch (direction) {
			case FACE_RIGHT: return getRightFaceCenter(index);	case FACE_LEFT: return getLeftFaceCenter(index);
			case FACE_UP:    return getUpFaceCenter(index);		case FACE_DOWN: return getDownFaceCenter(index);
			case FACE_FRONT: return getFrontFaceCenter(index);	case FACE_BACK: return getBackFaceCenter(index);
			default: printf("Bad direction specified\n");  return Vector3d(0,0,0);
		}
	}

	// cell space is defined as the position of the data
	//   cell space (0,0,0) corresponds to index[0][0][0],  (1,2,4) -> index[1][2][4]
	inline Vector3d worldSpaceToCellSpaceTransform(const Vector3d &worldSpacePosition) const { return ((worldSpacePosition-startPos) * invCellSize) - 0.5; }
	inline Vector3d cellSpaceToWorldSpaceTransform(const Vector3d &cellSpacePos) const { return ((cellSpacePos+0.5) * cellSize) + startPos; }

	inline bool isPointInVolume(const Vector3d &pt) const
		{ return startPos.x <= pt.x && startPos.y <= pt.y && startPos.z <= pt.z && endPos.x >= pt.x && endPos.y >= pt.y && endPos.z >= pt.z; }

	vector<Vector3ui> getSurroundingGridCenterIndices(const Vector3d &worldSpacePos) const {
		vector<Vector3ui> indices;

		Vector3d cellIdx = (worldSpacePos-startPos) * invCellSize;
		cellIdx.roundDown();

		for (int i=0; i<8; i++) {
			Vector3ui idx = Vector3ui(cellIdx.x, cellIdx.y, cellIdx.z);
			if (i >= 4) idx.x++;
			if (i%4 >= 2) idx.y++;
			if (i%2 == 1) idx.z++;

			if (idx.x >= 0 && idx.x < numCells.x && idx.y >= 0 && idx.y < numCells.y && idx.z >= 0 && idx.z < numCells.z)
				indices.push_back(idx);
		}

		return indices;
	}

	// first derivative estimated using finite difference - accuracy O(h^2)
	// x' = [x(i-1) - x(i+1)] / [2*|x|]
	inline Vector3<T> getCentralDifference12(const unsigned int x, const unsigned int y, const unsigned int z) const
		{ return Vector3<T>(getCentralDifference12_X(x,y,z), getCentralDifference12_Y(x,y,z), getCentralDifference12_Z(x,y,z)); }
	inline T getCentralDifference12_X(const unsigned int x, const unsigned int y, const unsigned int z) const
		{ return get((x > 0) ? x-1 : x,y,z)-get((x < numCells.x-1) ? x+1 : x,y,z) * invDoubleCellSize; }
	inline T getCentralDifference12_Y(const unsigned int x, const unsigned int y, const unsigned int z) const
		{ return get(x,(y > 0) ? y-1 : y,z)-get(x,(y < numCells.y-1) ? y+1 : y,z) * invDoubleCellSize; }
	inline T getCentralDifference12_Z(const unsigned int x, const unsigned int y, const unsigned int z) const
		{ return get(x,y,(z > 0) ? z-1 : z)-get(x,y,(z < numCells.z-1) ? z+1 : z) * invDoubleCellSize; }

	// x' = [x(i-1) + x(i+1) - 2*x(i)] / [2*|x|]
	inline Vector3<T> getCentralDifference22(const unsigned int x, const unsigned int y, const unsigned int z) const {
		unsigned int xMinus, yMinus, zMinus, xPlus, yPlus, zPlus;
		getNeighboringIndices_checkBounds(x, y, z, xMinus, yMinus, zMinus, xPlus, yPlus, zPlus);
		return (Vector3<T>(get(xMinus,y,z)+get(xPlus,y,z), get(x,yMinus,z)+get(x,yPlus,z), get(x,y,zMinus)+get(x,y,zPlus)) - 2*get(x,y,z))*invDoubleCellSize;
	}

	// to do: verify that we need to divide by 2*|x| and not 4... if its 4 need to change other stuff below as well
	// x' = [-x(i-2) + 2*x(i-1) - 2*x(i+1) + x(i+2)] / [2*|x|]
	inline Vector3<T> getCentralDifference32(const unsigned int x, const unsigned int y, const unsigned int z) const {
		unsigned int xMinus1, yMinus1, zMinus1, xMinus2, yMinus2, zMinus2, xPlus1, yPlus1, zPlus1, xPlus2, yPlus2, zPlus2;
		getNeighboringIndices_checkBounds(x, y, z, xMinus1, yMinus1, zMinus1, xMinus2, yMinus2, zMinus2, xPlus1, yPlus1, zPlus1, xPlus2, yPlus2, zPlus2);
		return Vector3<T>(2*(get(xMinus1,y,z)-get(xPlus1,y,z)) + get(xPlus2,y,z) - get(xMinus2,y,z),
						  2*(get(x,yMinus1,z)-get(x,yPlus1,z)) + get(x,yPlus2,z) - get(x,yMinus2,z),
						  2*(get(x,y,zMinus1)-get(x,y,zPlus1)) + get(x,y,zPlus2) - get(x,y,zMinus2)) * invDoubleCellSize;
	}

	// x' = [6*x(i) - 4*x(i-1) - 4*(x+1) + x(i-2) + x(i+2)] / [2*|x|]
	inline Vector3<T> getCentralDifference42(const unsigned int x, const unsigned int y, const unsigned int z) const {
		unsigned int xMinus1, yMinus1, zMinus1, xMinus2, yMinus2, zMinus2, xPlus1, yPlus1, zPlus1, xPlus2, yPlus2, zPlus2;
		getNeighboringIndices_checkBounds(x, y, z, xMinus1, yMinus1, zMinus1, xMinus2, yMinus2, zMinus2, xPlus1, yPlus1, zPlus1, xPlus2, yPlus2, zPlus2);
		return (6*get(x,y,z) + Vector3<T>(get(xPlus2,y,z) + get(xMinus2,y,z) - 4*(get(xMinus1,y,z)+get(xPlus1,y,z)),
										  get(x,yPlus2,z) + get(x,yMinus2,z) - 4*(get(x,yMinus1,z)+get(x,yPlus1,z)),
										  get(x,y,zPlus2) + get(x,y,zMinus2) - 4*(get(x,y,zMinus1)+get(x,y,zPlus1)))) * invDoubleCellSize;
	}

	// to do:
	/****************************************************************************
	* Central differences with accuracy O(h^4)  
	****************************************************************************/
	//BZ_DECLARE_DIFF(central14) { return 8.0 * (A.shift(1,dim) - A.shift(-1,dim)) - (A.shift(2,dim) - A.shift(-2,dim)); }
	//BZ_DECLARE_DIFF(central24) { return -30.0 * (*A) + 16.0 * (A.shift(1,dim) + A.shift(-1,dim)) - (A.shift(2,dim) + A.shift(-2,dim)); }
	//BZ_DECLARE_DIFF(central34) { return -13.0 * (A.shift(1,dim) - A.shift(-1,dim)) + 8.0 * (A.shift(2,dim) - A.shift(-2,dim)) - (A.shift(3,dim) - A.shift(-3,dim)); }
	//BZ_DECLARE_DIFF(central44) { return 56.0 * (*A) - 39.0 * (A.shift(1,dim) + A.shift(-1,dim)) + 12.0 * (A.shift(2,dim) + A.shift(-2,dim)) - (A.shift(3,dim) + A.shift(-3,dim)); }
	/****************************************************************************
	 * Backward differences with accuracy O(h)    (forward differences just shift in the positive direction)
	 ****************************************************************************/
	//BZ_DECLARE_DIFF(backward11) { return (*A) - A.shift(-1,dim); }
	//BZ_DECLARE_DIFF(backward21) { return (*A) - 2.0 * A.shift(-1,dim) + A.shift(-2,dim); }
	//BZ_DECLARE_DIFF(backward31) { return (*A) - 3.0 * A.shift(-1,dim) + 3.0 * A.shift(-2,dim) - A.shift(-3,dim); }
	//BZ_DECLARE_DIFF(backward41) { return (*A) - 4.0 * A.shift(-1,dim) + 6.0 * A.shift(-2,dim) - 4.0 * A.shift(-3,dim) + A.shift(-4,dim); }
	/****************************************************************************
	 * Backward differences with accuracy O(h^2)
	 ****************************************************************************/
	//BZ_DECLARE_DIFF(backward12) { return 3.0 * (*A) - 4.0 * A.shift(-1,dim) + A.shift(-2,dim); 
	//BZ_DECLARE_DIFF(backward22) { return 2.0 * (*A) - 5.0 * A.shift(-1,dim) + 4.0 * A.shift(-2,dim) - A.shift(-3,dim); }
	//BZ_DECLARE_DIFF(backward32) { return 5.0 * (*A) - 18.0 * A.shift(-1,dim) + 24.0 * A.shift(-2,dim) - 14.0 * A.shift(-3,dim) + 3.0 * A.shift(-4,dim); }
	//BZ_DECLARE_DIFF(backward42) { return 3.0 * (*A) - 14.0 * A.shift(-1,dim) + 26.0 * A.shift(-2,dim) - 24.0 * A.shift(-3,dim) + 11.0 * A.shift(-4,dim) - 2.0 * A.shift(-5,dim); }

	inline Vector3<T> getGradient(Vector3d worldPos, bool useCubic) const {
		Vector3<T> grad;
		if (useCubic) {
			grad.x = getCubic(worldPos + Vector3d(-cellSize.x, 0, 0)) + getCubic(worldPos + Vector3d(cellSize.x, 0, 0));
			grad.y = getCubic(worldPos + Vector3d(0, -cellSize.y, 0)) + getCubic(worldPos + Vector3d(0, cellSize.y, 0));
			grad.z = getCubic(worldPos + Vector3d(0, 0, -cellSize.z)) + getCubic(worldPos + Vector3d(0, 0, cellSize.z));
		}
		else {
			grad.x = getLinear(worldPos + Vector3d(-cellSize.x, 0, 0)) + getLinear(worldPos + Vector3d(cellSize.x, 0, 0));
			grad.y = getLinear(worldPos + Vector3d(0, -cellSize.y, 0)) + getLinear(worldPos + Vector3d(0, cellSize.y, 0));
			grad.z = getLinear(worldPos + Vector3d(0, 0, -cellSize.z)) + getLinear(worldPos + Vector3d(0, 0, cellSize.z));
		}
		return grad*invDoubleCellSize;
	}
	inline Vector3<T> getNormal(Vector3d worldPos, bool useCubic = false) const { return getGradient(worldPos, useCubic).unit(); }

	inline Vector3<T> getGradient(const unsigned int x, const unsigned int y, const unsigned int z) const {
		unsigned int xMinus, yMinus, zMinus, xPlus, yPlus, zPlus;
		getNeighboringIndices_checkBounds(x, y, z, xMinus, yMinus, zMinus, xPlus, yPlus, zPlus);

		Vector3<T> grad;
		grad.x = get(xMinus,y,z)+get(xPlus,y,z);
		grad.y = get(x,yMinus,z)+get(x,yPlus,z);
		grad.z = get(x,y,zMinus)+get(x,y,zPlus);
		return grad*invDoubleCellSize;
	}
	inline Vector3<T> getNormal(const unsigned int x, const unsigned int y, const unsigned int z) const { return getGradient(x,y,z).unit(); }

	//  laplacian = divergence (gradient (x))
	//  to do:  this is wrong, need to do a grid of the gradients then do the divergence on that grid
	/*inline T getLaplacian(const unsigned int x, const unsigned int y, const unsigned int z) const {
		Vector3<T> grad = getGradient(x,y,z);
		return grad.x + grad.y + grad.z;
	}*/

	T getCurvature(const unsigned int x, const unsigned int y, const unsigned int z) const {
		unsigned int xMinus, yMinus, zMinus, xPlus, yPlus, zPlus;
		getNeighboringIndices_checkBounds(x, y, z, xMinus, yMinus, zMinus, xPlus, yPlus, zPlus);

		T phi_xMinus = get(xMinus, y, z);
		T phi_yMinus = get(x, yMinus, z);
		T phi_zMinus = get(x, y, zMinus);
		T phi_xPlus = get(xPlus, y, z);
		T phi_yPlus = get(x, yPlus, z);
		T phi_zPlus = get(x, y, zPlus);

		// first derivatives estimated using finite difference
		//   x' = [x(i-1) - x(i+1)] / [2*|x|]
		T phi_x = (phi_xMinus - phi_xPlus) * invCellSize.x * 0.5;
		T phi_y = (phi_yMinus - phi_yPlus) * invCellSize.y * 0.5;
		T phi_z = (phi_zMinus - phi_zPlus) * invCellSize.z * 0.5;

		// second derivatives
		//   x'' = [x(i-1) - 2*x(i) + x(i+1)] / [|x|^sq]
		T twoPhi = 2.0*getValue(x, y, z);
		T phi_xx = (phi_xMinus - twoPhi + phi_xPlus) * invCellSize.x * invCellSize.x;
		T phi_yy = (phi_yMinus - twoPhi + phi_yPlus) * invCellSize.y * invCellSize.y;
		T phi_zz = (phi_zMinus - twoPhi + phi_zPlus) * invCellSize.z * invCellSize.z;

		T phi_xy = (get(xMinus, yMinus, z) + get(xPlus, yPlus, z) - get(xMinus, yPlus, z) - 
						 get(xPlus, yMinus, z)) * 0.25 * invCellSize.x * invCellSize.y;
		T phi_xz = (get(xMinus, y, zMinus) + get(xPlus, y, zPlus) - get(xMinus, y, zPlus) - 
						 get(xPlus, y, zMinus)) * 0.25 * invCellSize.x * invCellSize.z;
		T phi_yz = (get(x, yMinus, zMinus) + get(x, yPlus, zPlus) - get(x, yMinus, zPlus) - 
						 get(x, yPlus, zMinus)) * 0.25 * invCellSize.y * invCellSize.z;

		// Eq. 9 from Osher and Fedkiw 2002 "Level Set Methods and Dynamic Implicit Surfaces"
		T phi_x2 = phi_x*phi_x;
		T phi_y2 = phi_y*phi_y;
		T phi_z2 = phi_z*phi_z;

		T gradPhi = sqrt(phi_x2 + phi_y2 + phi_z2);
		T gradPhiCubed = gradPhi*gradPhi*gradPhi;
		if (gradPhiCubed == 0) gradPhiCubed = 1;

		return (phi_x2*phi_yy - 2.0*phi_x*phi_y*phi_xy + phi_y2*phi_xx + phi_x2*phi_zz - 2.0*phi_x*phi_z*phi_xz +
				phi_z2*phi_xx + phi_y2*phi_zz - 2.0*phi_y*phi_z*phi_yz + phi_z2*phi_yy) / gradPhiCubed;
	}

	inline void getNeighboringIndices_checkBounds(const unsigned int x, const unsigned int y, const unsigned int z,
												  unsigned int &xMinus, unsigned int &yMinus, unsigned int &zMinus,
												  unsigned int &xPlus, unsigned int &yPlus, unsigned int &zPlus) const {
		xMinus = (x > 0) ? x-1 : x;				yMinus = (y > 0) ? y-1 : y;				zMinus = (z > 0) ? z-1 : z;
		xPlus = (x < numCells.x-1) ? x+1 : x;	yPlus = (y < numCells.y-1) ? y+1 : y;	zPlus = (z < numCells.z-1) ? z+1 : z;
	}

	inline void getNeighboringIndices_checkBounds(const unsigned int x, const unsigned int y, const unsigned int z,
												  unsigned int &xMinus1, unsigned int &yMinus1, unsigned int &zMinus1,
												  unsigned int &xMinus2, unsigned int &yMinus2, unsigned int &zMinus2,
												  unsigned int &xPlus1, unsigned int &yPlus1, unsigned int &zPlus1,
												  unsigned int &xPlus2, unsigned int &yPlus2, unsigned int &zPlus2) const {
		if (x > 1)	  { xMinus1 = x-1;  xMinus2 = x-2; }
		else if (x > 0) xMinus1 = xMinus2 = x-1;
		else			xMinus1 = xMinus2 = x;

		if (y > 1)	  { yMinus1 = y-1;  yMinus2 = y-2; }
		else if (y > 0) yMinus1 = yMinus2 = y-1;
		else			yMinus1 = yMinus2 = y;

		if (z > 1)	  { zMinus1 = z-1;  zMinus2 = z-2; }
		else if (z > 0) zMinus1 = zMinus2 = z-1;
		else			zMinus1 = zMinus2 = z;

		if (x < numCells.x-2)	 { xPlus1 = x+1;  xPlus2 = x+2; }
		else if (x < numCells.x-1) xPlus1 = xPlus2 = x+1;
		else					   xPlus1 = xPlus2 = x;

		if (y < numCells.y-2)	 { yPlus1 = y+1;  yPlus2 = y+2; }
		else if (y < numCells.y-1) yPlus1 = yPlus2 = y+1;
		else					   yPlus1 = yPlus2 = y;

		if (z < numCells.z-2)	 { zPlus1 = z+1;  zPlus2 = z+2; }
		else if (z < numCells.z-1) zPlus1 = zPlus2 = z+1;
		else					   zPlus1 = zPlus2 = z;
	}

	inline void getSurroundingValues(Vector3d &cellSpacePos, T &val000, T &val001, T &val010, T &val011, T &val100, T &val101, T &val110, T &val111) const {
		int i0, i1, j0, j1, k0, k1;
		getBounds(cellSpacePos, i0, i1, j0, j1, k0, k1);

		if (i0 == i1) {
			if (j0 == j1) {
				if (k0 == k1) val000 = val001 = val010 = val011 = val100 = val101 = val110 = val111 = get(i0,j0,k0);
				else { val000 = val010 = val100 = val110 = get(i0,j0,k0);	val001 = val011 = val101 = val111 = get(i0,j0,k1); }
			}
			else {
				if (k0 == k1) { val000 = val001 = val100 = val101 = get(i0,j0,k0);	val010 = val011 = val110 = val111 = get(i0,j1,k0); }
				else { val000 = val100 = get(i0,j0,k0);	val001 = val101 = get(i0,j0,k1);	val010 = val110 = get(i0,j1,k0);	val011 = val111 = get(i0,j1,k1); }
			}
		}
		else {
			if (j0 == j1) {
				if (k0 == k1) { val000 = val001 = val010 = val011 = get(i0,j0,k0);	val100 = val101 = val110 = val111 = get(i1,j0,k0); }
				else { val000 = val010 = get(i0,j0,k0);	val001 = val011 = get(i0,j0,k1);	val100 = val110 = get(i1,j0,k0);	val101 = val111 = get(i1,j0,k1); }
			}
			else {
				val000 = get(i0,j0,k0);		val001 = get(i0,j0,k1);		val010 = get(i0,j1,k0);		val011 = get(i0,j1,k1);
				val100 = get(i1,j0,k0);		val101 = get(i1,j0,k1);		val110 = get(i1,j1,k0);		val111 = get(i1,j1,k1);
			}
		}
	}

	inline Vector3d getSurroundingIndices(const Vector3d &worldSpacePos, int &i0, int &i1, int &j0, int &j1, int &k0, int &k1) const {
		Vector3d cellSpacePos = worldSpaceToCellSpaceTransform(worldSpacePos);
		getBounds(cellSpacePos, i0, i1, j0, j1, k0, k1);
		return cellSpacePos;
	}

	virtual inline double getNumberOfBytes() const { return 0; }
	
	inline void writeType(FILE *file) const { int type = getType();  fwrite(&type, sizeof(type), 1, file); }

	inline void writeGridProperties(FILE *file) const {
		fwrite(numCells.data, sizeof(numCells.x), 3, file);
		fwrite(startPos.data, sizeof(startPos.x), 3, file);
		fwrite(cellSize.data, sizeof(cellSize.x), 3, file);
	}
	inline void readGridProperties(FILE *file) {
		fread(numCells.data, sizeof(numCells.x), 3, file);
		fread(startPos.data, sizeof(startPos.x), 3, file);
		fread(cellSize.data, sizeof(cellSize.x), 3, file);

		setGridPosition(startPos, cellSize);
		rebuildGrid(numCells);
	}

	inline void writeOutsideValue(FILE *file) const { TemplateFileIO::writeBinary(file, outsideVal); }
	inline void readOutsideValue(FILE *file) { TemplateFileIO::readBinary(file, &outsideVal); }

	virtual inline void writeData(FILE *file) const { printf("Unimplemented Grid3D::writeData specialization\n"); }
	virtual inline void readData(FILE *file) { printf("Unimplemented Grid3D::readData specialization\n"); }

	inline void writeAll(FILE *file) const {
		writeType(file);
		writeGridProperties(file);
		writeOutsideValue(file);
		writeData(file);
	}
	inline void writeAll(const string filename) const {
		FILE *file = fopen(filename.c_str(), "wb");
		if (file == NULL) {
			printf("Could not write grid to %s\n", filename.c_str());
			return;
		}
		writeAll(file);
		fclose(file);
	}

protected:
	unsigned int numCellsYZ;
	double numCellsYZInv;

	inline void getBounds(Vector3d &cellSpacePos, int &i0, int &i1, int &j0, int &j1, int &k0, int &k1) const {
		static double tolerance = 1e-8;

		MathUtils::getBoundingIntegers(cellSpacePos.x, i0, i1, tolerance);
		MathUtils::getBoundingIntegers(cellSpacePos.y, j0, j1, tolerance);
		MathUtils::getBoundingIntegers(cellSpacePos.z, k0, k1, tolerance);

		MathUtils::clamp2(i0, i1, 0, numCellsMinus1.x);
		MathUtils::clamp2(j0, j1, 0, numCellsMinus1.y);
		MathUtils::clamp2(k0, k1, 0, numCellsMinus1.z);

		cellSpacePos.clamp(Vector3d(0,0,0), Vector3d(numCellsMinus1.x, numCellsMinus1.y, numCellsMinus1.z));
	}

	inline T get_checkBounds(const int x, const int y, const int z, const bool clampIndex) const {
		if (x < 0) { if (clampIndex) x = 0;  else return outsideVal; }
		else if (x >= (int)numCells.x) { if (clampIndex) x = numCells.x-1;  else return outsideVal; }

		if (y < 0) { if (clampIndex) y = 0;  else return outsideVal; }
		else if (y >= (int)numCells.y) { if (clampIndex) y = numCells.y-1;  else return outsideVal; }

		if (z < 0) { if (clampIndex) z = 0;  else return outsideVal; }
		else if (z >= (int)numCells.z) { if (clampIndex) z = numCells.z-1;  else return outsideVal; }

		return get(x,y,z);
	}

	inline T triLinearInterpolation(const Vector3d &cellSpacePos, const int i0, const int i1, const int j0, const int j1, const int k0, const int k1) const {
		// short cicuit: potential to be faster and more accurate
		if (i0 == i1) {
			if (j0 == j1) {
				if (k0 == k1) return get(i0,j0,k0);
				else return lerp(k1-cellSpacePos.z, get(i0,j0,k0), get(i0,j0,k1));
			}
			else {
				if (k0 == k1) return lerp(j1-cellSpacePos.y, get(i0,j0,k0), get(i0,j1,k0));
				else return bilinearLerp(get(i0,j0,k0), get(i0,j0,k1), get(i0,j1,k0), get(i0,j1,k1), j1-cellSpacePos.y, k1-cellSpacePos.z);
			}
		}
		else {
			if (j0 == j1) {
				if (k0 == k1) return lerp(i1-cellSpacePos.x, get(i0,j0,k0), get(i1,j0,k0));
				else return bilinearLerp(get(i0,j0,k0), get(i0,j0,k1), get(i1,j0,k0), get(i1,j0,k1), i1-cellSpacePos.x, k1-cellSpacePos.z);
			}
			else return trilinearLerp(get(i0, j0, k0), get(i0, j0, k1), get(i0, j1, k0), get(i0, j1, k1),
									  get(i1, j0, k0), get(i1, j0, k1), get(i1, j1, k0), get(i1, j1, k1),
									  i1-cellSpacePos.x, j1-cellSpacePos.y, k1-cellSpacePos.z);
		}
	}

	inline T catmullRomInterpolation(const Vector3d &cellSpacePos, const int i0, const int i1, const int j0, const int j1, const int k0, const int k1,
							  bool preserveMonotonicity, bool clampToLocalBounds) const {
		// short cicuit: potential to be faster and more accurate
		if (i0 == i1) {
			if (j0 == j1) {
				if (k0 == k1) return get(i0, j0, k0);
				else return catmullRomInterpolate(get(i0,j0,k0 == 0 ? 0 : k0-1), get(i0,j0,k0),
												  get(i0,j0,k1), get(i0,j0,k1 == numCellsMinus1.z ? k1 : k1+1),
												  cellSpacePos.z-k0, preserveMonotonicity, clampToLocalBounds);
			}
			else {
				if (k0 == k1) return catmullRomInterpolate(get(i0,j0 == 0 ? 0 : j0-1,k0), get(i0,j0,k0),
														   get(i0,j1,k0), get(i0,j1 == numCellsMinus1.y ? j1 : j1+1,k0),
														   cellSpacePos.y-j0, preserveMonotonicity, clampToLocalBounds);
				else {
					int jm = j0==0 ? 0 : j0-1, jp = j1==numCellsMinus1.y ? j1 : j1+1;
					int km = k0==0 ? 0 : k0-1, kp = k1==numCellsMinus1.z ? k1 : k1+1;
					return bicubicCatmullRomInterpolate(get(i0,jm,km), get(i0,jm,k0), get(i0,jm,k1), get(i0,jm,kp),
														get(i0,j0,km), get(i0,j0,k0), get(i0,j0,k1), get(i0,j0,kp),
														get(i0,j1,km), get(i0,j1,k0), get(i0,j1,k1), get(i0,j1,kp),
														get(i0,jp,km), get(i0,jp,k0), get(i0,jp,k1), get(i0,jp,kp),
														cellSpacePos.y-j0, cellSpacePos.z-k0, preserveMonotonicity, clampToLocalBounds);
				}
			}
		}
		else {
			if (j0 == j1) {
				if (k0 == k1) return catmullRomInterpolate(get(i0 == 0 ? 0 : i0-1,j0,k0), get(i0,j0,k0), 
														   get(i1,j0,k0), get(i1 == numCellsMinus1.x ? i1 : i1+1,j0,k0),
														   cellSpacePos.x-i0, preserveMonotonicity, clampToLocalBounds);
				else {
					int im = i0==0 ? 0 : i0-1, ip = i1==numCellsMinus1.x ? i1 : i1+1;
					int km = k0==0 ? 0 : k0-1, kp = k1==numCellsMinus1.z ? k1 : k1+1;
					return bicubicCatmullRomInterpolate(get(im,j0,km), get(im,j0,k0), get(im,j0,k1), get(im,j0,kp),
														get(i0,j0,km), get(i0,j0,k0), get(i0,j0,k1), get(i0,j0,kp),
														get(i1,j0,km), get(i1,j0,k0), get(i1,j0,k1), get(i1,j0,kp),
														get(ip,j0,km), get(ip,j0,k0), get(ip,j0,k1), get(ip,j0,kp),
														cellSpacePos.x-i0, cellSpacePos.z-k0, preserveMonotonicity, clampToLocalBounds);
				}
			}
			else {
				if (k0 == k1) {
					int im = i0==0 ? 0 : i0-1, ip = i1==numCellsMinus1.x ? i1 : i1+1;
					int jm = j0==0 ? 0 : j0-1, jp = j1==numCellsMinus1.y ? j1 : j1+1;
					return bicubicCatmullRomInterpolate(get(im,jm,k0), get(im,j0,k0), get(im,j1,k0), get(im,jp,k0),
														get(i0,jm,k0), get(i0,j0,k0), get(i0,j1,k0), get(i0,jp,k0),
														get(i1,jm,k0), get(i1,j0,k0), get(i1,j1,k0), get(i1,jp,k0),
														get(ip,jm,k0), get(ip,j0,k0), get(ip,j1,k0), get(ip,jp,k0),
														cellSpacePos.x-i0, cellSpacePos.y-j0, preserveMonotonicity, clampToLocalBounds);
				}
				else {
					int im = i0==0 ? 0 : i0-1, ip = i1==numCellsMinus1.x ? i1 : i1+1;
					int jm = j0==0 ? 0 : j0-1, jp = j1==numCellsMinus1.y ? j1 : j1+1;
					int km = k0==0 ? 0 : k0-1, kp = k1==numCellsMinus1.z ? k1 : k1+1;
					return tricubicCatmullRomInterpolate(
							get(im,jm,km), get(im,jm,k0), get(im,jm,k1), get(im,jm,kp), get(im,j0,km), get(im,j0,k0), get(im,j0,k1), get(im,j0,kp),
							get(im,j1,km), get(im,j1,k0), get(im,j1,k1), get(im,j1,kp), get(im,jp,km), get(im,jp,k0), get(im,jp,k1), get(im,jp,kp),
							get(i0,jm,km), get(i0,jm,k0), get(i0,jm,k1), get(i0,jm,kp), get(i0,j0,km), get(i0,j0,k0), get(i0,j0,k1), get(i0,j0,kp),
							get(i0,j1,km), get(i0,j1,k0), get(i0,j1,k1), get(i0,j1,kp), get(i0,jp,km), get(i0,jp,k0), get(i0,jp,k1), get(i0,jp,kp),
							get(i1,jm,km), get(i1,jm,k0), get(i1,jm,k1), get(i1,jm,kp), get(i1,j0,km), get(i1,j0,k0), get(i1,j0,k1), get(i1,j0,kp),
							get(i1,j1,km), get(i1,j1,k0), get(i1,j1,k1), get(i1,j1,kp), get(i1,jp,km), get(i1,jp,k0), get(i1,jp,k1), get(i1,jp,kp),
							get(ip,jm,km), get(ip,jm,k0), get(ip,jm,k1), get(ip,jm,kp), get(ip,j0,km), get(ip,j0,k0), get(ip,j0,k1), get(ip,j0,kp),
							get(ip,j1,km), get(ip,j1,k0), get(ip,j1,k1), get(ip,j1,kp), get(ip,jp,km), get(ip,jp,k0), get(ip,jp,k1), get(ip,jp,kp),
							cellSpacePos.x-i0, cellSpacePos.y-j0, cellSpacePos.z-k0, preserveMonotonicity, clampToLocalBounds);
				}
			}
		}
		
		
	}
};

#endif