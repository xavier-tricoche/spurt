#ifndef GRID_2D_REGULAR_BASE_H
#define GRID_2D_REGULAR_BASE_H

#include "Math/VectorN.h"

#include "Math/DirectionConstants.h"

#include "Math/MathUtils.h"
using namespace MathUtils;

#include <vector>
using namespace std;

template <class T> class Grid2D_Regular_Base {
public:
	Grid2D_Regular_Base() { }

	Vector2ui numCells;
	Vector2i numCellsMinus1;
	Vector2d startPos, endPos, cellSize, invCellSize, halfCellSize, invHalfCellSize, numCellsInv, doubleCellSize, invDoubleCellSize;

	T outsideVal;

	virtual inline int getType() const = 0;
	virtual inline void deletePointers() = 0;
	virtual inline void clear() = 0;

	void setGridPosition(Vector2d startPosition, Vector2d cellSizes) {
		startPos = startPosition;
		cellSize = cellSizes;
		endPos = getCellMaxPos(numCells.x-1, numCells.y-1);
		invCellSize = 1.0 / cellSize;
		halfCellSize = cellSize*0.5;
		invHalfCellSize = 1.0 / halfCellSize;
		doubleCellSize = 2.0 * cellSize;
		invDoubleCellSize = 1.0 / doubleCellSize;
	}

	void setGridPositionUsingBoundingBox(Vector2d startPosition, Vector2d endPosition)
		{ setGridPosition(startPosition, (endPosition-startPosition)/Vector2d(numCells.x, numCells.y)); }

	virtual void rebuildGrid(Vector2ui numberOfCells, T outsideValue) {
		numCells = numberOfCells;

		if (numCells.x == 0) numCells.x++;
		if (numCells.y == 0) numCells.y++;

		outsideVal = outsideValue;

		numCellsInv = 1.0/Vector2d(numCells.x, numCells.y);
		numCellsMinus1 = numCells-1;
	}

	virtual inline void resetAllToValue(const T &defaultValue) { }

	virtual inline unsigned int getCellIndex(const unsigned int x, const unsigned int y) const { return x*numCells.y + y; }
	inline int getCellIndex(const Vector2d &worldSpacePos) const {
		Vector2d index = worldSpaceToCellSpaceTransform(worldSpacePos);
		Vector2i idx(round(index.x), round(index.y));
		return getCellIndex(idx.x, idx.y);
	}
	inline int getCellIndex_checkBounds(const Vector2d &worldSpacePos, const bool clampToNearest = false) const {
		Vector2d index = worldSpaceToCellSpaceTransform(worldSpacePos);
		Vector2i idx(round(index.x), round(index.y));

		if (areIndicesInBounds(idx.x, idx.y)) return getCellIndex(idx.x, idx.y);
		else if (clampToNearest) return getCellIndex(clamp(idx.x, 0, numCellsMinus1.x), clamp(idx.y, 0, numCellsMinus1.y));
		else return -1;
	}

	virtual inline unsigned int getCellIndexComponentX(const unsigned int index) const { return (unsigned int)(index*numCellsInv.y); }
	virtual inline unsigned int getCellIndexComponentY(const unsigned int index) const { return index%numCells.y; }
	virtual inline Vector2ui getCellIndices(const unsigned int index) const { return Vector2ui(getCellIndexComponentX(index), getCellIndexComponentY(index)); }
	inline Vector2i getCellIndices(const Vector2d &worldSpacePos) const {
		Vector2d index = worldSpaceToCellSpaceTransform(worldSpacePos);
		return Vector2i(round(index.x), round(index.y));
	}
	inline Vector2i getCellIndices_clampToBounds(const Vector2d &worldSpacePos) const {
		Vector2d index = worldSpaceToCellSpaceTransform(worldSpacePos);
		return Vector2i(clamp(round(index.x), 0, numCellsMinus1.x), clamp(round(index.y), 0, numCellsMinus1.y));
	}

	inline int getNeighborCellIndex(const int index, const int direction) const {
		switch (direction) {
			case GRID_RIGHT: return getRightCellIndex(index);	case GRID_LEFT: return getLeftCellIndex(index);
			case GRID_UP:	 return getUpCellIndex(index);		case GRID_DOWN: return getDownCellIndex(index);
			default: printf("Bad direction specified\n");  return -1;
		}
	}
	inline int getRightCellIndex(const int index) const { return index+numCells.y; }
	inline int getLeftCellIndex(const int index) const { return index-numCells.y; }
	inline int getUpCellIndex(const int index) const { return index+1; }
	inline int getDownCellIndex(const int index) const { return index-1; }

	inline int getNeighborCellIndex_checkBounds(const int index, const int direction) const {
		switch (direction) {
			case GRID_RIGHT: return getRightCellIndex_checkBounds(index);	case GRID_LEFT: return getLeftCellIndex_checkBounds(index);
			case GRID_UP: return getUpCellIndex_checkBounds(index);		case GRID_DOWN: return getDownCellIndex_checkBounds(index);
			default: printf("Bad direction specified\n");  return -1;
		}
	}
	inline int getRightCellIndex_checkBounds(const int index) const { if ((int)getCellIndexComponentX(index) < numCellsMinus1.x) return index+numCells.y;  else return -1; }
	inline int getLeftCellIndex_checkBounds(const int index) const { if ((int)getCellIndexComponentX(index) > 0) return index-numCells.y;  else return -1; }
	inline int getUpCellIndex_checkBounds(const int index) const { if ((int)getCellIndexComponentY(index) < numCellsMinus1.y) return index+1;  else return -1; }
	inline int getDownCellIndex_checkBounds(const int index) const { if ((int)getCellIndexComponentY(index)) return index-1;  else return -1; }

	inline bool areIndicesInBounds(const int x, const int y) const { return x >= 0 && x <= numCellsMinus1.x && y >= 0 && y <= numCellsMinus1.y; }
	inline bool areIndicesInBounds(Vector2i indices) const { return areIndicesInBounds(indices.x, indices.y); }

	inline void clampIndices(int &x, int &y) const {
		x = clamp(x, 0, numCellsMinus1.x);
		y = clamp(y, 0, numCellsMinus1.y);
	}
	inline void clampIndices(Vector2i &indices) const { clampIndices(indices.x, indices.y); }
	
	virtual inline void set(const unsigned int x, const unsigned int y, const T &val) = 0;
	virtual inline T set_returnOld(const unsigned int x, const unsigned int y, const T &val) = 0;

	virtual inline T get(const unsigned int x, const unsigned int y) const = 0;
	virtual inline T* getPtr(const unsigned int x, const unsigned int y) = 0;

	virtual inline void get2x2Block(const unsigned int lowX, const unsigned int lowY, T &val00, T &val01, T &val10, T &val11) const { }

	virtual inline bool entryExists(const unsigned int x, const unsigned int y) const { return true; }
	
	inline T getNearest(const Vector2d &worldSpacePos) const { Vector2i idx = getCellIndices_clampToBounds(worldSpacePos);  return get(idx.x, idx.y); }
	
	inline T getLinear(const Vector2d &worldSpacePos) const {
		int i0, i1, j0, j1;
		Vector2d cellSpacePos = worldSpaceToCellSpaceTransform(worldSpacePos);
		getBounds(cellSpacePos, i0, i1, j0, j1);

		return biLinearInterpolation(cellSpacePos, i0, i1, j0, j1);
	}
	
	inline T getCubic(const Vector2d &worldSpacePos) const {
		int i0, i1, j0, j1;
		Vector2d cellSpacePos = worldSpaceToCellSpaceTransform(worldSpacePos);
		getBounds(cellSpacePos, i0, i1, j0, j1);

		return catmullRomInterpolation(cellSpacePos, i0, i1, j0, j1, true, true);
	}

	inline void findIndexBoundsOfValues(Vector2i &minIndex, Vector2i &maxIndex, const T &lowVal, const T &highVal) const {
		minIndex.set(numCells.x+1, numCells.y+1);
		maxIndex.set(-1,-1);
		for (unsigned int i=0; i<numCells.x; i++)  for (unsigned int j=0; j<numCells.y; j++) {
			T val = get(i,j);
			if (val >= lowVal && val <= highVal) { minIndex = minIndex.minVector(i,j);  maxIndex = maxIndex.maxVector(i,j); }
		}
	}

	inline Vector2d getBoundingBoxMin() const { return startPos; }
	inline Vector2d getBoundingBoxMax() const { return endPos; }
	inline Vector2d getBoundingBoxCenter() const { return (endPos+startPos)*0.5; }
	inline Vector2d getBoundingBoxRadius() const { return (endPos-startPos)*0.5; }
	inline Vector2d getBoundingBoxDimensions() const { return endPos-startPos; }

	inline Vector2d getPositionOfData(const unsigned int x, const unsigned int y) const { return getCellCenter(x,y); }
	
	inline Vector2d getCellCenter(const unsigned int x, const unsigned int y) const { return cellSpaceToWorldSpaceTransform(Vector2d(x,y)); }

	inline Vector2d getCellMaxPos(const unsigned int x, const unsigned int y) const { return getCellMaxPos(getCellMinPos(x,y)); }
	inline Vector2d getCellMaxPos(const Vector2d &cellMinPos) const { return cellMinPos+cellSize; }
	
	inline Vector2d getCellMinPos(const unsigned int x, const unsigned int y) const { return (Vector2d(x,y) * cellSize) + startPos; }
	inline Vector2d getCellMinPos(const Vector2d &cellMaxPos) const { return cellMaxPos-cellSize; }

	inline Vector2d getLeftFaceCenter(const unsigned int x, const unsigned int y, const unsigned int z) const { return cellSpaceToWorldSpaceTransform(Vector2d((double)x-0.5,y)); }
	inline Vector2d getRightFaceCenter(const unsigned int x, const unsigned int y, const unsigned int z) const { return cellSpaceToWorldSpaceTransform(Vector2d((double)x+0.5,y)); }
	inline Vector2d getDownFaceCenter(const unsigned int x, const unsigned int y, const unsigned int z) const { return cellSpaceToWorldSpaceTransform(Vector2d(x,(double)y-0.5)); }
	inline Vector2d getUpFaceCenter(const unsigned int x, const unsigned int y, const unsigned int z) const { return cellSpaceToWorldSpaceTransform(Vector2d(x,(double)y+0.5)); }
	inline Vector2d getFaceCenter(const unsigned int x, const unsigned int y, const int direction) const {
		switch (direction) {
			case GRID_RIGHT: return getRightFaceCenter(index);	case GRID_LEFT: return getLeftFaceCenter(index);
			case GRID_UP:    return getUpFaceCenter(index);		case GRID_DOWN: return getDownFaceCenter(index);
			default: printf("Bad direction specified\n");  return Vector2d(0,0,0);
		}
	}

	// cell space is defined as the position of the data
	//   cell space (0,0) corresponds to index[0][0],  (1,2) -> index[1][2]
	inline Vector2d worldSpaceToCellSpaceTransform(const Vector2d &worldSpacePosition) const { return ((worldSpacePosition-startPos) * invCellSize) - 0.5; }
	inline Vector2d cellSpaceToWorldSpaceTransform(const Vector2d &cellSpacePos) const { return ((cellSpacePos+0.5) * cellSize) + startPos; }

	inline bool isPointInVolume(const Vector2d &pt) const
		{ return startPos.x <= pt.x && startPos.y <= pt.y && endPos.x >= pt.x && endPos.y >= pt.y; }

	vector<Vector2ui> getSurroundingGridCenterIndices(const Vector2d &worldSpacePos) const {
		vector<Vector2ui> indices;

		Vector2d cellIdx = (worldSpacePos-startPos) * invCellSize;
		cellIdx.roundDown();

		for (int i=0; i<4; i++) {
			Vector2ui idx = Vector2ui(cellIdx.x, cellIdx.y);
			if (i >= 2) idx.x++;
			if (i%2 >= 1) idx.y++;

			if (idx.x >= 0 && idx.x < numCells.x && idx.y >= 0 && idx.y < numCells.y)
				indices.push_back(idx);
		}

		return indices;
	}

	inline Vector2<T> getGradient(Vector2d worldPos) const {
		Vector2<T> grad;
		grad.x = getCubic(worldPos + Vector2d(-cellSize.x, 0)) + getCubic(worldPos + Vector2d(cellSize.x, 0));
		grad.y = getCubic(worldPos + Vector2d(0, -cellSize.y)) + getCubic(worldPos + Vector2d(0, cellSize.y));
		return grad*invDoubleCellSize;
	}
	inline Vector2<T> getNormal(Vector2d worldPos) const { return getGradient(worldPos).unit(); }

	inline Vector2<T> getGradient(const unsigned int x, const unsigned int y) const {
		unsigned int xMinus, yMinus, xPlus, yPlus;
		getNeighboringIndices_checkBounds(x, y, z, xMinus, yMinus, xPlus, yPlus);

		Vector2<T> grad;
		grad.x = get(xMinus,y)+get(xPlus,y);
		grad.y = get(x,yMinus)+get(x,yPlus);
		return grad*invDoubleCellSize;
	}
	inline Vector2<T> getNormal(const unsigned int x, const unsigned int y) const { return getGradient(x,y).unit(); }

	inline void getNeighboringIndices_checkBounds(const unsigned int x, const unsigned int y,
												  unsigned int &xMinus, unsigned int &yMinus, unsigned int &xPlus, unsigned int &yPlus) const {
		xMinus = (x > 0) ? x-1 : x;				yMinus = (y > 0) ? y-1 : y;
		xPlus = (x < numCells.x-1) ? x+1 : x;	yPlus = (y < numCells.y-1) ? y+1 : y;
	}

	inline void getNeighboringIndices_checkBounds(const unsigned int x, const unsigned int y,
												  unsigned int &xMinus1, unsigned int &yMinus1, unsigned int &xMinus2, unsigned int &yMinus2,
												  unsigned int &xPlus1, unsigned int &yPlus1, unsigned int &xPlus2, unsigned int &yPlus2) const {
		if (x > 1)	  { xMinus1 = x-1;  xMinus2 = x-2; }
		else if (x > 0) xMinus1 = xMinus2 = x-1;
		else			xMinus1 = xMinus2 = x;

		if (y > 1)	  { yMinus1 = y-1;  yMinus2 = y-2; }
		else if (y > 0) yMinus1 = yMinus2 = y-1;
		else			yMinus1 = yMinus2 = y;

		if (x < numCells.x-2)	 { xPlus1 = x+1;  xPlus2 = x+2; }
		else if (x < numCells.x-1) xPlus1 = xPlus2 = x+1;
		else					   xPlus1 = xPlus2 = x;

		if (y < numCells.y-2)	 { yPlus1 = y+1;  yPlus2 = y+2; }
		else if (y < numCells.y-1) yPlus1 = yPlus2 = y+1;
		else					   yPlus1 = yPlus2 = y;
	}

	inline void getSurroundingValues(Vector2d &cellSpacePos, T &val00, T &val01, T &val10, T &val11) const {
		int i0, i1, j0, j1;
		getBounds(cellSpacePos, i0, i1, j0, j1);

		if (i0 == i1) {
			if (j0 == j1) val00 = val01 = val10 = val11 = get(i0,j0);
			else { val00 = val10 = get(i0,j0);  val01 = val11 = get(i0,j1); }
		}
		else {
			if (j0 == j1) { val00 = val01 = get(i0,j0);  val10 = val11 = get(i1,j0); }
			else { val00 = get(i0,j0);  val01 = get(i0,j0);  val10 = get(i0,j1);  val11 = get(i0,j1); }
		}
	}

	inline Vector2d getSurroundingIndices(const Vector2d &worldSpacePos, int &i0, int &i1, int &j0, int &j1) const {
		Vector2d cellSpacePos = worldSpaceToCellSpaceTransform(worldSpacePos);
		getBounds(cellSpacePos, i0, i1, j0, j1);
		return cellSpacePos;
	}


	// ----- derivatives below mostly taken from TEEM library -----

	// first derivative estimated using finite difference - accuracy O(h^2)
	// x' = [x(i-1) - x(i+1)] / [2*|x|]
	inline Vector3<T> getCentralDifference12(const unsigned int x, const unsigned int y) const
		{ return Vector2<T>(getCentralDifference12_X(x,y), getCentralDifference12_Y(x,y)); }
	inline T getCentralDifference12_X(const unsigned int x, const unsigned int y) const
		{ return get((x > 0) ? x-1 : x,y)-get((x < numCells.x-1) ? x+1 : x,y) * invDoubleCellSize; }
	inline T getCentralDifference12_Y(const unsigned int x, const unsigned int y) const
		{ return get(x,(y > 0) ? y-1 : y)-get(x,(y < numCells.y-1) ? y+1 : y) * invDoubleCellSize; }

	inline Vector3<T> getCentralDifference12(const Vector2d pos) const
		{ return Vector2<T>(getCentralDifference12_X(pos), getCentralDifference12_Y(pos))*invDoubleCellSize; }
	inline T getCentralDifference12_X(const Vector2d pos) const {
		double minus = pos.x-cellSize.x, plus = pos.x+cellSize.x;

		if (minus < 0) minus = 0;
		else if (plus > endPos.x) plus = endPos.x;
		
		return getLinear(Vector2d(minus, pos.y)) - getLinear(Vector2d(plus, pos.y));
	}
	inline T getCentralDifference12_Y(const Vector2d pos) const {
		double minus = pos.y-cellSize.y, plus = pos.y+cellSize.y;

		if (minus < 0) minus = 0;
		else if (plus > endPos.y) plus = endPos.y;
		
		return getLinear(Vector2d(pos.x, minus)) - getLinear(Vector2d(pos.y, plus));
	}

	// x' = [x(i-1) + x(i+1) - 2*x(i)] / [2*|x|]
	inline Vector3<T> getCentralDifference22(const unsigned int x, const unsigned int y) const {
		unsigned int xMinus, yMinus, xPlus, yPlus;
		getNeighboringIndices_checkBounds(x, y, z, xMinus, yMinus, xPlus, yPlus);
		return (Vector2<T>(get(xMinus,y)+get(xPlus,y), get(x,yMinus)+get(x,yPlus)) - 2*get(x,y))*invDoubleCellSize;
	}

	// to do: verify that we need to divide by 2*|x| and not 4... if its 4 need to change other stuff below as well
	// x' = [-x(i-2) + 2*x(i-1) - 2*x(i+1) + x(i+2)] / [2*|x|]
	inline Vector3<T> getCentralDifference32(const unsigned int x, const unsigned int y) const {
		unsigned int xMinus1, yMinus1, xMinus2, yMinus2, xPlus1, yPlus1, xPlus2, yPlus2;
		getNeighboringIndices_checkBounds(x, y, xMinus1, yMinus1, xMinus2, yMinus2, xPlus1, yPlus1, xPlus2, yPlus2);
		return Vector2<T>(2*(get(xMinus1,y)-get(xPlus1,y)) + get(xPlus2,y) - get(xMinus2,y),
						  2*(get(x,yMinus1)-get(x,yPlus1)) + get(x,yPlus2) - get(x,yMinus2)) * invDoubleCellSize;
	}

	// x' = [6*x(i) - 4*x(i-1) - 4*(x+1) + x(i-2) + x(i+2)] / [2*|x|]
	inline Vector3<T> getCentralDifference42(const unsigned int x, const unsigned int y) const {
		unsigned int xMinus1, yMinus1, xMinus2, yMinus2, xPlus1, yPlus1, xPlus2, yPlus2;
		getNeighboringIndices_checkBounds(x, y, xMinus1, yMinus1, xMinus2, yMinus2, xPlus1, yPlus1, xPlus2, yPlus2);
		return (6*get(x,y,z) + Vector3<T>(get(xPlus2,y) + get(xMinus2,y) - 4*(get(xMinus1,y)+get(xPlus1,y)),
										  get(x,yPlus2) + get(x,yMinus2) - 4*(get(x,yMinus1)+get(x,yPlus1)))) * invDoubleCellSize;
	}


	virtual inline double getNumberOfBytes() const { return 0; }

protected:
	inline void getBounds(Vector2d &cellSpacePos, int &i0, int &i1, int &j0, int &j1) const {
		static double tolerance = 1e-8;

		MathUtils::getBoundingIntegers(cellSpacePos.x, i0, i1, tolerance);
		MathUtils::getBoundingIntegers(cellSpacePos.y, j0, j1, tolerance);

		MathUtils::clamp2(i0, i1, 0, numCellsMinus1.x);
		MathUtils::clamp2(j0, j1, 0, numCellsMinus1.y);

		cellSpacePos.clamp(Vector2d(0,0), Vector2d(numCellsMinus1.x, numCellsMinus1.y));
	}

	inline T get_checkBounds(const int x, const int y, const bool clampIndex) const {
		if (x < 0) { if (clampIndex) x = 0;  else return outsideVal; }
		else if (x >= (int)numCells.x) { if (clampIndex) x = numCells.x-1;  else return outsideVal; }

		if (y < 0) { if (clampIndex) y = 0;  else return outsideVal; }
		else if (y >= (int)numCells.y) { if (clampIndex) y = numCells.y-1;  else return outsideVal; }

		return get(x,y);
	}

	inline T biLinearInterpolation(const Vector2d &cellSpacePos, const int i0, const int i1, const int j0, const int j1) const {
		// short cicuit: potential to be faster and more accurate
		if (i0 == i1) {
			if (j0 == j1) return get(i0,j0);
			else return lerp(j1-cellSpacePos.y, get(i0,j0), get(i0,j1));
		}
		else {
			if (j0 == j1) return lerp(i1-cellSpacePos.x, get(i0,j0), get(i1,j0));
			else return bilinearLerp(get(i0, j0), get(i0, j1), get(i1, j0), get(i1, j1), i1-cellSpacePos.x, j1-cellSpacePos.y);
		}
	}

	inline T catmullRomInterpolation(const Vector2d &cellSpacePos, const int i0, const int i1, const int j0, const int j1,
							  bool preserveMonotonicity, bool clampToLocalBounds) const {
		// short cicuit: potential to be faster and more accurate
		if (i0 == i1) {
			if (j0 == j1) return get(i0, j0);
			else return catmullRomInterpolate(get(i0,j0 == 0 ? 0 : j0-1), get(i0,j0),
											  get(i0,j1), get(i0,j1 == numCellsMinus1.y ? j1 : j1+1),
											  cellSpacePos.y-j0, preserveMonotonicity, clampToLocalBounds);
		}
		else {
			if (j0 == j1) return catmullRomInterpolate(get(i0 == 0 ? 0 : i0-1,j0), get(i0,j0), 
													   get(i1,j0), get(i1 == numCellsMinus1.x ? i1 : i1+1,j0),
													   cellSpacePos.x-i0, preserveMonotonicity, clampToLocalBounds);
			else {
				int im = i0==0 ? 0 : i0-1, ip = i1==numCellsMinus1.x ? i1 : i1+1;
				int jm = j0==0 ? 0 : j0-1, jp = j1==numCellsMinus1.y ? j1 : j1+1;
				return bicubicCatmullRomInterpolate(get(im,jm), get(im,j0), get(im,j1), get(im,jp),
													get(i0,jm), get(i0,j0), get(i0,j1), get(i0,jp),
													get(i1,jm), get(i1,j0), get(i1,j1), get(i1,jp),
													get(ip,jm), get(ip,j0), get(ip,j1), get(ip,jp),
													cellSpacePos.x-i0, cellSpacePos.y-j0, preserveMonotonicity, clampToLocalBounds);
			}
		}
	}
};

#endif