#ifndef GRID_3D_STAGGERED_H
#define GRID_3D_STAGGERED_H

#include "Grid3D_Regular.h"

template <class T> class Grid3D_Staggered {
public:
	Grid3D_Staggered() { gridX = gridY = gridZ = NULL; }

	Grid3D_Staggered(Vector3ui numberOfCells, T outsideValue, Vector3d startPosition = Vector3d(0,0,0), Vector3d cellSizes = Vector3d(1,1,1)) {
		gridX = gridY = gridZ = NULL;

		Vector3ui numCells = numberOfCells;

		if (numCells.x == 0) numCells.x++;
		if (numCells.y == 0) numCells.y++;
		if (numCells.z == 0) numCells.z++;

		gridX = new Grid3D_Regular<T>(numberOfCells+Vector3ui(1,0,0), outsideValue);
		gridY = new Grid3D_Regular<T>(numberOfCells+Vector3ui(0,1,0), outsideValue);
		gridZ = new Grid3D_Regular<T>(numberOfCells+Vector3ui(0,0,1), outsideValue);

		setGridPosition(startPosition, cellSizes);
	}

	Grid3D_Staggered(const Grid3D_Staggered<T> &grid, bool copyGridsData = false) {
		gridX = new Grid3D_Regular<T>(*grid.gridX, copyGridsData);
		gridY = new Grid3D_Regular<T>(*grid.gridY, copyGridsData);
		gridZ = new Grid3D_Regular<T>(*grid.gridZ, copyGridsData);
	}

	Grid3D_Regular<T> *gridX, *gridY, *gridZ;

	~Grid3D_Staggered() {
		if (gridX != NULL) { delete gridX;  gridX = NULL; }
		if (gridY != NULL) { delete gridY;  gridY = NULL; }
		if (gridZ != NULL) { delete gridZ;  gridZ = NULL; }
	}

	void clear(bool freePointers) {
		if (gridX != NULL) {
			gridX->clear(freePointers);
			delete gridX;  gridX = NULL;
		}
		if (gridY != NULL) {
			gridY->clear(freePointers);
			delete gridY;  gridY = NULL;
		}
		if (gridZ != NULL) {
			gridZ->clear(freePointers);
			delete gridZ;  gridZ = NULL;
		}
		totalCells = 0;
	}

	void setGridPosition(Vector3d startPosition, Vector3d cellSizes) {
		gridX->setGridPosition(startPosition-Vector3d(cellSizes.x*0.5,0,0), cellSizes);
		gridY->setGridPosition(startPosition-Vector3d(0,cellSizes.y*0.5,0), cellSizes);
		gridZ->setGridPosition(startPosition-Vector3d(0,0,cellSizes.z*0.5), cellSizes);
	}

	void setGridPositionUsingBoundingBox(Vector3d startPosition, Vector3d endPosition)
		{ setGridPosition(startPosition, (endPosition-startPosition)/Vector3d(numCells.x, numCells.y, numCells.z)); }

	inline void resetAllToValue(T defaultValue) {
		gridX->resetAllToValue(defaultValue);
		gridY->resetAllToValue(defaultValue);
		gridZ->resetAllToValue(defaultValue);
	}

	inline void copyData(const Grid3D_Staggered<T> &grid) {
		gridX->copyData(*grid.gridX);
		gridY->copyData(*grid.gridY);
		gridZ->copyData(*grid.gridZ);
	}

	inline bool areIndicesInBoundsX(int x, int y, int z) const { return gridX->areIndicesInBounds(x,y,z); }
	inline bool areIndicesInBoundsX(Vector3i indices) const { return areIndicesInBoundsX(indices.x, indices.y, indices.z); }
	inline bool isIndexInBoundsX(int index) const { return gridX->isIndexInBounds(index); }

	inline bool areIndicesInBoundsY(int x, int y, int z) const { return gridY->areIndicesInBounds(x,y,z); }
	inline bool areIndicesInBoundsY(Vector3i indices) const { return areIndicesInBoundsY(indices.x, indices.y, indices.z); }
	inline bool isIndexInBoundsY(int index) const { return gridY->isIndexInBounds(index); }

	inline bool areIndicesInBoundsZ(int x, int y, int z) const { return gridZ->areIndicesInBounds(x,y,z); }
	inline bool areIndicesInBoundsZ(Vector3i indices) const { return areIndicesInBoundsZ(indices.x, indices.y, indices.z); }
	inline bool isIndexInBoundsZ(int index) const { return gridZ->isIndexInBounds(index); }

	inline void clampIndicesX(int &x, int &y, int &z) const { gridX->clampIndices(x,y,z); }
	inline void clampIndicesX(Vector3i &indices) const { clampIndicesX(indices.x, indices.y, indices.z); }
	inline void clampIndexX(int &index) const { gridX->clampIndex(index); }

	inline void clampIndicesY(int &x, int &y, int &z) const { gridY->clampIndices(x,y,z); }
	inline void clampIndicesY(Vector3i &indices) const { clampIndicesY(indices.x, indices.y, indices.z); }
	inline void clampIndexY(int &index) const { gridY->clampIndex(index); }

	inline void clampIndicesZ(int &x, int &y, int &z) { gridZ->clampIndices(x,y,z); }
	inline void clampIndicesZ(Vector3i &indices) { clampIndicesZ(indices.x, indices.y, indices.z); }
	inline void clampIndexZ(int &index) { gridZ->clampIndex(index); }

	inline unsigned int getCellIndexX(unsigned int x, unsigned int y, unsigned int z) const { return gridX->getCellIndex(x,y,z); }
	inline unsigned int getCellIndexX(Vector3i indices) const { return getCellIndexX(indices.x, indices.y, indices.z); }
	inline int getCellIndexX(Vector3d worldSpacePos, bool clampToNearest = false) const { return gridX->getCellIndex(worldSpacePos, clampToNearest); }
	inline Vector3ui getCellIndicesX(unsigned int index) const { return gridX->getCellIndices(index); }
	inline Vector3i getCellIndicesX(Vector3d worldSpacePos, bool clampToNearest = false) const { gridX->getCellIndices(worldSpacePos, clampToNearest); }

	inline unsigned int getCellIndexY(unsigned int x, unsigned int y, unsigned int z) const { return gridY->getCellIndex(x,y,z); }
	inline unsigned int getCellIndexY(Vector3i indices) const { return getCellIndexX(indices.x, indices.y, indices.z); }
	inline int getCellIndexY(Vector3d worldSpacePos, bool clampToNearest = false) const { return gridY->getCellIndex(worldSpacePos, clampToNearest); }
	inline Vector3ui getCellIndicesY(unsigned int index) const { return gridY->getCellIndices(index); }
	inline Vector3i getCellIndicesY(Vector3d worldSpacePos, bool clampToNearest = false) const { gridY->getCellIndices(worldSpacePos, clampToNearest); }

	inline unsigned int getCellIndexZ(unsigned int x, unsigned int y, unsigned int z) const { return gridZ->getCellIndex(x,y,z); }
	inline unsigned int getCellIndexZ(Vector3i indices) const { return getCellIndexX(indices.x, indices.y, indices.z); }
	inline int getCellIndexZ(Vector3d worldSpacePos, bool clampToNearest = false) const { return gridZ->getCellIndex(worldSpacePos, clampToNearest); }
	inline Vector3ui getCellIndicesZ(unsigned int index) const { return gridZ->getCellIndices(index); }
	inline Vector3i getCellIndicesZ(Vector3d worldSpacePos, bool clampToNearest = false) const { gridZ->getCellIndices(worldSpacePos, clampToNearest); }

	inline void setX(const int index, const T &val) { gridX->set(index, val); }
	inline void setX(const int x, const int y, const int z, const T &val) { gridX->set(x,y,z,val); }

	inline void setY(int index, const T &val) { gridY->set(index, val); }
	inline void setY(const int x, const int y, const int z, const T &val) { gridY->set(x,y,z,val); }

	inline void setZ(int index, const T &val) { gridZ->set(index, val); }
	inline void setZ(const int x, const int y, const int z, const T &val) { gridZ->set(x,y,z,val); }

	inline T getX(int index) const { return gridX->get(index); }
	inline T getX(const int x, const int y, const int z) const { return gridX->get(x,y,z); }
	inline T getXLinear(Vector3d worldSpacePos) const { return gridX->getLinear(worldSpacePos); }
	inline T getXCubic(Vector3d worldSpacePos) const { return gridX->getCubic(worldSpacePos); }

	inline T getY(int index) const { return gridY->get(index); }
	inline T getY(const int x, const int y, const int z) const { return gridY->get(x,y,z); }
	inline T getYLinear(Vector3d worldSpacePos) const { return gridY->getLinear(worldSpacePos); }
	inline T getYCubic(Vector3d worldSpacePos) const { return gridY->getCubic(worldSpacePos); }

	inline T getZ(int index) const { return gridZ->get(index); }
	inline T getZ(const int x, const int y, const int z) const { return gridZ->get(x,y,z); }
	inline T getZLinear(Vector3d worldSpacePos) const { return gridZ->getLinear(worldSpacePos); }
	inline T getZCubic(Vector3d worldSpacePos) const { return gridZ->getCubic(worldSpacePos); }

	inline T* getPtrX(int index) { return gridX->getPtr(index); }
	inline T* getPtrX(const int x, const int y, const int z) { return gridX->getPtr(x,y,z); }

	inline T* getPtrY(int index) { return gridY->getPtr(index); }
	inline T* getPtrY(const int x, const int y, const int z) { return gridY->getPtr(x,y,z); }

	inline T* getPtrZ(int index) { return gridZ->getPtr(index); }
	inline T* getPtrZ(const int x, const int y, const int z) { return gridZ->getPtr(x,y,z); }

	void advectLinear(double timestep, const VectorField &vectorField) const {
		gridX->advectLinear(timestep, vectorField);
		gridY->advectLinear(timestep, vectorField);
		gridZ->advectLinear(timestep, vectorField);
	}

	void advectCubic(double timestep, const VectorField &vectorField) const {
		gridX->advectCubic(timestep, vectorField);
		gridY->advectCubic(timestep, vectorField);
		gridZ->advectCubic(timestep, vectorField);
	}

	void extrapolate(const Grid3D_Regular<double> &signedDistances, const double maxSignedDistance, Grid3D_Staggered<bool> *isKnown, const T &unknownDefaultValue) {
		gridX->extrapolate_resampleSignedDistance(signedDistances, maxSignedDistance, isKnown->gridX, unknownDefaultValue);
		gridY->extrapolate_resampleSignedDistance(signedDistances, maxSignedDistance, isKnown->gridY, unknownDefaultValue);
		gridZ->extrapolate_resampleSignedDistance(signedDistances, maxSignedDistance, isKnown->gridZ, unknownDefaultValue);
	}
	void extrapolate(const Grid3D_Staggered<double> &signedDistances, const double maxSignedDistance, Grid3D_Staggered<bool> *isKnown, const T &unknownDefaultValue) {
		gridX->extrapolate(*signedDistances.gridX, maxSignedDistance, isKnown->gridX, unknownDefaultValue);
		gridY->extrapolate(*signedDistances.gridY, maxSignedDistance, isKnown->gridY, unknownDefaultValue);
		gridZ->extrapolate(*signedDistances.gridZ, maxSignedDistance, isKnown->gridZ, unknownDefaultValue);
	}

	Vector3<T> getMaxValues() const { return Vector3<T>(gridX->getMaxValue(), gridY->getMaxValue(), gridZ->getMaxValue()); }
	Vector3<T> getMinValues() const { return Vector3<T>(gridX->getMinValue(), gridY->getMinValue(), gridZ->getMinValue()); }

};

#endif