#ifndef GRID_3D_REGULAR_H
#define GRID_3D_REGULAR_H

#include "Grid3D_Regular_Base.h"
#include "Grid3D_Pow2_Base.h"
#include "StlVectorUtils.h"
#include <algorithm>

#pragma warning(disable:4996) // D_SCL_SECURE_NO_WARNINGS for std::copy

#include "CUDA/CudaWrapper.h"
#ifdef COMPILE_WITH_CUDA
	extern "C" void loadGridToCuda(float *volumeData, cudaExtent volumeSize);
	extern "C" void cudaFreeGrid();
#endif

template <class T> class Grid3D_Regular : public Grid3D_Regular_Base<T> {
public:
	Grid3D_Regular() { totalCells = 0;  loadedToCuda = false; }

	Grid3D_Regular(Vector3ui numberOfCells, T outsideValue) {
		totalCells = 0;
		loadedToCuda = false;
		outsideVal = outsideValue;
		rebuildGrid(numberOfCells);
		setGridPosition(Vector3d(0,0,0), Vector3d(1,1,1));
	}

	Grid3D_Regular(Vector3ui numberOfCells, T outsideValue, Vector3d startPosition, Vector3d cellSizes) {
		totalCells = 0;
		loadedToCuda = false;
		outsideVal = outsideValue;
		rebuildGrid(numberOfCells);
		setGridPosition(startPosition, cellSizes);
	}

	template <class S> Grid3D_Regular(const Grid3D_Regular<S> &grid, T outsideVal) {
		totalCells = 0;
		loadedToCuda = false;
		outsideVal = outsideValue;
		rebuildGrid(grid.numCells);
		setGridPosition(grid.startPos, grid.cellSize);
	}

	Grid3D_Regular(const Grid3D_Regular<T> &grid, bool copyGridData = false) {
		totalCells = 0;
		loadedToCuda = false;
		outsideVal = grid.outsideVal;
		rebuildGrid(grid.numCells);
		setGridPosition(grid.startPos, grid.cellSize);
		if (copyGridData) copyData(grid);
	}

	unsigned int totalCells;

	~Grid3D_Regular() { clear(); }

	inline void clear() {
		data.clear();
		totalCells = 0;
	}

	inline void deletePointers() { data.deleteElements(); }

	inline int getType() const { return GRID_REGULAR; }

	virtual void rebuildGrid(Vector3ui numberOfCells) {
		int prevSize = totalCells;

		Grid3D_Regular_Base::rebuildGrid(numberOfCells);

		totalCells = numCells.x*numCells.y*numCells.z;
		if (totalCells == 0) totalCells++;

		if (totalCells != prevSize) {
			data.clear();
			data.resize(totalCells);
		}
	}

	inline void resetAllToValue(const T &defaultValue) { fill(data.vec.begin(), data.vec.end(), defaultValue); }
	inline void copyData(const Grid3D_Regular<T> &grid) {
		if (totalCells != grid.totalCells) printf("Error Grid3D_Regular::copyData: array sizes need to match\n");
		else copy(grid.data.vec.begin(), grid.data.vec.end(), data.vec.begin());
	}

	inline bool isIndexInBounds(const int index) const { return index >= 0 && index < (int)totalCells; }

	inline void clampIndex(int &index) const { clamp(index, 0, (int)totalCells-1); }

	inline void set(const unsigned int index, const T &val) { data.vec[index] = val; }
	inline void set(const unsigned int x, const unsigned int y, const unsigned int z, const T &val) { data.vec[getCellIndex(x,y,z)] = val; }

	inline T set_returnOld(const unsigned int index, const T &val) { T old = data[index];  data.vec[index] = val;  return old; }
	inline T set_returnOld(const unsigned int x, const unsigned int y, const unsigned int z, const T &val) { return set_returnOld(getCellIndex(x,y,z), val); }

	inline T get(const unsigned int index) const { return data.vec[index]; }
	inline T get(const unsigned int x, const unsigned int y, const unsigned int z) const { return data.vec[getCellIndex(x,y,z)]; }
	
	inline T* getPtr(const unsigned int index) { return &data.vec[index]; }
	inline T* getPtr(const unsigned int x, const unsigned int y, const unsigned int z) { return &data.vec[getCellIndex(x,y,z)]; }

	// y,z are fixed, x varies
	inline void get2x2BlockX(const unsigned int lowX, const unsigned int y, const unsigned int z, T &val0, T &val1) const {
		unsigned int idx = getCellIndex(lowX, y, z);
		val0 = get(idx);		val1 = get(getRightCellIndex(idx));
	}

	inline void get2x2BlockY(const unsigned int x, const unsigned int lowY, const unsigned int z, T &val0, T &val1) const {
		unsigned int idx = getCellIndex(x, lowY, z);
		val0 = get(idx);		val1 = get(getUpCellIndex(idx));
	}

	inline void get2x2BlockZ(const unsigned int x, const unsigned int y, const unsigned int lowZ, T &val0, T &val1) const {
		unsigned int idx = getCellIndex(x, y, lowZ);
		val0 = get(idx);		val1 = get(getFrontCellIndex(idx));
	}

	// x is fixed, y and z vary
	inline void get2x2BlockYZ(const unsigned int x, const unsigned int lowY, const unsigned int lowZ,
							  T &val00, T &val01, T &val10, T &val11) const {
		unsigned int idx = getCellIndex(x, lowY, lowZ);
		val00 = get(idx);		val01 = get(getFrontCellIndex(idx));

		idx = getUpCellIndex(idx);
		val10 = get(idx);		val11 = get(getFrontCellIndex(idx));
	}

	inline void get2x2BlockXZ(const unsigned int lowX, const unsigned int y, const unsigned int lowZ,
							  T &val00, T &val01, T &val10, T &val11) const {
		unsigned int idx = getCellIndex(lowX, y, lowZ);
		val00 = get(idx);		val01 = get(getFrontCellIndex(idx));

		idx = getRightCellIndex(idx);
		val10 = get(idx);		val11 = get(getFrontCellIndex(idx));
	}

	inline void get2x2BlockXY(const unsigned int lowX, const unsigned int lowY, const unsigned int z,
							  T &val00, T &val01, T &val10, T &val11) const {
		unsigned int idx = getCellIndex(lowX, lowY, z);
		val00 = get(idx);		val01 = get(getUpCellIndex(idx));

		idx = getRightCellIndex(idx);
		val10 = get(idx);		val11 = get(getUpCellIndex(idx));
	}

	inline void get2x2x2Block(const unsigned int lowX, const unsigned int lowY, const unsigned int lowZ,
							  T &val000, T &val001, T &val010, T &val011, T &val100, T &val101, T &val110, T &val111) const {
		unsigned int idx = getCellIndex(lowX, lowY, lowZ);
		val000 = get(idx);
		val001 = get(getFrontCellIndex(idx));

		idx = getUpCellIndex(idx);
		val010 = get(idx);
		val011 = get(getFrontCellIndex(idx));

		idx = getRightCellIndex(idx);
		val110 = get(idx);
		val111 = get(getFrontCellIndex(idx));

		idx = getDownCellIndex(idx);
		val100 = get(idx);
		val101 = get(getFrontCellIndex(idx));
	}

	inline bool entryExists(const unsigned int x, const unsigned int y, const unsigned int z) const { return true; }

	inline unsigned int size() const { return totalCells; }

	/*void advectLinear(const double timestep, const VectorField &vectorField) {
		StlVector<T> nextVal;  nextVal.resize(totalCells);
		for (unsigned int i=0; i<numCells.x; i++)  for (unsigned int j=0; j<numCells.y; j++)  for (unsigned int k=0; k<numCells.z; k++) {
			Vector3d worldSpacePos = getPositionOfData(i,j,k);
			Vector3d advectedPos = vectorField.advect(timestep, worldSpacePos);

			if ((worldSpacePos-advectedPos).magnitudeSquared() < 0.000001) {
				unsigned int idx = getCellIndex(i,j,k);
				nextVal.vec[idx] = get(idx);
			}
			else nextVal.vec[getCellIndex(i,j,k)] = getLinear(worldSpaceToCellSpaceTransform(advectedPos));
		}
		data.deleteElements();
		data.clear();
		data = nextVal;
	}

	void advectCubic(const double timestep, const VectorField &vectorField) {
		StlVector<T> nextVal;  nextVal.resize(totalCells);
		for (unsigned int i=0; i<numCells.x; i++)  for (unsigned int j=0; j<numCells.y; j++)  for (unsigned int k=0; k<numCells.z; k++) {
			Vector3d worldSpacePos = getPositionOfData(i,j,k);
			Vector3d advectedPos = vectorField.advect(timestep, worldSpacePos);

			if ((worldSpacePos-advectedPos).magnitudeSquared() < 0.000001) {
				unsigned int idx = getCellIndex(i,j,k);
				nextVal.vec[idx] = get(idx);
			}
			else nextVal.vec[getCellIndex(i,j,k)] = getCubic(worldSpaceToCellSpaceTransform(advectedPos));
		}
		data.deleteElements();
		data.clear();
		data = nextVal;
	}*/

	void advectLinear(const Grid3D_Regular<Vector3d> &advectedWorldPositions) {
		StlVector<T> nextVal;  nextVal.resize(totalCells);
		for (unsigned int i=0; i<numCells.x; i++)  for (unsigned int j=0; j<numCells.y; j++)  for (unsigned int k=0; k<numCells.z; k++) {
			unsigned int idx = getCellIndex(i,j,k);
			Vector3d advectedPos = advectedWorldPositions.get(idx);
			Vector3d worldSpacePos = getPositionOfData(i,j,k);

			if ((worldSpacePos-advectedPos).magnitudeSquared() < 0.000001) nextVal.vec[idx] = get(idx);
			else nextVal.vec[i] = getLinear(worldSpaceToCellSpaceTransform(advectedPos));
		}
		data.deleteElements();
		data.clear();
		data = nextVal;
	}

	void advectCubic(const Grid3D_Regular<Vector3d> &advectedWorldPositions) {
		StlVector<T> nextVal;  nextVal.resize(totalCells);
		for (unsigned int i=0; i<numCells.x; i++)  for (unsigned int j=0; j<numCells.y; j++)  for (unsigned int k=0; k<numCells.z; k++) {
			unsigned int idx = getCellIndex(i,j,k);
			Vector3d advectedPos = advectedWorldPositions.get(idx);
			Vector3d worldSpacePos = getPositionOfData(i,j,k);

			if ((worldSpacePos-advectedPos).magnitudeSquared() < 0.000001) nextVal.vec[idx] = get(idx);
			else nextVal.vec[i] = getCubic(worldSpaceToCellSpaceTransform(advectedPos));
		}
		data.deleteElements();
		data.clear();
		data = nextVal;
	}

	// uses the fast marching method to extrapolate known values into the narrow band defined by max sign dist.
	void extrapolate(const Grid3D_Regular<double> &signedDistances, const double maxSignedDistance, Grid3D_Regular<bool> *isKnown, const T &unknownDefaultValue) {
		Grid3D_Regular<T> *extrapolatedVal = new Grid3D_Regular<T>(numCells, unknownDefaultValue, startPos, cellSize);
		
		BinaryHeap<HeapNode<double, Vector3ui> > minMaxHeap;
		minMaxHeap.setAsMinHeap();

		for (unsigned int i=0; i<numCells.x; i++)  for (unsigned int j=0; j<numCells.y; j++)  for (unsigned int k=0; k<numCells.z; k++) {
			if (isKnown->get(i,j,k)) {
				extrapolatedVal->set(i,j,k,get(i,j,k));

				HeapNode<double, Vector3ui> newEntry(signedDistances.get(i,j,k), Vector3ui(i,j,k));
				minMaxHeap.add(newEntry);
			}
			else extrapolatedVal->set(i,j,k,unknownDefaultValue);
		}

		// for each sorted signed distance, move its known value outward to the unknown region
		while (minMaxHeap.size() > 0) {
			HeapNode<double, Vector3ui> sdt = minMaxHeap.removeTop();

			int x = sdt.userData.x;
			int y = sdt.userData.y;
			int z = sdt.userData.z;

			double signDist = signedDistances.get(x,y,z);
			double val = extrapolatedVal->get(x,y,z);

			for (int iter=0; iter<6; iter++) {
				bool inRange;
				Vector3i idx(x,y,z);

				if		(iter==0) { inRange = x > 0;  idx.x--; }
				else if (iter==1) { inRange = x < (int)numCells.x-1;  idx.x++; }
				else if (iter==2) { inRange = y > 0;  idx.y--; }
				else if (iter==3) { inRange = y < (int)numCells.y-1;  idx.y++; }
				else if (iter==4) { inRange = z > 0;  idx.z--; }
				else/*(iter==5)*/ { inRange = z < (int)numCells.z-1;  idx.z++; }

				if (inRange && !isKnown->get(idx.x, idx.y, idx.z)) {
					double neighborDist = signedDistances.get(idx.x, idx.y, idx.z);
					if (signDist < neighborDist && neighborDist < maxSignedDistance) {
						extrapolatedVal->set(idx.x, idx.y, idx.z, extrapolatedVal->get(x,y,z));
						isKnown->set(idx.x, idx.y, idx.z, true);

						HeapNode<double, Vector3ui> newEntry(neighborDist, Vector3ui(idx.x, idx.y, idx.z));
						minMaxHeap.add(newEntry);
					}
				}
			}
		}

		minMaxHeap.clear();

		copyData(*extrapolatedVal);

		delete extrapolatedVal;  extrapolatedVal = NULL;
	}
	// the signed distance grid provided may not match this grid, so need to resample it first at this grid's data positions
	void extrapolate_resampleSignedDistance(const Grid3D_Regular<double> &signedDistances, const double maxSignedDistance, Grid3D_Regular<bool> *isKnown, const T &unknownDefaultValue) {
		Grid3D_Regular<double> *resampledSignedDistance = new Grid3D_Regular<double>(numCells, maxSignedDistance, startPos, cellSize);
		
		for (unsigned int i=0; i<numCells.x; i++)  for (unsigned int j=0; j<numCells.y; j++)  for (unsigned int k=0; k<numCells.z; k++) {
			resampledSignedDistance->set(i,j,k, signedDistances.getCubic(getPositionOfData(i,j,k)));
		}

		extrapolate(*resampledSignedDistance, maxSignedDistance, isKnown, unknownDefaultValue);

		delete resampledSignedDistance;  resampledSignedDistance = NULL;
	}

	T getMaxValue() const { T maxVal = data[0];  for (unsigned int i=1; i<totalCells; i++) if (maxVal < data[i]) maxVal = data[i];  return maxVal; }
	T getMinValue() const { T minVal = data[0];  for (unsigned int i=1; i<totalCells; i++) if (minVal > data[i]) minVal = data[i];  return minVal; }

	inline Grid3D_Regular< Vector3<T> >* getGradient() const {
		Grid3D_Regular< Vector3<T> > *grads = new Grid3D_Regular< Vector3<T> >(numCells, Vector3<T>(outsideVal,outsideVal,outsideVal), startPos, cellSize);
		for (unsigned int i=0; i<numCells.x; i++)
		for (unsigned int j=0; j<numCells.y; j++)
		for (unsigned int k=0; k<numCells.z; k++) grads->set(i,j,k, Grid3D_Regular_Base::getGradient(i,j,k));
		return grads;
	}
	inline Grid3D_Regular<T>* getGradientMagnitudes() const {
		Grid3D_Regular<T> *gradMag = new Grid3D_Regular<T>(numCells, outsideVal, startPos, cellSize);
		for (unsigned int i=0; i<numCells.x; i++)
		for (unsigned int j=0; j<numCells.y; j++)
		for (unsigned int k=0; k<numCells.z; k++) gradMag->set(i,j,k, Grid3D_Regular_Base::getGradient(i,j,k).magnitude());
		return gradMag;
	}
	inline Grid3D_Regular<T>* getGradientMagnitudesSquared() const {
		Grid3D_Regular<T> *gradMag = new Grid3D_Regular<T>(numCells, outsideVal, startPos, cellSize);
		for (unsigned int i=0; i<numCells.x; i++)
		for (unsigned int j=0; j<numCells.y; j++)
		for (unsigned int k=0; k<numCells.z; k++) gradMag->set(i,j,k, Grid3D_Regular_Base::getGradient(i,j,k).magnitudeSquared());
		return gradMag;
	}

	inline Grid3D_Regular<T>* subSampleNearestNeighbor(const unsigned int subSampleSizeX, const unsigned int subSampleSizeY, const unsigned int subSampleSizeZ) const {
		Vector3ui numberOfCells(subSampleSizeX,subSampleSizeY,subSampleSizeZ);
		numberOfCells.clamp(Vector3ui(1,1,1), numCells);
		Grid3D_Regular<T> *sample = new Grid3D_Regular<T>(numberOfCells, outsideVal);
		sample->setGridPositionUsingBoundingBox(startPos, endPos);
		for (unsigned int i=0; i<sample->numCells.x; i++)
		for (unsigned int j=0; j<sample->numCells.y; j++)
		for (unsigned int k=0; k<sample->numCells.z; k++) {
			Vector3d pos = sample->getPositionOfData(i,j,k);
			sample->set(i,j,k, getNearest(pos));
		}
		return sample;
	}
	inline Grid3D_Regular<T>* subSampleNearestNeighbor(const Vector3d sampleDistance) const
		{ Vector3ui numOfSamples = getBoundingBoxDimensions()/sampleDistance;  subSampleNearestNeighbor(numOfSamples.x, numOfSamples.y, numOfSamples.z); }

	inline Grid3D_Regular<T>* subSampleLinear(const unsigned int subSampleSizeX, const unsigned int subSampleSizeY, const unsigned int subSampleSizeZ) const {
		Vector3ui numberOfCells(subSampleSizeX,subSampleSizeY,subSampleSizeZ);
		numberOfCells.clamp(Vector3ui(1,1,1), numCells);
		Grid3D_Regular<T> *sample = new Grid3D_Regular<T>(numberOfCells, outsideVal);
		sample->setGridPositionUsingBoundingBox(startPos, endPos);
		for (unsigned int i=0; i<sample->numCells.x; i++)
		for (unsigned int j=0; j<sample->numCells.y; j++)
		for (unsigned int k=0; k<sample->numCells.z; k++) {
			Vector3d pos = sample->getPositionOfData(i,j,k);
			sample->set(i,j,k, getLinear(pos));
		}
		return sample;
	}
	inline Grid3D_Regular<T>* subSampleLinear(const Vector3d sampleDistance) const
		{ Vector3ui numOfSamples = getBoundingBoxDimensions()/sampleDistance;  subSampleLinear(numOfSamples.x, numOfSamples.y, numOfSamples.z); }

	inline Grid3D_Regular<T>* subSampleCubic(const unsigned int subSampleSizeX, const unsigned int subSampleSizeY, const unsigned int subSampleSizeZ) const {
		Vector3ui numberOfCells(subSampleSizeX,subSampleSizeY,subSampleSizeZ);
		numberOfCells.clamp(Vector3ui(1,1,1), numCells);
		Grid3D_Regular<T> *sample = new Grid3D_Regular<T>(numberOfCells, outsideVal);
		sample->setGridPositionUsingBoundingBox(startPos, endPos);
		for (unsigned int i=0; i<sample->numCells.x; i++)
		for (unsigned int j=0; j<sample->numCells.y; j++)
		for (unsigned int k=0; k<sample->numCells.z; k++) {
			Vector3d pos = sample->getPositionOfData(i,j,k);
			sample->set(i,j,k, getCubic(pos));
		}
		return sample;
	}
	inline Grid3D_Regular<T>* subSampleCubic(const Vector3d sampleDistance) const
		{ Vector3ui numOfSamples = getBoundingBoxDimensions()/sampleDistance;  subSampleCubic(numOfSamples.x, numOfSamples.y, numOfSamples.z); }

	inline Grid3D_Regular<T>* getTrimedMatrix(const float minVal, const float maxVal) const {
		Vector3ui startIdx(0,0,0), endIdx(numCells-1);

		bool allOutsideRange = true;
		for (unsigned int i=0; i<numCells.x && allOutsideRange; i++) {
			for (unsigned int j=0; j<numCells.y && allOutsideRange; j++)
			for (unsigned int k=0; k<numCells.z && allOutsideRange; k++) {
				T val = get(i,j,k);
				if (val > minVal && val < maxVal) allOutsideRange = false;
			}
			if (allOutsideRange) startIdx.x++;
		}

		allOutsideRange = true;
		for (unsigned int i=endIdx.x; i>startIdx.x && allOutsideRange; i--) {
			for (unsigned int j=0; j<numCells.y && allOutsideRange; j++)
			for (unsigned int k=0; k<numCells.z && allOutsideRange; k++) {
				T val = get(i,j,k);
				if (val > minVal && val < maxVal) allOutsideRange = false;
			}
			if (allOutsideRange) endIdx.x--;
		}

		allOutsideRange = true;
		for (unsigned int j=0; j<numCells.y && allOutsideRange; j++) {
			for (unsigned int i=startIdx.x; i<=endIdx.x && allOutsideRange; i++)
			for (unsigned int k=0; k<numCells.z && allOutsideRange; k++) {
				T val = get(i,j,k);
				if (val > minVal && val < maxVal) allOutsideRange = false;
			}
			if (allOutsideRange) startIdx.y++;
		}

		allOutsideRange = true;
		for (unsigned int j=endIdx.y; j>startIdx.y && allOutsideRange; j--) {
			for (unsigned int i=startIdx.x; i<=endIdx.x && allOutsideRange; i++)
			for (unsigned int k=0; k<numCells.z && allOutsideRange; k++) {
				T val = get(i,j,k);
				if (val > minVal && val < maxVal) allOutsideRange = false;
			}
			if (allOutsideRange) endIdx.y--;
		}

		allOutsideRange = true;
		for (unsigned int k=0; k<numCells.z && allOutsideRange; k++) {
			for (unsigned int i=startIdx.x; i<=endIdx.x && allOutsideRange; i++)
			for (unsigned int j=startIdx.y; j<=endIdx.y && allOutsideRange; j++) {
				T val = get(i,j,k);
				if (val > minVal && val < maxVal) allOutsideRange = false;
			}
			if (allOutsideRange) startIdx.z++;
		}

		allOutsideRange = true;
		for (unsigned int k=endIdx.z; k>startIdx.z && allOutsideRange; k--) {
			for (unsigned int i=startIdx.x; i<=endIdx.x && allOutsideRange; i++)
			for (unsigned int j=startIdx.y; j<=endIdx.y && allOutsideRange; j++) {
				T val = get(i,j,k);
				if (val > minVal && val < maxVal) allOutsideRange = false;
			}
			if (allOutsideRange) startIdx.z++;
		}

		Vector3ui numberOfCells = (endIdx-startIdx)+1;
		Vector3d startPosition = startPos + cellSize*Vector3d(numberOfCells.x, numberOfCells.y, numberOfCells.z);

		Grid3D_Regular<T> *trimed = new Grid3D_Regular<T>(numberOfCells, outsideVal, startPosition, cellSize);

		for (unsigned int i=0; i<numberOfCells.x; i++)
		for (unsigned int j=0; j<numberOfCells.y; j++)
		for (unsigned int k=0; k<numberOfCells.z; k++) {
			trimed->set(i,j,k, get(i+startIdx.x, j+startIdx.y, k+startIdx.z));
		}

		return trimed;
	}

	inline Grid3D_Regular<T>* getGaussianBlurredMatrix(const double sigma) const {
		Grid3D_Regular<T> *blur = new Grid2D_Regular<T>(numCells, outsideVal, startPos, cellSize);

		int radius = ceil(3.0*sigma);

		double linearScale = 1.0 / pow(M_2_PI*sigma*sigma, 1.5);
		double exponentialScale = -1.0 / (2.0*sigma*sigma);

		float blurNormalizationScale = 0;
		for (int i=-radius; i<=radius; i++)
		for (int j=-radius; j<=radius; j++)
		for (int k=-radius; k<=radius; k++) {
			blurNormalizationScale += linearScale * exp(exponentialScale * (i*i + j*j + k*k));
		}
		blurNormalizationScale = 1.0 / blurNormalizationScale;

		for (unsigned int i=0; i<numberOfCells.x; i++)
		for (unsigned int j=0; j<numberOfCells.y; j++)
		for (unsigned int k=0; k<numberOfCells.z; k++) {
			T val;

			Vector3i idx;
			bool firstIter = true;
			for (int x=(int)i-radius; x<=(int)i+radius; x++) {
				idx.x = MathUtils::clamp(x, 0, numberOfCells.x-1);

				for (int y=(int)j-radius; y<=(int)j+radius; y++) {
					idx.y = MathUtils::clamp(y, 0, numberOfCells.y-1);

					for (int z=(int)k-radius; z<=(int)k+radius; z++) {
						idx.z = MathUtils::clamp(z, 0, numberOfCells.z-1);
						
						Vector3i diff = idx-Vector2i(i,j,k);
						
						if (firstIter) {
							val = get(i,j,k) * linearScale*exp(exponentialScale * diff.magnitudeSquared());
							firstIter = false;
						}
						else val += get(i,j,k) * linearScale*exp(exponentialScale * diff.magnitudeSquared());
					}
				}
			}

			blur->set(i,j,k,val*blurNormalizationScale);
		}

		return blur;
	}

	inline double getAverageNumberOfBytesPerElement() const { return TemplateSize::getAverageNumberOfBytesPerElement(data); }
	inline double getNumberOfBytes() const { return TemplateSize::getNumberOfBytes(data); }
	inline unsigned int getNumberOfNonNullEntries() const {
		unsigned int count = 0;
		for (unsigned int i=0; i<data.size(); i++) {
			if (data[i] != outsideVal) count++;
		}
		return count;
	}

	inline void writeData(FILE *file) const { TemplateFileIO::writeBinary(file, data); }
	inline void readData(FILE *file) { TemplateFileIO::readBinary(file, &data); }

#ifdef COMPILE_WITH_CUDA
	inline void loadToCuda() {
		if (loadedToCuda) return;

		float *volumeData = new float[size()];

		unsigned int count = 0;
		for (unsigned int k=0; k<numCells.z; k++)
		for (unsigned int j=0; j<numCells.y; j++)
		for (unsigned int i=0; i<numCells.x; i++) {
			volumeData[count] = get(i,j,k);
			count++;
		}

		loadGridToCuda(volumeData, make_cudaExtent(numCells.x, numCells.y, numCells.z));

		delete[] volumeData;
		volumeData = NULL;

		loadedToCuda = true;
	}
#endif

	inline void resetCuda() {
		loadedToCuda = false;

		#ifdef COMPILE_WITH_CUDA
			cudaFreeGrid();
		#endif
	}

	StlVector<T> data;
	bool loadedToCuda;

protected:

};

template <class T> class Grid3D_Pow2 : public Grid3D_Regular<T> {
public:
	Grid3D_Pow2() { totalCells = 0; }

	Grid3D_Pow2(Vector3ui numberOfCells, T outsideValue) {
		totalCells = 0;
		outsideVal = outsideValue;
		rebuildGrid(numberOfCells);
		setGridPosition(Vector3d(0,0,0), Vector3d(1,1,1));
	}

	Grid3D_Pow2(Vector3ui numberOfCells, T outsideValue, Vector3d startPosition, Vector3d cellSizes) {
		totalCells = 0;
		outsideVal = outsideValue;
		rebuildGrid(numberOfCells);
		setGridPosition(startPosition, cellSizes);
	}

	void rebuildGrid(Vector3ui numberOfCells) { Grid3D_Regular::rebuildGrid(getPowerOf2Indices(numberOfCells)); }

	inline int getType() { return GRID_POW2; };

	Vector3ui initialNumCells;

	inline unsigned int getCellIndex(const unsigned int x, const unsigned int y, const unsigned int z) const { return MathUtils::interleaveBits(x,y,z); }

	inline unsigned int getCellIndexComponentX(const unsigned int index) const { return MathUtils::unInterleaveBitCompenentX(index); }
	inline unsigned int getCellIndexComponentY(const unsigned int index) const { return MathUtils::unInterleaveBitCompenentY(index); }
	inline unsigned int getCellIndexComponentZ(const unsigned int index) const { return MathUtils::unInterleaveBitCompenentZ(index); }
	inline Vector3ui getCellIndices(const unsigned int index) const { unsigned int x,y,z;  MathUtils::unInterleaveBits(index,x,y,z);  return Vector3ui(x,y,z); }

private:
	inline Vector3ui getPowerOf2Indices(const Vector3ui originalNumberOfCells) {
		initialNumCells = originalNumberOfCells;
		return Vector3ui(MathUtils::roundUpPowerOf2(originalNumberOfCells.x), MathUtils::roundUpPowerOf2(originalNumberOfCells.y), MathUtils::roundUpPowerOf2(originalNumberOfCells.z));
	}
};


template <> inline bool* Grid3D_Regular<bool>::getPtr(const unsigned int x, const unsigned int y, const unsigned int z) { return NULL; }


#endif