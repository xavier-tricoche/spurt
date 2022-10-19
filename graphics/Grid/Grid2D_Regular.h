#ifndef GRID_2D_REGULAR_H
#define GRID_2D_REGULAR_H

#include "Grid2D_Regular_Base.h"
#include "Grid3D_Pow2_Base.h"
#include "StlVectorUtils.h"
#include <algorithm>

#pragma warning(disable:4996) // D_SCL_SECURE_NO_WARNINGS for std::copy

#include "CUDA/CudaWrapper.h"
#ifdef COMPILE_WITH_CUDA
	extern "C" void loadGridToCuda(float *volumeData, cudaExtent volumeSize);
#endif

template <class T> class Grid2D_Regular : public Grid2D_Regular_Base<T> {
public:
	Grid2D_Regular() { totalCells = 0;  loadedToCuda = false; }

	Grid2D_Regular(Vector2ui numberOfCells, T outsideValue) {
		totalCells = 0;
		loadedToCuda = false;
		outsideVal = outsideValue;
		rebuildGrid(numberOfCells);
		setGridPosition(Vector2d(0,0), Vector2d(1,1));
	}

	Grid2D_Regular(Vector2ui numberOfCells, T outsideValue, Vector2d startPosition, Vector2d cellSizes) {
		totalCells = 0;
		loadedToCuda = false;
		outsideVal = outsideValue;
		rebuildGrid(numberOfCells);
		setGridPosition(startPosition, cellSizes);
	}

	template <class S> Grid2D_Regular(const Grid2D_Regular<S> &grid, T outsideVal) {
		totalCells = 0;
		loadedToCuda = false;
		outsideVal = outsideValue;
		rebuildGrid(grid.numCells);
		setGridPosition(grid.startPos, grid.cellSize);
	}

	Grid2D_Regular(const Grid2D_Regular<T> &grid, bool copyGridData = false) {
		totalCells = 0;
		loadedToCuda = false;
		outsideVal = grid.outsideVal;
		rebuildGrid(grid.numCells);
		setGridPosition(grid.startPos, grid.cellSize);
		if (copyGridData) copyData(grid);
	}

	unsigned int totalCells;

	~Grid2D_Regular() { clear(); }

	inline void clear() {
		data.clear();
		totalCells = 0;
	}

	inline void deletePointers() { data.deleteElements(); }

	inline int getType() const { return GRID2D_REGULAR; }

	virtual void rebuildGrid(Vector2ui numberOfCells) {
		int prevSize = totalCells;

		Grid2D_Regular_Base::rebuildGrid(numberOfCells, outsideVal);

		totalCells = numCells.x*numCells.y;
		if (totalCells == 0) totalCells++;

		if (totalCells != prevSize) {
			data.clear();
			data.resize(totalCells);
		}
	}

	inline void resetAllToValue(const T &defaultValue) { fill(data.vec.begin(), data.vec.end(), defaultValue); }
	inline void copyData(const Grid2D_Regular<T> &grid) {
		if (totalCells != grid.totalCells) printf("Error Grid2D_Regular::copyData: array sizes need to match\n");
		else copy(grid.data.vec.begin(), grid.data.vec.end(), data.vec.begin());
	}

	inline bool isIndexInBounds(const int index) const { return index >= 0 && index < (int)totalCells; }

	inline void clampIndex(int &index) const { clamp(index, 0, (int)totalCells-1); }

	inline void set(const unsigned int index, const T &val) { data.vec[index] = val; }
	inline void set(const unsigned int x, const unsigned int y, const T &val) { data.vec[getCellIndex(x,y)] = val; }

	inline T set_returnOld(const unsigned int index, const T &val) { T old = data[index];  data.vec[index] = val;  return old; }
	inline T set_returnOld(const unsigned int x, const unsigned int y, const T &val) { return set_returnOld(getCellIndex(x,y), val); }

	inline T get(const unsigned int index) const { return data.vec[index]; }
	inline T get(const unsigned int x, const unsigned int y) const { return data.vec[getCellIndex(x,y)]; }
	
	inline T* getPtr(const unsigned int index) { return &data.vec[index]; }
	inline T* getPtr(const unsigned int x, const unsigned int y) { return &data.vec[getCellIndex(x,y)]; }

	inline bool entryExists(const unsigned int x, const unsigned int y) const { return true; }

	inline unsigned int size() const { return totalCells; }

	/*void advectLinear(const double timestep, const VectorField &vectorField) {
		StlVector<T> nextVal;  nextVal.resize(totalCells);
		for (unsigned int i=0; i<numCells.x; i++)  for (unsigned int j=0; j<numCells.y; j++)  for (unsigned int k=0; k<numCells.z; k++) {
			Vector2d worldSpacePos = getPositionOfData(i,j,k);
			Vector2d advectedPos = vectorField.advect(timestep, worldSpacePos);

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
			Vector2d worldSpacePos = getPositionOfData(i,j,k);
			Vector2d advectedPos = vectorField.advect(timestep, worldSpacePos);

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

	void advectLinear(const Grid2D_Regular<Vector2d> &advectedWorldPositions) {
		StlVector<T> nextVal;  nextVal.resize(totalCells);
		for (unsigned int i=0; i<numCells.x; i++)  for (unsigned int j=0; j<numCells.y; j++)  for (unsigned int k=0; k<numCells.z; k++) {
			unsigned int idx = getCellIndex(i,j,k);
			Vector2d advectedPos = advectedWorldPositions.get(idx);
			Vector2d worldSpacePos = getPositionOfData(i,j,k);

			if ((worldSpacePos-advectedPos).magnitudeSquared() < 0.000001) nextVal.vec[idx] = get(idx);
			else nextVal.vec[i] = getLinear(worldSpaceToCellSpaceTransform(advectedPos));
		}
		data.deleteElements();
		data.clear();
		data = nextVal;
	}

	void advectCubic(const Grid2D_Regular<Vector2d> &advectedWorldPositions) {
		StlVector<T> nextVal;  nextVal.resize(totalCells);
		for (unsigned int i=0; i<numCells.x; i++)  for (unsigned int j=0; j<numCells.y; j++)  for (unsigned int k=0; k<numCells.z; k++) {
			unsigned int idx = getCellIndex(i,j,k);
			Vector2d advectedPos = advectedWorldPositions.get(idx);
			Vector2d worldSpacePos = getPositionOfData(i,j,k);

			if ((worldSpacePos-advectedPos).magnitudeSquared() < 0.000001) nextVal.vec[idx] = get(idx);
			else nextVal.vec[i] = getCubic(worldSpaceToCellSpaceTransform(advectedPos));
		}
		data.deleteElements();
		data.clear();
		data = nextVal;
	}


	T getMaxValue() const { T maxVal = data[0];  for (unsigned int i=1; i<totalCells; i++) if (maxVal < data[i]) maxVal = data[i];  return maxVal; }
	T getMinValue() const { T minVal = data[0];  for (unsigned int i=1; i<totalCells; i++) if (minVal > data[i]) minVal = data[i];  return minVal; }

	inline Grid2D_Regular< Vector2<T> >* getGradient() const {
		Grid2D_Regular< Vector2<T> > *grads = new Grid2D_Regular< Vector2<T> >(numCells, Vector2<T>(outsideVal,outsideVal,outsideVal), startPos, cellSize);
		for (unsigned int i=0; i<numCells.x; i++)
		for (unsigned int j=0; j<numCells.y; j++)
		for (unsigned int k=0; k<numCells.z; k++) grads->set(i,j,k, Grid2D_Regular_Base::getGradient(i,j,k));
		return grads;
	}
	inline Grid2D_Regular<T>* getGradientMagnitudes() const {
		Grid2D_Regular<T> *gradMag = new Grid2D_Regular<T>(numCells, outsideVal, startPos, cellSize);
		for (unsigned int i=0; i<numCells.x; i++)
		for (unsigned int j=0; j<numCells.y; j++)
		for (unsigned int k=0; k<numCells.z; k++) gradMag->set(i,j,k, Grid2D_Regular_Base::getGradient(i,j,k).magnitude());
		return gradMag;
	}
	inline Grid2D_Regular<T>* getGradientMagnitudesSquared() const {
		Grid2D_Regular<T> *gradMag = new Grid2D_Regular<T>(numCells, outsideVal, startPos, cellSize);
		for (unsigned int i=0; i<numCells.x; i++)
		for (unsigned int j=0; j<numCells.y; j++)
		for (unsigned int k=0; k<numCells.z; k++) gradMag->set(i,j,k, Grid2D_Regular_Base::getGradient(i,j,k).magnitudeSquared());
		return gradMag;
	}

	inline Grid2D_Regular<T>* subSampleNearestNeighbor(const unsigned int subSampleSizeX, const unsigned int subSampleSizeY, const unsigned int subSampleSizeZ) const {
		Vector2ui numberOfCells(subSampleSizeX,subSampleSizeY,subSampleSizeZ);
		numberOfCells.clamp(Vector2ui(1,1,1), numCells);
		Grid2D_Regular<T> *sample = new Grid2D_Regular<T>(numberOfCells, outsideVal);
		sample->setGridPositionUsingBoundingBox(startPos, endPos);
		for (unsigned int i=0; i<sample->numCells.x; i++)
		for (unsigned int j=0; j<sample->numCells.y; j++)
		for (unsigned int k=0; k<sample->numCells.z; k++) {
			Vector2d pos = sample->getPositionOfData(i,j,k);
			sample->set(i,j,k, getNearest(pos));
		}
		return sample;
	}
	inline Grid2D_Regular<T>* subSampleNearestNeighbor(const Vector2d sampleDistance) const
		{ Vector2ui numOfSamples = getBoundingBoxDimensions()/sampleDistance;  subSampleNearestNeighbor(numOfSamples.x, numOfSamples.y, numOfSamples.z); }

	inline Grid2D_Regular<T>* subSampleLinear(const unsigned int subSampleSizeX, const unsigned int subSampleSizeY, const unsigned int subSampleSizeZ) const {
		Vector2ui numberOfCells(subSampleSizeX,subSampleSizeY,subSampleSizeZ);
		numberOfCells.clamp(Vector2ui(1,1,1), numCells);
		Grid2D_Regular<T> *sample = new Grid2D_Regular<T>(numberOfCells, outsideVal);
		sample->setGridPositionUsingBoundingBox(startPos, endPos);
		for (unsigned int i=0; i<sample->numCells.x; i++)
		for (unsigned int j=0; j<sample->numCells.y; j++)
		for (unsigned int k=0; k<sample->numCells.z; k++) {
			Vector2d pos = sample->getPositionOfData(i,j,k);
			sample->set(i,j,k, getLinear(pos));
		}
		return sample;
	}
	inline Grid2D_Regular<T>* subSampleLinear(const Vector2d sampleDistance) const
		{ Vector2ui numOfSamples = getBoundingBoxDimensions()/sampleDistance;  subSampleLinear(numOfSamples.x, numOfSamples.y, numOfSamples.z); }

	inline Grid2D_Regular<T>* subSampleCubic(const unsigned int subSampleSizeX, const unsigned int subSampleSizeY, const unsigned int subSampleSizeZ) const {
		Vector2ui numberOfCells(subSampleSizeX,subSampleSizeY,subSampleSizeZ);
		numberOfCells.clamp(Vector2ui(1,1,1), numCells);
		Grid2D_Regular<T> *sample = new Grid2D_Regular<T>(numberOfCells, outsideVal);
		sample->setGridPositionUsingBoundingBox(startPos, endPos);
		for (unsigned int i=0; i<sample->numCells.x; i++)
		for (unsigned int j=0; j<sample->numCells.y; j++)
		for (unsigned int k=0; k<sample->numCells.z; k++) {
			Vector2d pos = sample->getPositionOfData(i,j,k);
			sample->set(i,j,k, getCubic(pos));
		}
		return sample;
	}
	inline Grid2D_Regular<T>* subSampleCubic(const Vector2d sampleDistance) const
		{ Vector2ui numOfSamples = getBoundingBoxDimensions()/sampleDistance;  subSampleCubic(numOfSamples.x, numOfSamples.y, numOfSamples.z); }

	inline Grid2D_Regular<T>* getTrimedMatrix(const float minVal, const float maxVal) const {
		Vector2ui startIdx(0,0,0), endIdx(numCells-1);

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

		Vector2ui numberOfCells = (endIdx-startIdx)+1;
		Vector2d startPosition = startPos + cellSize*Vector2d(numberOfCells.x, numberOfCells.y, numberOfCells.z);

		Grid2D_Regular<T> *trimed = new Grid2D_Regular<T>(numberOfCells, outsideVal, startPosition, cellSize);

		for (unsigned int i=0; i<numberOfCells.x; i++)
		for (unsigned int j=0; j<numberOfCells.y; j++)
		for (unsigned int k=0; k<numberOfCells.z; k++) {
			trimed->set(i,j,k, get(i+startIdx.x, j+startIdx.y, k+startIdx.z));
		}

		return trimed;
	}

	inline Grid2D_Regular<T>* getGaussianBlurredMatrix(const double sigma) const {
		Grid2D_Regular<T> *blur = new Grid2D_Regular<T>(numCells, outsideVal, startPos, cellSize);

		int radius = ceil(3.0*sigma);

		double linearScale = 1.0 / (M_2_PI*sigma*sigma);
		double exponentialScale = -1.0 / (2.0*sigma*sigma);

		float blurNormalizationScale = 0;
		for (int i=-radius; i<=radius; i++)
		for (int j=-radius; j<=radius; j++) {
			blurNormalizationScale += linearScale * exp(exponentialScale * (i*i + j*j));
		}
		blurNormalizationScale = 1.0 / blurNormalizationScale;

		for (unsigned int i=0; i<numberOfCells.x; i++)
			for (unsigned int j=0; j<numberOfCells.y; j++) {
			T val;

			Vector2i idx;
			bool firstIter = true;
			for (int x=(int)i-radius; x<=(int)i+radius; x++) {
				idx.x = MathUtils::clamp(x, 0, numberOfCells.x-1);

				for (int y=(int)j-radius; y<=(int)j+radius; y++) {
					idx.y = MathUtils::clamp(y, 0, numberOfCells.y-1);
						
					Vector2i diff = idx-Vector2i(i,j,k);
					
					if (firstIter) {
						val = get(i,j) * linearScale*exp(exponentialScale * diff.magnitudeSquared());
						firstIter = false;
					}
					else val += get(i,j) * linearScale*exp(exponentialScale * diff.magnitudeSquared());
				}
			}

			blur->set(i,j,val*blurNormalizationScale);
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

	StlVector<T> data;
	bool loadedToCuda;

protected:

};

template <class T> class Grid3D_Pow2 : public Grid2D_Regular<T> {
public:
	Grid3D_Pow2() { totalCells = 0; }

	Grid3D_Pow2(Vector2ui numberOfCells, T outsideValue) {
		totalCells = 0;
		outsideVal = outsideValue;
		rebuildGrid(numberOfCells);
		setGridPosition(Vector2d(0,0,0), Vector2d(1,1,1));
	}

	Grid3D_Pow2(Vector2ui numberOfCells, T outsideValue, Vector2d startPosition, Vector2d cellSizes) {
		totalCells = 0;
		outsideVal = outsideValue;
		rebuildGrid(numberOfCells);
		setGridPosition(startPosition, cellSizes);
	}

	void rebuildGrid(Vector2ui numberOfCells) { Grid2D_Regular::rebuildGrid(getPowerOf2Indices(numberOfCells)); }

	inline int getType() { return GRID_POW2; };

	Vector2ui initialNumCells;

	inline unsigned int getCellIndex(const unsigned int x, const unsigned int y) const { return MathUtils::interleaveBits(x,y); }

	inline unsigned int getCellIndexComponentX(const unsigned int index) const { return MathUtils::unInterleaveBitCompenentX(index); }
	inline unsigned int getCellIndexComponentY(const unsigned int index) const { return MathUtils::unInterleaveBitCompenentY(index); }
	inline unsigned int getCellIndexComponentZ(const unsigned int index) const { return MathUtils::unInterleaveBitCompenentZ(index); }
	inline Vector2ui getCellIndices(const unsigned int index) const { unsigned int x,y;  MathUtils::unInterleaveBits(index,x,y);  return Vector2ui(x,y); }

private:
	inline Vector2ui getPowerOf2Indices(const Vector2ui originalNumberOfCells) {
		initialNumCells = originalNumberOfCells;
		return Vector2ui(MathUtils::roundUpPowerOf2(originalNumberOfCells.x), MathUtils::roundUpPowerOf2(originalNumberOfCells.y), MathUtils::roundUpPowerOf2(originalNumberOfCells.z));
	}
};


template <> inline bool* Grid2D_Regular<bool>::getPtr(const unsigned int x, const unsigned int y) { return NULL; }


#endif