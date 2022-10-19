#ifndef GRID_3D_REGULAR_BLOCK_H
#define GRID_3D_REGULAR_BLOCK_H

#include "Grid3D_Regular_Base.h"
#include "Grid3D_Regular.h"

#include "StlVectorUtils.h"

template <class T> class Grid3D_Regular_Block : public Grid3D_Regular_Base<T> {
public:
	Grid3D_Regular_Block(bool duplicateBorders = false) {
		duplicateDataAlongBorders = duplicateBorders;
		resetBlockSize(Vector3ui(8,8,8));
	}

	inline void resetBlockSize(Vector3ui blockSizes) {
		blockSize = blockSizes;
		invBlockSize = 1.0 / Vector3d(blockSize);
	}

	Grid3D_Regular_Block(Vector3ui numberOfCells, const T &outsideValue, bool duplicateBorders = false, Vector3ui blockSizes = Vector3ui(8,8,8)) {
		duplicateDataAlongBorders = duplicateBorders;
		outsideVal = outsideValue;
		resetBlockSize(blockSizes);
		rebuildGrid(numberOfCells);
		setGridPosition(Vector3d(0,0,0), Vector3d(1,1,1));
	}

	Grid3D_Regular_Block(Vector3ui numberOfCells, const T &outsideValue, Vector3d startPosition, Vector3d cellSizes, bool duplicateBorders = false, Vector3ui blockSizes = Vector3ui(8,8,8)) {
		duplicateDataAlongBorders = duplicateBorders;
		outsideVal = outsideValue;
		resetBlockSize(blockSizes);
		rebuildGrid(numberOfCells);
		setGridPosition(startPosition, cellSizes);
	}

	~Grid3D_Regular_Block() { clear(); }

	inline void clear() {
		blockDataPointer.clear();

		for (unsigned int i=0; i<blocks.size(); i++) {
			blocks[i]->clear();
			delete blocks[i];  blocks[i] = NULL;
		}
		blocks.clear();
	}

	inline void deletePointers() { for (unsigned int i=0; i<blocks.size(); i++) blocks[i]->deletePointers(); }

	inline int getType() const { return GRID_REGULAR_BLOCK; }

	void rebuildGrid(Vector3ui numberOfCells) {
		Grid3D_Regular_Base::rebuildGrid(numberOfCells);

		clear();

		Vector3d numBlocksf = Vector3d(numCells) / Vector3d(blockSize);
		numBlocksf.roundUp();
		numBlocks = numBlocksf;

		blockDataPointer.rebuildGrid(numBlocks);
		blockDataPointer.outsideVal = numBlocks.x*numBlocks.y*numBlocks.z;
		blockDataPointer.resetAllToValue(blockDataPointer.outsideVal);
		blockDataPointer.setGridPosition(startPos, numCells);
	}

	inline void set(const unsigned int x, const unsigned int y, const unsigned int z, const T &val) {
		Vector3ui blockIndex = getBlockNumber(x,y,z);
		Vector3ui indexWithinBlock = getIndexWithinBlock(x,y,z, blockIndex);

		unsigned int idx = blockDataPointer.get(blockIndex.x, blockIndex.y, blockIndex.z);
		if (idx == blockDataPointer.outsideVal) addBlock(blockIndex, idx);
		
		setValue(blockIndex, indexWithinBlock, idx, val);
	}

	inline T set_returnOld(const unsigned int x, const unsigned int y, const unsigned int z, const T &val) {
		Vector3ui blockIndex = getBlockNumber(x,y,z);
		Vector3ui indexWithinBlock = getIndexWithinBlock(x,y,z, blockIndex);

		unsigned int idx = blockDataPointer.get(blockIndex.x, blockIndex.y, blockIndex.z);
		if (idx == blockDataPointer.outsideVal) addBlock(blockIndex, idx);
		
		T oldVal = blocks[idx]->get(indexWithinBlock.x, indexWithinBlock.y, indexWithinBlock.z);
		setValue(blockIndex, indexWithinBlock, idx, val);
		return oldVal;
	}

	inline T get(const unsigned int x, const unsigned int y, const unsigned int z) const {
		Vector3ui blockIndex = getBlockNumber(x,y,z);
		Vector3ui indexWithinBlock = getIndexWithinBlock(x,y,z, blockIndex);

		unsigned int idx = blockDataPointer.get(blockIndex.x, blockIndex.y, blockIndex.z);
		if (idx == blockDataPointer.outsideVal) return outsideVal;

		return blocks[idx]->get(indexWithinBlock.x, indexWithinBlock.y, indexWithinBlock.z);
	}

	inline T* getPtr(const unsigned int x, const unsigned int y, const unsigned int z) {
		Vector3ui blockIndex = getBlockNumber(x,y,z);
		Vector3ui indexWithinBlock = getIndexWithinBlock(x,y,z, blockIndex);

		unsigned int idx = blockDataPointer.get(blockIndex.x, blockIndex.y, blockIndex.z);
		if (idx == blockDataPointer.outsideVal) return NULL;

		return blocks[idx]->getPtr(indexWithinBlock.x, indexWithinBlock.y, indexWithinBlock.z);
	}

	inline void get2x2x2Block(const unsigned int lowX, const unsigned int lowY, const unsigned int lowZ,
							T &val000, T &val001, T &val010, T &val011, T &val100, T &val101, T &val110, T &val111) const {
		if (duplicateDataAlongBorders) {
			Vector3ui blockIndex = getBlockNumber(lowX, lowY, lowZ);
			Vector3ui iwb = getIndexWithinBlock(lowX, lowY, lowZ, blockIndex);
			unsigned int idx = blockDataPointer.get(blockIndex.x, blockIndex.y, blockIndex.z);

			blocks[idx]->get2x2x2Block(iwb.x, iwb.y, iwb.z, val000, val001, val010, val011, val100, val101, val110, val111);
		}
		else {
			Vector3ui lowBlock = getBlockNumber(lowX, lowY, lowZ);
			Vector3ui highBlock = getBlockNumber(lowX+1, lowY+1, lowZ+1);

			if (lowBlock.x == highBlock.x) {
				if (lowBlock.y == highBlock.y) {
					if (lowBlock.z == highBlock.z) {
						Vector3ui iwb = getIndexWithinBlock(lowX, lowY, lowZ, lowBlock);
						unsigned int idx = blockDataPointer.get(lowBlock.x, lowBlock.y, lowBlock.z);

						blocks[idx]->get2x2x2Block(iwb.x, iwb.y, iwb.z, val000, val001, val010, val011, val100, val101, val110, val111);

						return;
					}
					else {
						Vector3ui iwb = getIndexWithinBlock(lowX, lowY, lowZ, lowBlock);
						unsigned int lowIdx = blockDataPointer.get(lowBlock.x, lowBlock.y, lowBlock.z);
						unsigned int highIdx = blockDataPointer.get(highBlock.x, highBlock.y, highBlock.z);

						blocks[lowIdx]->get2x2BlockXY(iwb.x, iwb.y, iwb.z, val000, val010, val100, val110);
						blocks[highIdx]->get2x2BlockXY(iwb.x, iwb.y, 0, val001, val011, val101, val111);

						return;
					}
				}
				else {
					if (lowBlock.z == highBlock.z) {
						Vector3ui iwb = getIndexWithinBlock(lowX, lowY, lowZ, lowBlock);
						unsigned int lowIdx = blockDataPointer.get(lowBlock.x, lowBlock.y, lowBlock.z);
						unsigned int highIdx = blockDataPointer.get(highBlock.x, highBlock.y, highBlock.z);

						blocks[lowIdx]->get2x2BlockXZ(iwb.x, iwb.y, iwb.z, val000, val001, val100, val101);
						blocks[highIdx]->get2x2BlockXZ(iwb.x, 0, iwb.z, val010, val011, val110, val111);

						return;
					}
				}
			}
			else {
				if (lowBlock.y == highBlock.y) {
					if (lowBlock.z == highBlock.z) {
						Vector3ui iwb = getIndexWithinBlock(lowX, lowY, lowZ, lowBlock);
						unsigned int lowIdx = blockDataPointer.get(lowBlock.x, lowBlock.y, lowBlock.z);
						unsigned int highIdx = blockDataPointer.get(highBlock.x, highBlock.y, highBlock.z);

						blocks[lowIdx]->get2x2BlockYZ(iwb.x, iwb.y, iwb.z, val000, val001, val010, val011);
						blocks[highIdx]->get2x2BlockYZ(0, iwb.y, iwb.z, val100, val101, val110, val111);

						return;
					}
				}
			}

			val000 = get(lowX, lowY, lowZ);
			val001 = get(lowX, lowY, lowZ+1);
			val010 = get(lowX, lowY+1, lowZ);
			val011 = get(lowX, lowY+1, lowZ+1);
			val100 = get(lowX+1, lowY, lowZ);
			val101 = get(lowX+1, lowY, lowZ+1);
			val110 = get(lowX+1, lowY+1, lowZ);
			val111 = get(lowX+1, lowY+1, lowZ+1);
		}
	}

	inline bool entryExists(const unsigned int x, const unsigned int y, const unsigned int z) const {
		Vector3ui blockIndex = getBlockNumber(x,y,z);
		Vector3ui indexWithinBlock = getIndexWithinBlock(x,y,z, blockIndex);

		unsigned int idx = blockDataPointer.get(blockIndex.x, blockIndex.y, blockIndex.z);
		if (idx == blockDataPointer.outsideVal) return false;

		T val = blocks[idx]->get(indexWithinBlock.x, indexWithinBlock.y, indexWithinBlock.z);
		return val != outsideVal;
	}

	inline bool entryExists(const unsigned int x, const unsigned int y, const unsigned int z, T **dataPtr) {
		// to do: make more efficient
		if (entryExists(x,y,z)) {
			*dataPtr = getPtr(x,y,z);
			return true;
		}
		else {
			*dataPtr = NULL;
			return false;
		}
	}

	inline bool getDataIndex(const unsigned int x, const unsigned int y, const unsigned int z, unsigned int &index) const {
		Vector3ui blockIndex = getBlockNumber(x,y,z);
		Vector3ui indexWithinBlock = getIndexWithinBlock(x,y,z, blockIndex);

		unsigned int idx = blockDataPointer.get(blockIndex.x, blockIndex.y, blockIndex.z);
		if (idx == blockDataPointer.outsideVal) return false;

		T val = blocks[idx]->get(indexWithinBlock.x, indexWithinBlock.y, indexWithinBlock.z);
		if (val == outsideVal) return false;

		blocks[idx]->getDataIndex(indexWithinBlock.x, indexWithinBlock.y, indexWithinBlock.z, index);
		index += idx*blocks[idx]->size();
		return true;
	}

	inline double getNumberOfBytes() const {
		double total = blockDataPointer.getNumberOfBytes();
		for (unsigned int i=0; i<blocks.size(); i++) total += blocks[i]->getNumberOfBytes();
		return total;
	}

	inline unsigned int getNumberOfNonNullEntries() const {
		unsigned int count = 0;
		for (unsigned int i=0; i<blocks.size(); i++) count += blocks[i]->getNumberOfNonNullEntries();
		return count;
	}

	inline void writeData(FILE *file) const {
		TemplateFileIO::writeBinary(file, blockSize);
		TemplateFileIO::writeBinary(file, duplicateDataAlongBorders);

		blockDataPointer.writeGridProperties(file);
		blockDataPointer.writeOutsideValue(file);
		blockDataPointer.writeData(file);

		TemplateFileIO::writeBinary(file, blocks.size());
		for (unsigned int i=0; i<blocks.size(); i++) {
			blocks[i]->writeGridProperties(file);
			blocks[i]->writeOutsideValue(file);
			blocks[i]->writeData(file);
		}
	}
	inline void readData(FILE *file) {
		TemplateFileIO::readBinary(file, &blockSize);
		resetBlockSize(blockSize);
		rebuildGrid(numCells);

		TemplateFileIO::readBinary(file, &duplicateDataAlongBorders);

		blockDataPointer.readGridProperties(file);
		blockDataPointer.readOutsideValue(file);
		blockDataPointer.readData(file);

		unsigned int size;
		TemplateFileIO::readBinary(file, &size);
		blocks.resize(size);

		TemplateFileIO::readBinary(file, &size); // for whatever reason, and additional 4-bytes are being written somewhere
												 // this seems to always be 0
		for (unsigned int i=0; i<blocks.size(); i++) {
			blocks[i] = new Grid3D_Regular<T>();
			blocks[i]->readGridProperties(file);
			blocks[i]->readOutsideValue(file);
			blocks[i]->readData(file);
		}
	}

	inline unsigned int getNumberOfBlocks() const { return (unsigned int)blocks.size(); }
	inline Vector3ui getBlockSize() const { return blockSize; }

	inline double getSparsityRatio() const {
		double val = (double)(blocks.size()*blockSize.x*blockSize.y*blockSize.z)/(double)(numCells.x*numCells.y*numCells.z);
		if (val > 1) return 1;
		else return val;
	}
	inline double getActualSparsityRatio() const {
		double entriesPerBlock = blockSize.x*blockSize.y*blockSize.z;
		double entries = 0;
		for (unsigned int i=0; i<blocks.size(); i++) {
			for (unsigned int j=0; j<blocks[i]->size(); j++) {
				if (blocks[i]->get(j) != outsideVal) entries++;
			}
		}
		return entries/(double)(numCells.x*numCells.y*numCells.z);
	}

	vector<Grid3D_Regular<T>*> blocks;
	Grid3D_Regular<unsigned int> blockDataPointer;

protected:
	Vector3ui blockSize, numBlocks;
	Vector3d invBlockSize;

	bool duplicateDataAlongBorders;

	inline Vector3ui getBlockNumber(const unsigned int x, const unsigned int y, const unsigned int z) const
		{ return Vector3ui(x*invBlockSize.x, y*invBlockSize.y, z*invBlockSize.z); }
	inline Vector3ui getIndexWithinBlock(const unsigned int x, const unsigned int y, const unsigned int z, const Vector3ui block) const
		{ return Vector3ui(x,y,z) - block*blockSize; }

	inline void addBlock(const Vector3ui &blockIndex, unsigned int &blockVectorIdx) {
		Vector3ui blockMinIdx = blockIndex*blockSize;
		Vector3d startPos = getCellMinPos(blockMinIdx.x, blockMinIdx.y, blockMinIdx.z);

		Vector3ui numCells = blockSize;
		if (duplicateDataAlongBorders) numCells = numCells+1;

		blocks.push_back(new Grid3D_Regular<T>(numCells, outsideVal, startPos, cellSize));
		blockVectorIdx = (unsigned int)blocks.size()-1;

		blockDataPointer.set(blockIndex.x, blockIndex.y, blockIndex.z, blockVectorIdx);

		blocks[blockVectorIdx]->resetAllToValue(outsideVal);
	}

	inline void setValue(const Vector3ui &blockIndex, const Vector3ui &indexWithinBlock, const unsigned int blockVectorIdx, const T &val) {
		blocks[blockVectorIdx]->set(indexWithinBlock.x, indexWithinBlock.y, indexWithinBlock.z, val);

		if (duplicateDataAlongBorders) {
			// overlap is only in the positive direction, so only need to set additional values on the negative end

			vector<Vector3ui> overlappingBlocks, iwb;
			Vector3i borders(0,0,0);
			if (indexWithinBlock.x == 0 && blockIndex.x > 0) {
				borders.x = -1;
				overlappingBlocks.push_back(Vector3ui(blockIndex.x-1, blockIndex.y, blockIndex.z));
				iwb.push_back(Vector3ui(blockSize.x, indexWithinBlock.y, indexWithinBlock.z));
			}
			if (indexWithinBlock.y == 0 && blockIndex.y > 0) {
				borders.y = -1;
				overlappingBlocks.push_back(Vector3ui(blockIndex.x, blockIndex.y-1, blockIndex.z));
				iwb.push_back(Vector3ui(indexWithinBlock.x, blockSize.y, indexWithinBlock.z));
			}
			if (indexWithinBlock.z == 0 && blockIndex.z > 0) {
				borders.z = -1;
				overlappingBlocks.push_back(Vector3ui(blockIndex.x, blockIndex.y, blockIndex.z-1));
				iwb.push_back(Vector3ui(indexWithinBlock.x, indexWithinBlock.y, blockSize.z));
			}

			if (borders.x != 0 && borders.y != 0) {
				overlappingBlocks.push_back(Vector3ui(blockIndex.x-1, blockIndex.y-1, blockIndex.z));
				iwb.push_back(Vector3ui(blockSize.x, blockSize.y, indexWithinBlock.z));
			}
			if (borders.x != 0 && borders.z != 0) {
				overlappingBlocks.push_back(Vector3ui(blockIndex.x-1, blockIndex.y, blockIndex.z-1));
				iwb.push_back(Vector3ui(blockSize.x, indexWithinBlock.y, blockSize.z));
			}
			if (borders.y != 0 && borders.z != 0) {
				overlappingBlocks.push_back(Vector3ui(blockIndex.x, blockIndex.y-1, blockIndex.z-1));
				iwb.push_back(Vector3ui(indexWithinBlock.x, blockSize.y, blockSize.z));
			}
			if (borders.x != 0 && borders.y != 0 && borders.z != 0) {
				overlappingBlocks.push_back(Vector3ui(blockIndex.x-1, blockIndex.y-1, blockIndex.z-1));
				iwb.push_back(Vector3ui(blockSize.x, blockSize.y, blockSize.z));
			}

			for (unsigned int i=0; i<overlappingBlocks.size(); i++) {
				Vector3ui blockIndex = overlappingBlocks[i];
				unsigned int idx = blockDataPointer.get(blockIndex.x, blockIndex.y, blockIndex.z);

				if (idx == blockDataPointer.outsideVal) addBlock(blockIndex, idx);

				blocks[idx]->set(iwb[i].x, iwb[i].y, iwb[i].z, val);
			}
		}
	}
};

template <> inline bool* Grid3D_Regular_Block<bool>::getPtr(const unsigned int x, const unsigned int y, const unsigned int z) { return NULL; }
template <> inline bool Grid3D_Regular_Block<bool>::entryExists(const unsigned int x, const unsigned int y, const unsigned int z, bool **dataPtr) { return entryExists(x,y,z); }



#endif