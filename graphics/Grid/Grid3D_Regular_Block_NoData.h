#ifndef GRID_3D_REGULAR_BLOCK_NODATA_H
#define GRID_3D_REGULAR_BLOCK_NODATA_H

#include "Grid3D_Regular_Base.h"
#include "Grid3D_Regular.h"

#include "StlVectorUtils.h"

template <class T> class Grid3D_Regular_Block_NoData : public Grid3D_Regular_Base<T> {
public:
	Grid3D_Regular_Block_NoData() {
		resetBlockSize(Vector3ui(8,8,8));
	}

	inline void resetBlockSize(Vector3ui blockSizes) {
		blockSize = blockSizes;
		if (blockSize.x > numCells.x) blockSize.x = numCells.x;
		if (blockSize.y > numCells.y) blockSize.y = numCells.y;
		if (blockSize.z > numCells.z) blockSize.z = numCells.z;

		invBlockSize = 1.0 / Vector3d(blockSize);
	}

	Grid3D_Regular_Block_NoData(Vector3ui numberOfCells, const T &outsideValue, Vector3ui blockSizes = Vector3ui(8,8,8)) {
		outsideVal = outsideValue;
		resetBlockSize(blockSizes);
		rebuildGrid(numberOfCells);
		setGridPosition(Vector3d(0,0,0), Vector3d(1,1,1));
	}

	Grid3D_Regular_Block_NoData(Vector3ui numberOfCells, const T &outsideValue, Vector3d startPosition, Vector3d cellSizes, Vector3ui blockSizes = Vector3ui(8,8,8)) {
		outsideVal = outsideValue;
		resetBlockSize(blockSizes);
		rebuildGrid(numberOfCells);
		setGridPosition(startPosition, cellSizes);
	}

	~Grid3D_Regular_Block_NoData() { clear(); }

	inline void clear() {
		blockDataPointer.clear();

		for (unsigned int i=0; i<blocks.size(); i++) {
			blocks[i]->clear();
			delete blocks[i];  blocks[i] = NULL;
		}
		blocks.clear();
	}

	inline void deletePointers() { for (unsigned int i=0; i<blocks.size(); i++) blocks[i]->deletePointers(); }

	inline int getType() const { return GRID_REGULAR_BLOCK_NODATA; }

	void rebuildGrid(Vector3ui numberOfCells) {
		Grid3D_Regular_Base::rebuildGrid(numberOfCells);

		resetBlockSize(blockSize);

		clear();

		Vector3d numBlocksf = Vector3d(numCells) / Vector3d(blockSize);
		numBlocksf.roundUp();
		numBlocks = numBlocksf;

		blockDataPointer.rebuildGrid(numBlocks);
		blockDataPointer.outsideVal = numBlocks.x*numBlocks.y*numBlocks.z;
		blockDataPointer.resetAllToValue(blockDataPointer.outsideVal);
		blockDataPointer.setGridPositionUsingBoundingBox(startPos, endPos);
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
		
		bool oldVal = blocks[idx]->get(indexWithinBlock.x, indexWithinBlock.y, indexWithinBlock.z);
		setValue(blockIndex, indexWithinBlock, idx, val);
		return oldVal ? (T)1 : (T)0;
	}

	inline T get(const unsigned int x, const unsigned int y, const unsigned int z) const {
		Vector3ui blockIndex = getBlockNumber(x,y,z);
		Vector3ui indexWithinBlock = getIndexWithinBlock(x,y,z, blockIndex);

		unsigned int idx = blockDataPointer.get(blockIndex.x, blockIndex.y, blockIndex.z);
		if (idx == blockDataPointer.outsideVal) return outsideVal;

		return blocks[idx]->get(indexWithinBlock.x, indexWithinBlock.y, indexWithinBlock.z) ? (T)1 : (T)0;
	}

	inline T* getPtr(const unsigned int x, const unsigned int y, const unsigned int z) {
		return NULL;
	}

	inline void get2x2x2Block(const unsigned int lowX, const unsigned int lowY, const unsigned int lowZ,
							T &val000, T &val001, T &val010, T &val011, T &val100, T &val101, T &val110, T &val111) const {
	}

	inline bool entryExists(const unsigned int x, const unsigned int y, const unsigned int z) const {
		Vector3ui blockIndex = getBlockNumber(x,y,z);
		Vector3ui indexWithinBlock = getIndexWithinBlock(x,y,z, blockIndex);

		unsigned int idx = blockDataPointer.get(blockIndex.x, blockIndex.y, blockIndex.z);
		if (idx == blockDataPointer.outsideVal) return false;

		return blocks[idx]->get(indexWithinBlock.x, indexWithinBlock.y, indexWithinBlock.z);
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
		else if (!blocks[idx]->get(indexWithinBlock.x, indexWithinBlock.y, indexWithinBlock.z)) return false;

		blocks[idx]->getDataIndex(indexWithinBlock.x, indexWithinBlock.y, indexWithinBlock.z, index);
		index += idx*blocks[idx]->size();
		return true;
	}

	inline double getNumberOfBytes() const {
		double total = blockDataPointer.getNumberOfBytes();
		for (unsigned int i=0; i<blocks.size(); i++) total += (double)(blockSize.x*blockSize.y*blockSize.z)/8.0;
		return total;
	}

	inline unsigned int getNumberOfNonNullEntries() const {
		unsigned int count = 0;
		for (unsigned int i=0; i<blocks.size(); i++) {
			for (unsigned int j=0; j<blocks[i]->size(); j++) {
				if (blocks[i]->get(j)) count++;
			}
		}
		return count;
	}

	inline void writeData(FILE *file) const {
		TemplateFileIO::writeBinary(file, blockSize);

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

		blockDataPointer.readGridProperties(file);
		blockDataPointer.readOutsideValue(file);
		blockDataPointer.readData(file);

		unsigned int size;
		TemplateFileIO::readBinary(file, &size);
		blocks.resize(size);

		TemplateFileIO::readBinary(file, &size); // for whatever reason, and additional 4-bytes are being written somewhere
												 // this seems to always be 0
		for (unsigned int i=0; i<blocks.size(); i++) {
			blocks[i] = new Grid3D_Regular<bool>();
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
				if (blocks[i]->get(j)) entries++;
			}
		}
		return entries/(double)(numCells.x*numCells.y*numCells.z);
	}

	void copyDataFrom(const Grid3D_Regular_Block<T> &block) {
		outsideVal = block.outsideVal;
		resetBlockSize(block.getBlockSize());
		rebuildGrid(block.numCells);
		setGridPosition(block.startPos, block.cellSize);

		blocks.resize(block.getNumberOfBlocks());
		for (unsigned int i=0; i<blocks.size(); i++) {
			blocks[i] = new Grid3D_Regular<bool>(block.blocks[i]->numCells, false, block.blocks[i]->startPos, block.blocks[i]->cellSize);
			blocks[i]->resetAllToValue(false);

			T nullVal = block.blocks[i]->outsideVal;

			for (unsigned int j=0; j<blocks[i]->size(); j++) {
				if (block.blocks[i]->get(j) != nullVal) blocks[i]->set(j, 1);
			}
		}

		for (unsigned int i=0; i<blockDataPointer.size(); i++) {
			blockDataPointer.set(i, block.blockDataPointer.get(i));
		}
	}

	vector<Grid3D_Regular<bool>*> blocks;
	Grid3D_Regular<unsigned int> blockDataPointer;

protected:
	Vector3ui blockSize, numBlocks;
	Vector3d invBlockSize;

	inline Vector3ui getBlockNumber(const unsigned int x, const unsigned int y, const unsigned int z) const
		{ return Vector3ui(x*invBlockSize.x, y*invBlockSize.y, z*invBlockSize.z); }
	inline Vector3ui getIndexWithinBlock(const unsigned int x, const unsigned int y, const unsigned int z, const Vector3ui block) const
		{ return Vector3ui(x,y,z) - block*blockSize; }

	inline void addBlock(const Vector3ui &blockIndex, unsigned int &blockVectorIdx) {
		Vector3ui blockMinIdx = blockIndex*blockSize;
		Vector3d startPos = getCellMinPos(blockMinIdx.x, blockMinIdx.y, blockMinIdx.z);

		blocks.push_back(new Grid3D_Regular<bool>(blockSize, false, startPos, cellSize));
		blockVectorIdx = (unsigned int)blocks.size()-1;

		blockDataPointer.set(blockIndex.x, blockIndex.y, blockIndex.z, blockVectorIdx);

		blocks[blockVectorIdx]->resetAllToValue(false);
	}

	inline void setValue(const Vector3ui &blockIndex, const Vector3ui &indexWithinBlock, const unsigned int blockVectorIdx, const T &val) {
		blocks[blockVectorIdx]->set(indexWithinBlock.x, indexWithinBlock.y, indexWithinBlock.z, true);
	}
};

template <> inline bool* Grid3D_Regular_Block_NoData<bool>::getPtr(const unsigned int x, const unsigned int y, const unsigned int z) { return NULL; }
template <> inline bool Grid3D_Regular_Block_NoData<bool>::entryExists(const unsigned int x, const unsigned int y, const unsigned int z, bool **dataPtr) { return entryExists(x,y,z); }



#endif