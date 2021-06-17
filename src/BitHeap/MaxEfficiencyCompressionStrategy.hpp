#ifndef MAXEFFICIENCYCOMPRESSIONSTRATEGY_HPP
#define MAXEFFICIENCYCOMPRESSIONSTRATEGY_HPP

#include "BitHeap/CompressionStrategy.hpp"
#include "BitHeap/BitHeap.hpp"

namespace flopoco
{

class BitHeap;

	class MaxEfficiencyCompressionStrategy : public CompressionStrategy
	{
	public:

		/**
		 * A basic constructor for a compression strategy
		 */
		MaxEfficiencyCompressionStrategy(BitHeap *bitheap);


	private:

	    /* struct that is used for the range computation */
        struct RangeEntry {
            int range;
            int weight;
            bool isSet;
        };

		/**
		 *	@brief starts the compression algorithm. It will call maxEfficiencyAlgorithm()
		 */
		void compressionAlgorithm();

		/**
		 * generates the compressor tree
		 */
		void maxEfficiencyAlgorithm();

        /**
         * computes the number of bits that are needed to represent a number in two's complement
         */
        int reqBitsForRange2Complement(int min, int max);

        /**
         * return true if the first remainder is more efficient than the second
         */
        bool isRemainderMoreEfficient(int rem, int remToCompare);

        bool shouldPlacePseudoCompressors(vector<int> bitAmountStage);

        /**
         * chooses the better pseudo compressor for the column
         * @return the new range that this pseudo compressor creates and a bool value that is true if
         * the pseudo compressor creates inverted bits
         */
        pair<int,bool> placePseudoCompressor(int s, int column, int requiredBitsForRange, bool allowDeletion, bool useNegativeMSBValue = false);

        /**
         * creates a pseudo compressor that doe snot use sign extension for negative a negative output
         * but inverts the MSB. To be used when the sign extension is handled extra
         * @return the pseudo compressor with an inverted bit as MSB
         */
        BasicCompressor* createCompWithoutSignExtension(BasicCompressor* compressor);

        /**
         * computes the maximal range depending on the current maximal possible range.
         * assumes pseudo compressors are set at every column and every column has a height of one
         * @return new maximal range
         */
        int getMaxRangeForMaxValue(int maxValue, vector<int> currentRanges);

        /**
         * computes the maximal range for the next stage. Here the columns can be of different heights
         * @param maxValue the current maximal range
         * @param currentRanges the range that a pseudo compressor will create at each column
         * @param bitDistribution the bit distribution for the current stage
         * @param setPseudoComps the columns where pseudo comps are actually used
         * @param invertedRangeBits the columns where inverted bits are used a pseudo compressor is set
         *
         * @return new maximal range
         */
        int getMaxRangeForStage(int maxValue, vector<int> currentRanges, vector<int> bitDistribution, vector<bool> setPseudoComps, vector<bool> invertedRangeBits);

        /* recursive function to compute the range for the specified position in the actualRanges array*/
        int maxRangeForPosition(vector<RangeEntry> actualRanges, int currentPosition, int maxValue);

		vector<float> lowerBounds;

		int negativeSignExtension;
        bool needToApplyNegativeSignExtension;
	};

}
#endif
