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

        int placePseudoCompressor(int s, int column, int requiredBitsForRange, bool allowDeletion, bool useNegativeMSBValue = false);

        BasicCompressor* createCompWithoutSignExtension(BasicCompressor* compressor);

		vector<float> lowerBounds;

		int negativeSignExtension;
        bool needToApplyNegativeSignExtension;
	};

}
#endif
