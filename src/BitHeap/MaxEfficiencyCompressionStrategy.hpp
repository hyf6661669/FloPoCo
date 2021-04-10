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

		vector<float> lowerBounds;


	};

}
#endif
