#ifndef DIFFERENTIAL_COMPRESSION_HPP
#define DIFFERENTIAL_COMPRESSION_HPP

#include "Operator.hpp"

namespace flopoco {
	class DifferentialCompression {
	public:
		vector<mpz_class> subsampling;
		vector<mpz_class> diffs;
		int subsamplingIndexSize;
		int subsamplingWordSize;
		int diffWordSize; 
		int diffIndexSize; // also originalWin
		int originalWout;
		
				/**
          * Find a non-destructive compression of a table as a sum of a subsampling + offset
          * @param[in] values       The values of the table to compress
          * @param[in] wIn          The initial index width of the table to compress
          * @param[in] wOut         The output word size of the table to compress
          * @return                 A DifferentialCompression containing the values and the parameters of the compression for this table
          */
        static DifferentialCompression find_differential_compression(vector<mpz_class> const & values, int wIn, int wOut);

				/**
             * @brief Uncompress the table
             * @return a vector of mpz_class corresponding to the stored values
             */
            vector<mpz_class> getInitialTable() const;

            /**
             * @brief Compute the size in bits of the subsampling table content
             * @return the size as a size_t
             */
            size_t subsamplingStorageSize() const;

            /**
             * @brief Compute the size in bits of the diffs table content
             * @return the size as a size_t
             */
            size_t diffsStorageSize() const;

				string report() const;

		/**  Computes the parameters of the addition, then inserts the corresponding VHDL in Operator op. */  
		void insertAdditionVHDL(OperatorPtr op, string actualOutputName, string subsamplingOutName, string diffOutputName);
			
		
							/** This is a drop-in  replacement for Table::newUniqueInstance. 
					It will instantiate two tables and an adder, therefore will be suboptimal if the table output goes to a bit heap 
				 * @param[in] op            The Operator that will be the parent of this Table (usually "this")
				 * @param[in] actualInput   The actual input name
				 * @param[in] actualOutput  The actual input name
				 * @param[in] values        The vector of mpz_class values to be passed to the Table constructor
				 returns the report
				 */
				static string newUniqueInstance(OperatorPtr op,
																				string actualInput, string actualOutput,
																				vector<mpz_class> values, string name,
																				int wIn = -1, int wOut = -1,
																				int logicTable=0);				
			};
}
#endif
