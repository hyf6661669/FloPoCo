#ifndef TABLE_COMPRESSOR_HPP
#define TABLE_COMPRESSOR_HPP

#include "DifferentialCompression.hpp"

namespace flopoco {
    class TableCompressor {
    public:
        /**
          * Find a non-destructive compression of a table as a sum of a subsampling + offset
          * @param[in] values       The values of the table to compress
          * @param[in] wIn          The initial index width of the table to compress
          * @param[in] wOut         The output word size of the table to compress
          * @return                 A DifferentialCompression containing the values and the parameters of the compression for this table
          */
        static DifferentialCompression find_differential_compression(vector<mpz_class> const & values, int wIn, int wOut);

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
																				 int wIn = -1, int wOut = -1);

		};

}
#endif
