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
    };
}
#endif
