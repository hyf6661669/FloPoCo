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
            int diffIndexSize;
            int originalWout;

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
      };
}
#endif
