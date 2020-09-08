#include "DifferentialCompression.hpp"

namespace flopoco {
    vector<mpz_class> DifferentialCompression::getInitialTable() const {
        vector<mpz_class> reconstructedTable(1 << diffIndexSize);
        int stride = diffIndexSize - subsamplingIndexSize;
        int subsamplingShift = originalWout - subsamplingWordSize;

        for(int i = 0; i < 1 << diffIndexSize; i++)
            reconstructedTable[i] = (subsampling[i >> stride] << subsamplingShift) + diffs[i];

        return reconstructedTable;
    }

    size_t DifferentialCompression::subsamplingStorageSize() const
    {
        return subsamplingWordSize << subsamplingIndexSize;
    }

    size_t DifferentialCompression::diffsStorageSize() const
    {
        return diffWordSize << diffIndexSize;
    }
}
