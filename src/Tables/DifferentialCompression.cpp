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

	string DifferentialCompression::report() const
	{
		ostringstream t;
		t << "  Initial cost is:          " << originalWout << "x2^" << diffIndexSize << "=" << (originalWout << diffIndexSize) << endl;
		//t << "Initial estimated lut cost is :" << size_in_LUTs()<< endl;
		auto subsamplingCost = subsamplingWordSize << subsamplingIndexSize;
		t << "  Best subsampling cost is: " << subsamplingWordSize <<	"x2^" << subsamplingIndexSize << "=" << subsamplingCost<< endl;
		//auto lutinputs = getTarget()->lutInputs()<< endl;
		//auto lutcost = [lutinputs](int diffIndexSize, int wOut)->int {
		//								 auto effdiffIndexSize = ((diffIndexSize - lutinputs) > 0) ? diffIndexSize - lutinputs : 0;
		//								 return wOut << effwIn;
		//							 };

		// auto subsamplingLUTCost = lutcost(subsamplingIndexSize, subsamplingWordSize);
		// t << "Best subsampling LUT cost:" << subsamplingLUTCost<< endl;
		auto diffCost = diffWordSize << diffIndexSize;
		// auto diffLutCost = lutcost(diffIndexSize, diffWordSize);
		t << "  Best diff cost is:        " << diffWordSize << "x2^" << diffIndexSize << "=" << diffCost<< endl;
		// t << "Best diff LUT cost: "<< diffLutCost<< endl;
		t << "  Total compressed cost is: " << (diffCost + subsamplingCost)<< endl;
		// t << "Total LUT cost: " << (diffLutCost + subsamplingLUTCost)<< endl;
		
		// t << "Latex table line : & $" << wOut << "\\times 2^{" << diffIndexSize << "}$ & $" << (wOut << diffIndexSize) << "$ & $" <<
		// 	size_in_LUTs() << "$ & $" << diffWordSize << "\\times 2^{" << diffIndexSize << "} + " <<
		// 	subsamplingWordSize << "\\times 2^{" << subsamplingIndexSize << "}$ & $" <<
		// 	(subsamplingWordSize << subsamplingIndexSize) << "$ & $" << diffLutCost + subsamplingLUTCost <<
		// 	"$ \\\\"<< endl;
		return t.str();

	}
}
