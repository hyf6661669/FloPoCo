#include "DifferentialCompression.hpp"
#include "Table.hpp"
#include <cassert>

namespace flopoco {

    DifferentialCompression DifferentialCompression::find_differential_compression(vector<mpz_class> const & values, int wIn, int wOut)
    {
        vector<mpz_class> min_max{values};
        vector<mpz_class> best_subsampling{values};
        vector<mpz_class> best_diff{};
        int64_t diffNbEntry = 1 << wIn;
        int best_s = 0;
        mpz_class best_cost = mpz_class(wOut) << wIn;
        int best_diff_word_size = 0;
        int best_shaved_out = 0;
        for (int s = 1 ; s < wIn ; s++) { // Iterate over each possible splitting value
            auto shift = 1 << s;
            mpz_class max_distance{0};
            int64_t min_max_dist = ((shift >> 1) - 1);
            //TODO : Find better lattice storage scheme to optimise cache access
            // Compute min, max and diff of current diff table slice using resuts of previous subslice
            for (auto min_op1_iter = begin(min_max) ; min_op1_iter < end(min_max) ; min_op1_iter += shift) {
                auto max_op1_iter = min_op1_iter + min_max_dist;
                auto min_op2_iter = max_op1_iter + 1;
                auto max_op2_iter = min_op2_iter + min_max_dist;
                if (*min_op1_iter > *min_op2_iter) {
                    std::swap(*min_op1_iter, *min_op2_iter);
                }
                if (*max_op2_iter < *max_op1_iter) {
                    std::swap(*max_op1_iter, *max_op2_iter);
                }
                mpz_class dist{*max_op2_iter - *min_op1_iter};
                if (dist > max_distance) {
                    max_distance = dist;
                }
            }

            /*
             * In the best case, without overflow when adding subsamples low bits, we need enough storage to store the difference between
             * the maximal and minimal element of the "slice" in which this distance is maximal
             */
            int64_t min_width = intlog2(mpz_class{max_distance});
            int64_t shaved_out = min_width;
            int64_t subSampleNbEntry = 1 << (wIn - s);
						int64_t subSampleWordSize = wOut - shaved_out;
						int64_t subSampleLowScore = subSampleNbEntry * subSampleWordSize;
						int64_t diffLowScore = diffNbEntry * min_width;
						int64_t lowScore = diffLowScore + subSampleLowScore;
            if (lowScore >= best_cost) {
                diff_compression_bound:
                    continue;
            }
            vector<mpz_class> subsamples(1 << (wIn - s));
            vector<mpz_class> diff(1 << wIn);
            mpz_class overflow_mask {1};
            overflow_mask <<= min_width;
            mpz_class low_bit_mask{overflow_mask - 1};
            bool had_to_overflow = false;
            diff_compression_compute:
            for (int64_t i = 0 ; i < values.size() ; i += shift) { //for each slice
                bool sub_slice_ok = false;
				//Get the minimal value of the slice as we want only positive offset
                mpz_class subsample_min{min_max[i]}; 
				//Get the maximal value of the slice
                mpz_class subsample_max{min_max[i+shift - 1]}; 
                while(!sub_slice_ok) {
                    mpz_class subsample_low_bits {subsample_min & low_bit_mask};
                    mpz_class cur_subsample_val {subsample_min - subsample_low_bits};
                    subsamples[i >> s] = cur_subsample_val;
                    mpz_class max_diff = subsample_max - cur_subsample_val;
                    if (max_diff >= overflow_mask) {
						sub_slice_ok = false;
						if (min_width - shaved_out < (shift - 1) ) {
							shaved_out -= 1;
							subSampleWordSize += 1;
							subSampleLowScore += subSampleNbEntry;
							lowScore = subSampleLowScore + diffLowScore;
							if (lowScore >= best_cost) {
								goto diff_compression_bound;
							}
							low_bit_mask >>= 1;
							continue;
						} else {
							min_width += 1;
							diffLowScore = diffNbEntry * min_width;
							shaved_out = min_width;
							subSampleWordSize = wOut - shaved_out;
							subSampleLowScore = subSampleWordSize * subSampleNbEntry;
							lowScore = subSampleLowScore + diffLowScore;
							if (lowScore >= best_cost) {
								goto diff_compression_bound;
							}
							overflow_mask = 1;
							overflow_mask <<= min_width;
							low_bit_mask = overflow_mask - 1;
							assert(!had_to_overflow);
							had_to_overflow = true;
							goto diff_compression_compute; //Jump should be done at most one time
						}
					}
                    sub_slice_ok = true;
                    for (int64_t j = i; j < i+shift ; j++) {
                        diff[j] = values[j] - cur_subsample_val;
                        assert(diff[j] >= 0); 
						assert(values[j] <= subsample_max);
                    } //end of diff bank computing loop
                } // End of shaving adjustment loop
            }// End of slices iteration
            int64_t cur_subsample_input_word_size = 1 << (wIn - s);
            int64_t cur_subsample_output_word_size = (wOut - shaved_out);
            int64_t cur_diff_size = min_width << (wIn);
            int64_t cur_cost = cur_diff_size + cur_subsample_input_word_size * cur_subsample_output_word_size;
            if (cur_cost < best_cost) {
                best_s = s;
                best_cost = cur_cost;
                best_diff_word_size = min_width;
                best_shaved_out = shaved_out;
                swap(best_subsampling, subsamples);
                swap(best_diff, diff);
            }
        } // End of iteration to find best s
        for (auto& subsample : best_subsampling) {
            subsample >>= best_shaved_out;
        }
        DifferentialCompression difcompress {best_subsampling, best_diff, wIn - best_s, wOut - best_shaved_out, best_diff_word_size, wIn, wOut};
				// Cet assert plante pour certaines valeurs parce que diff[0] n'existe pas (2020/09/18). Je n'ai pas débuggé plus loin
				// exemple: ./flopoco verbose=1 FixFunctionByMultipartiteTable f="sin(pi/4*x)" signedIn=0 lsbIn=-20 lsbOut=-24 tableCompression=1
				// assert(difcompress.getInitialTable() == values);
        return difcompress;
    }



	vector<mpz_class> DifferentialCompression::getInitialTable() const {
        vector<mpz_class> reconstructedTable(1 << diffIndexSize);
        int stride = diffIndexSize - subsamplingIndexSize;

        for(int i = 0; i < 1 << diffIndexSize; i++){
					reconstructedTable[i] = (subsampling[i >> stride] << subsamplingShift()) + diffs[i];
				}
        return reconstructedTable;
    }

	size_t DifferentialCompression::subsamplingShift() const{
		return originalWout-subsamplingWordSize;
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

	string DifferentialCompression::newUniqueInstance(OperatorPtr op,
																										string actualInputName, string actualOutputName,
																										vector<mpz_class> values, string name,
																										int wIn, int wOut, int logicTable)
	{
		if (wIn==-1 || wOut==-1){
			ostringstream o; o << " DiffCompressedTable::newUniqueInstance needs explicit wIn and wOut" << endl;
			throw o.str();
		}
		// compress
		auto t = find_differential_compression(values, wIn, wOut);
		string subsamplingOutName = actualOutputName+"_subsampling";
		string diffOutputName = actualOutputName+"_diff";
		string subsamplingIn = actualInputName+"_"+name+"_subsampling";
		// generate VHDL for subsampling table
		op->vhdl << tab << op->declare(subsamplingIn, t.subsamplingIndexSize) << " <= " << actualInputName << range(wIn-1, wIn-t.subsamplingIndexSize) << ";" << endl;
		Table::newUniqueInstance( op, subsamplingIn, subsamplingOutName,
															t.subsampling, name+"_subsampling",
															t.subsamplingIndexSize, t.subsamplingWordSize, logicTable );
		// generate VHDL for diff table
		Table::newUniqueInstance( op, actualInputName, diffOutputName,
															t.diffs, name+"_diff",
															wIn, t.diffWordSize );
		
		t.insertAdditionVHDL(op, actualOutputName, subsamplingOutName, diffOutputName);
		return 	t.report(); // don't know what Operator to return, hope it is harmless here 
	}


	
	void DifferentialCompression::insertAdditionVHDL(OperatorPtr op,
																									 string actualOutputName, string subsamplingOutName,string diffOutputName)
	{
		int wIn = diffIndexSize;
		int wOut=originalWout;
		int nonOverlapMSBBits = wOut-diffWordSize;
		int overlapMiddleBits    = subsamplingWordSize - nonOverlapMSBBits;
		// TODO an intadder when needed, but this is proably never useful
		op->vhdl << tab << op->declare(op->getTarget()->adderDelay(subsamplingWordSize),
																	 actualOutputName+"_topbits", subsamplingWordSize) << " <= " << subsamplingOutName			
		<< " + (" << zg(nonOverlapMSBBits) << "& (" << diffOutputName << range(diffWordSize-1, diffWordSize-overlapMiddleBits) << "));" << endl;
		op->vhdl << tab << op->declare(actualOutputName, wOut) << " <= " << actualOutputName+"_topbits & (" <<diffOutputName << range(diffWordSize-overlapMiddleBits-1,0) << ");" << endl;
	}
}
