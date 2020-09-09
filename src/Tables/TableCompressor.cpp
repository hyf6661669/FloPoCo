#include "TableCompressor.hpp"
#include "Table.hpp"
#include <cassert>

namespace flopoco {

    DifferentialCompression TableCompressor::find_differential_compression(vector<mpz_class> const & values, int wIn, int wOut)
    {
        vector<mpz_class> min_max{values};
        vector<mpz_class> best_subsampling{values};
        vector<mpz_class> best_diff{};
        int diffNbEntry = 1 << wIn;
        int best_s = 0;
        int best_cost = wOut << wIn;
        int best_diff_word_size = 0;
        int best_shaved_out = 0;
        for (int s = 1 ; s < wIn ; s++) { // Iterate over each possible splitting value
            auto shift = 1 << s;
            mpz_class max_distance{0};
            int min_max_dist = ((shift >> 1) - 1);
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
            int min_width = intlog2(mpz_class{max_distance});
            int shaved_out = min_width;
            int subSampleNbEntry = 1 << (wIn - s);
            int subSampleWordSize = wOut - shaved_out;
            int subSampleLowScore = subSampleNbEntry * subSampleWordSize;
            int diffLowScore = diffNbEntry * min_width;
            int lowScore = diffLowScore + subSampleLowScore;
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
            for (int i = 0 ; i < values.size() ; i += shift) { //for each slice
                mpz_class subsample_min{min_max[i]}; //Get the minimal value of the slice as we want only positive offset
                mpz_class subsample_max{min_max[i+shift - 1]}; //Get the maximal value of the slice
                bool sub_slice_ok = false;
                while(!sub_slice_ok) {
                    mpz_class subsample_low_bits {subsample_min & low_bit_mask};
                    mpz_class cur_subsample_val {subsample_min - subsample_low_bits};
                    subsamples[i >> s] = cur_subsample_val;
                    mpz_class max_diff = subsample_max - cur_subsample_val;
                    if (max_diff >= overflow_mask) {
                            if (min_width - shaved_out < (shift - 1) ) {
                                sub_slice_ok = false;
                                shaved_out -= 1;
                                subSampleWordSize += 1;
                                subSampleLowScore += subSampleNbEntry;
                                lowScore = subSampleLowScore + diffLowScore;
                                if (lowScore >= best_cost) {
                                    goto diff_compression_bound;
                                }
                                low_bit_mask >>= 1;
                                break;
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
                    for (int j = i; j < i+shift ; j++) {
                        diff[j] = values[j] - cur_subsample_val;
                        assert(diff[j] >= 0);
                    } //end of diff bank computing loop
                } // End of shaving adjustment loop
            }// End of slices iteration
            int cur_subsample_input_word_size = 1 << (wIn - s);
            int cur_subsample_output_word_size = (wOut - shaved_out);
            int cur_diff_size = min_width << (wIn);
            int cur_cost = cur_diff_size + cur_subsample_input_word_size * cur_subsample_output_word_size;
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
        return difcompress;
    }


	string TableCompressor::newUniqueInstance(OperatorPtr op,
																								 string actualInput, string actualOutput,
																								 vector<mpz_class> values, string name,
																								 int wIn, int wOut)
	{
		if (wIn==-1 || wOut==-1){
			ostringstream o; o << " TableCompressor::newUniqueInstance needs explicit wIn and wOut" << endl;
			throw o.str();
		}
		// compress
		auto t = find_differential_compression(values, wIn, wOut);
		// generate VHDL for subsampling table
		string subsamplingIn = actualInput+"_"+name+"_subsampling";
		op->vhdl << tab << op->declare(subsamplingIn, t.subsamplingIndexSize) << " <= " << actualInput << range(wIn-1, wIn-t.subsamplingIndexSize) << ";" << endl;
		string subsamplingout = actualOutput+"_subsampling";
		string diffout = actualOutput+"_diff";
		Table::newUniqueInstance( op, subsamplingIn, subsamplingout,
															t.subsampling, name+"_subsampling",
															t.subsamplingIndexSize, t.subsamplingWordSize );
		// generate VHDL for diff table
		Table::newUniqueInstance( op, actualInput, diffout,
															t.diffs, name+"_diff",
															wIn, t.diffWordSize );
		for(int i=0; i<(1<<wIn); i++) cerr <<  "*************************** " << i << "  "  << t.diffs[i] << endl;

		int nonOverlapBits = wOut-t.diffWordSize;
		int overlapBits    = t.subsamplingWordSize - nonOverlapBits;
		// TODO an intadder when needed, but this is proably never useful
		op->vhdl << tab << op->declare(op->getTarget()->adderDelay(t.subsamplingWordSize),
																	 actualOutput+"_topbits", t.subsamplingWordSize) << " <= " << subsamplingout
			
#if 0			// Following line is with sign extension
			 			 << " + (" << rangeAssign(t.subsamplingWordSize-1, t.subsamplingWordSize-nonOverlapBits, diffout+of(t.diffWordSize-1)) << "& (" <<diffout << range(t.diffWordSize-1,t.diffWordSize-overlapBits) << "));" << endl;
#else
		<< " + (" << zg(nonOverlapBits) << "& (" <<diffout << range(t.diffWordSize-1,t.diffWordSize-overlapBits) << "));" << endl;
#endif
		op->vhdl << tab << op->declare(actualOutput, wOut) << " <= " << actualOutput+"_topbits & (" <<diffout << range(t.diffWordSize-overlapBits-1,0) << ");" << endl;
		return 	t.report(); // don't know what Operator to return, hope it is harmless here
	}

}
