/*
  A generic class for tables of values

  Author : Florent de Dinechin

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,
  2008-2010.
  All rights reserved.

 */

#include <algorithm>
#include <iostream>
#include <sstream>
#include <cassert>
#include <cstdlib>
#include "utils.hpp"
#include "Tables/Table.hpp"

using namespace std;
using std::begin;

namespace flopoco{

#if 0
	mpz_class Table::call_function(int x)
	{
		return function(x);
	}
#endif
	Table::DifferentialCompression Table::find_differential_compression(vector<mpz_class> const & values, int wIn, int wOut)
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
		Table::DifferentialCompression difcompress {best_subsampling, best_diff, wIn - best_s, wOut - best_shaved_out, best_diff_word_size, wIn, wOut};
		return difcompress;
	}

	Table::DifferentialCompression Table::compress() const
	{
		REPORT(INFO, "Performing differential compression on table.");
		REPORT(INFO, "Initial cost is " << wOut << "x2^" << wIn << "=" << (wOut << wIn));
		REPORT(INFO, "Initial estimated lut cost is :" << size_in_LUTs());
		DifferentialCompression ret = Table::find_differential_compression(values, wIn, wOut);
		REPORT(INFO, "Best compression split found: " << (wIn - ret.subsamplingIndexSize));
		auto subsamplingCost = ret.subsamplingWordSize << ret.subsamplingIndexSize;
		REPORT(INFO, "Best compression subsampling storage cost: " << ret.subsamplingWordSize <<
			   "x2^" << ret.subsamplingIndexSize << "=" << subsamplingCost);
		auto lutinputs = getTarget()->lutInputs();
		auto lutcost = [lutinputs](int wIn, int wOut)->int {
			auto effwIn = ((wIn - lutinputs) > 0) ? wIn - lutinputs : 0;
			return wOut << effwIn;
		};

		auto subsamplingLUTCost = lutcost(ret.subsamplingIndexSize, ret.subsamplingWordSize);
		REPORT(INFO, "Best subsampling LUT cost:" << subsamplingLUTCost);
		auto diffCost = ret.diffWordSize << ret.diffIndexSize;
		auto diffLutCost = lutcost(ret.diffIndexSize, ret.diffWordSize);
		REPORT(INFO, "Best compression diff cost: " << ret.diffWordSize << "x2^" <<
			   ret.diffIndexSize << "=" << diffCost);
		REPORT(INFO, "Best diff LUT cost: "<< diffLutCost);
		REPORT(INFO, "Total cost: " << (diffCost + subsamplingCost));
		REPORT(INFO, "Total LUT cost: " << (diffLutCost + subsamplingLUTCost));
		REPORT(INFO, "Latex table line : & $" << wOut << "\\times 2^{" << wIn << "}$ & $" << (wOut << wIn) << "$ & $" <<
			size_in_LUTs() << "$ & $" << ret.diffWordSize << "\\times 2^{" << ret.diffIndexSize << "} + " <<
			ret.subsamplingWordSize << "\\times 2^{" << ret.subsamplingIndexSize << "}$ & $" <<
			(ret.subsamplingWordSize << ret.subsamplingIndexSize) << "$ & $" << diffLutCost + subsamplingLUTCost <<
			"$ \\\\");
		return ret;
	}

	vector<mpz_class> Table::DifferentialCompression::getInitialTable() const {
		vector<mpz_class> reconstructedTable(1 << diffIndexSize);
		int stride = diffIndexSize - subsamplingIndexSize;
		int subsamplingShift = originalWout - subsamplingWordSize;

		for(int i = 0; i < 1 << diffIndexSize; i++)
			reconstructedTable[i] = (subsampling[i >> stride] << subsamplingShift) + diffs[i];

		return reconstructedTable;
	}

	size_t Table::DifferentialCompression::subsamplingStorageSize() const
	{
		return subsamplingWordSize << subsamplingIndexSize;
	}

	size_t Table::DifferentialCompression::diffsStorageSize() const
	{
		return diffWordSize << diffIndexSize;
	}


	Table::Table(OperatorPtr parentOp_, Target* target_, vector<mpz_class> _values, string _name, int _wIn, int _wOut, int _logicTable, int _minIn, int _maxIn) :
		Operator(parentOp_, target_)
	{
		srcFileName = "Table";
		setNameWithFreqAndUID(_name);
		setCopyrightString("Florent de Dinechin, Bogdan Pasca (2007-2020)");
		init(_values, _name, _wIn, _wOut,  _logicTable,  _minIn,  _maxIn);
	}


	void	Table::init(vector<mpz_class> _values, string _name,
										int _wIn, int _wOut, int _logicTable, int _minIn, int _maxIn)
	{

		values     = _values;
		wIn        = _wIn;
		wOut       = _wOut;
		minIn      = _minIn;
		maxIn      = _maxIn;
		//sanity checks: can't fill the table if there are no values to fill it with
		if(values.size() == 0)
			THROWERROR("Error in Table::init(): the set of values to be written in the table is empty" << endl);

		//set wIn
		if(wIn < 0){
			//set the value of wIn
			wIn = intlog2(values.size());
			REPORT(DEBUG, "WARNING: wIn value was not set, wIn=" << wIn << " was inferred from the vector of values");
		}
		else if(((unsigned)1<<wIn) < values.size()) {
			REPORT(DEBUG, "WARNING: wIn set to a value lower than the number of values which are to be written in the table.");
			//set the value of wIn
			wIn = intlog2(values.size());
		}

#if 0  // debug. Who never needs to print out the content of his Table ?
		for(unsigned int i=0; i<values.size(); i++)		{
			REPORT(0,"value["<<i<<"] = " << values[i]);
	}
#endif

		//determine the lowest and highest values stored in the table
		mpz_class maxValue = values[0], minValue = values[0];

		//this assumes that the values stored in the values array are all positive
		for(unsigned int i=0; i<values.size(); i++)		{
			if(values[i] < 0)
				THROWERROR("Error in table: value stored in table is negative: " << values[i] << endl);
			if(maxValue < values[i])
				maxValue = values[i];
			if(minValue > values[i])
				minValue = values[i];
		}

		//set wOut
		if(wOut < 0){
			//set the value of wOut
			wOut = intlog2(maxValue);
			REPORT(INFO, "WARNING: wOut value was not set, wOut=" << wOut << " was inferred from the vector of values");
		}
		else if(wOut < intlog2(maxValue))  {
#if 0 // This behaviour is probably a timebomb. Turns out it was fixing FixFunctionByMultipartiteTable, intconstdiv/BTCD, and probably others
			REPORT(INFO, "WARNING: wOut value was set too small. I'm fixing it up but, but I probably shouldn't");
			//set the value of wOut
			wOut = intlog2(maxValue);
#else
			THROWERROR("In Table::init(), some values require " << intlog2(maxValue) << " bits to be stored, but wOut was set to "<< wOut);
#endif
		}

		// if this is just a Table
		if(_name == "")
			setNameWithFreqAndUID(srcFileName+"_"+vhdlize(wIn)+"_"+vhdlize(wOut));

		//checks for logic table
		if(_logicTable == -1)
			logicTable = false;
		else if(_logicTable == 1)
			logicTable = true;
		else
			logicTable = (wIn <= getTarget()->lutInputs())  ||  (wOut * (mpz_class(1) << wIn) < getTarget()->sizeOfMemoryBlock()/2 );
		REPORT(DEBUG, "_logicTable=" << _logicTable << "  logicTable=" << logicTable);

		// Sanity check: the table is built using a RAM, but is underutilized
		if(!logicTable
				&& ((wIn <= getTarget()->lutInputs())  ||  (wOut*(mpz_class(1) << wIn) < 0.5*getTarget()->sizeOfMemoryBlock())))
			REPORT(0, "Warning: the table is built using a RAM block, but is underutilized");

		// Logic tables are shared by default, large tables are unique because they can have a multi-cycle schedule.
		if(logicTable)
			setShared();

		// Set up the IO signals -- this must come after the setShared()
		addInput("X", wIn, true);
		addOutput("Y", wOut, 1, true);

		//set minIn
		if(minIn < 0){
			REPORT(DEBUG, "WARNING: minIn not set, or set to an invalid value; will use the value determined from the values to be written in the table.");
			minIn = 0;
		}
		//set maxIn
		if(maxIn < 0){
			REPORT(DEBUG, "WARNING: maxIn not set, or set to an invalid value; will use the value determined from the values to be written in the table.");
			maxIn = values.size()-1;
		}

		if((1<<wIn) < maxIn){
			cerr << "ERROR: in Table constructor: maxIn too large\n";
			exit(EXIT_FAILURE);
		}

		//determine if this is a full table
		if((minIn==0) && (maxIn==(1<<wIn)-1))
			full=true;
		else
			full=false;

		//user warnings
		if(wIn > 12)
			REPORT(FULL, "WARNING: FloPoCo is building a table with " << wIn << " input bits, it will be large.");


		//create the code for the table
		REPORT(DEBUG,"Table.cpp: Filling the table");


		if(logicTable){
			int lutsPerBit;
			if(wIn < getTarget()->lutInputs())
				lutsPerBit = 1;
			else
				lutsPerBit = 1 << (wIn-getTarget()->lutInputs());
			REPORT(DETAILED, "Building a logic table that uses " << lutsPerBit << " LUTs per output bit");
		}

		cpDelay = getTarget()->tableDelay(wIn, wOut, logicTable);
		declare(cpDelay, "Y0", wOut);
		REPORT(0, "logicTable=" << logicTable << "   table delay is "<< cpDelay << "ns");

		vhdl << tab << "with X select Y0 <= " << endl;;

		for(unsigned int i=minIn.get_ui(); i<=maxIn.get_ui(); i++)
			vhdl << tab << tab << "\"" << unsignedBinary(values[i-minIn.get_ui()], wOut) << "\" when \"" << unsignedBinary(i, wIn) << "\"," << endl;
		vhdl << tab << tab << "\"";
		for(int i=0; i<wOut; i++)
			vhdl << "-";
		vhdl <<  "\" when others;" << endl;

		// TODO there seems to be several possibilities to make a BRAM; the following seems ineffective
		std::string tableAttributes;
		//set the table attributes
		if(getTarget()->getID() == "Virtex6")
			tableAttributes =  "attribute ram_extract: string;\nattribute ram_style: string;\nattribute ram_extract of Y0: signal is \"yes\";\nattribute ram_style of Y0: signal is ";
		else if(getTarget()->getID() == "Virtex5")
			tableAttributes =  "attribute rom_extract: string;\nattribute rom_style: string;\nattribute rom_extract of Y0: signal is \"yes\";\nattribute rom_style of Y0: signal is ";
		else
			tableAttributes =  "attribute ram_extract: string;\nattribute ram_style: string;\nattribute ram_extract of Y0: signal is \"yes\";\nattribute ram_style of Y0: signal is ";

		if((logicTable == 1) || (wIn <= getTarget()->lutInputs())){
			//logic
			if(getTarget()->getID() == "Virtex6")
				tableAttributes += "\"pipe_distributed\";";
			else
				tableAttributes += "\"distributed\";";
		}else{
			//block RAM
			tableAttributes += "\"block\";";
		}
		getSignalByName("Y0") -> setTableAttributes(tableAttributes);
		schedule();
		vhdl << declare("Y1", wOut) << " <= Y0; -- for the possible blockram register" << endl;

		if(!logicTable && getTarget()->registerLargeTables()){ // force a register so that a blockRAM can be infered
			setSequential();
			int cycleY0=getCycleFromSignal("Y0");
			getSignalByName("Y1") -> setSchedule(cycleY0+1, 0);
			getSignalByName("Y0") -> updateLifeSpan(1) ;
		}

		vhdl << tab << "Y <= Y1;" << endl;
	}


	Table::Table(OperatorPtr parentOp, Target* target) :
		Operator(parentOp, target){
		setCopyrightString("Florent de Dinechin, Bogdan Pasca (2007, 2018)");
	}

	mpz_class Table::val(int x){
		if(x<minIn || x>maxIn) {
			THROWERROR("Error in table: input index " << x
								 << " out of range ["<< minIn << " " << maxIn << "]" << endl);
		}
		return values[x];
	}


	int Table::size_in_LUTs() const {
		return wOut*int(intpow2(wIn-getTarget()->lutInputs()));
	}


	OperatorPtr Table::newUniqueInstance(OperatorPtr op,
																			 string actualInput, string actualOutput,
																			 vector<mpz_class> values, string name,
																			 int wIn, int wOut){
		op->schedule();
		op->inPortMap("X", actualInput);
		op->outPortMap("Y", actualOutput);
		Table* t = new Table(op, op->getTarget(), values, name, wIn, wOut);
		auto diffcompress = t->compress();
		auto decompress = diffcompress.getInitialTable();
		assert(decompress == values);
		op->vhdl << op->instance(t, name, false);
		return t;
	}



}
