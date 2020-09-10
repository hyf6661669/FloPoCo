#include "DiffCompressedTable.hpp"
#include "Table.hpp"
#include <cassert>

namespace flopoco {

	DiffCompressedTable::DiffCompressedTable(OperatorPtr parentOp_, Target* target_, vector<mpz_class> _values, string _name, int _wIn, int _wOut, int _logicTable, int _minIn, int _maxIn) :
		Table(parentOp_, target_)
	{
		srcFileName = "DiffCompressedTable";
		setNameWithFreqAndUID(_name);
		setCopyrightString("Florent de Dinechin, Luc Forget, Maxime Christ (2020)");
		Table::init(_values, _name, _wIn, _wOut,  _logicTable,  _minIn,  _maxIn);

		auto t = DifferentialCompression::find_differential_compression(values, wIn, wOut);
		// generate VHDL for subsampling table
		string subsamplingIn = "X_subsampling";
		vhdl << tab << declare(subsamplingIn, t.subsamplingIndexSize) << " <= " << "X" << range(wIn-1, wIn-t.subsamplingIndexSize) << ";" << endl;
		string subsamplingOut = "Y_subsampling";
		string diffout = "Y_diff";
		Table::newUniqueInstance(this, subsamplingIn, subsamplingOut,
														 t.subsampling, getName()+"_subsampling",
															t.subsamplingIndexSize, t.subsamplingWordSize, _logicTable );
		// generate VHDL for diff table
		Table::newUniqueInstance(this, "X", diffout,
															t.diffs, getName()+"_diff",
														 wIn, t.diffWordSize,
														 _logicTable);

		int nonOverlapBits = wOut-t.diffWordSize;
		int overlapBits    = t.subsamplingWordSize - nonOverlapBits;
		// TODO an intadder when needed, but this is proably never useful since this addition is so small
		vhdl << tab << declare(getTarget()->adderDelay(t.subsamplingWordSize),
													 "Y_topbits", t.subsamplingWordSize) << " <= " << subsamplingOut			
		<< " + (" << zg(nonOverlapBits) << "& (" <<diffout << range(t.diffWordSize-1,t.diffWordSize-overlapBits) << "));" << endl;
		vhdl << tab << declare("fullOut", wOut) << " <= " << "Y_topbits & (" <<diffout << range(t.diffWordSize-overlapBits-1,0) << ");" << endl;		
		vhdl << tab << "Y <= fullOut;" << endl;
	}

		OperatorPtr DiffCompressedTable::newUniqueInstance(OperatorPtr op,
																				 string actualInput, string actualOutput,
																				 vector<mpz_class> values, string name,
																				 int wIn, int wOut, int logicTable){
		op->schedule();
		op->inPortMap("X", actualInput);
		op->outPortMap("Y", actualOutput);
		DiffCompressedTable* t = new DiffCompressedTable(op, op->getTarget(), values, name, wIn, wOut,logicTable);
		op->vhdl << op->instance(t, name, false);
		return t;
	}

}
