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
		string diffOut = "Y_diff";
		Table::newUniqueInstance(this, subsamplingIn, subsamplingOut,
														 t.subsampling, getName()+"_subsampling",
															t.subsamplingIndexSize, t.subsamplingWordSize, _logicTable );
		// generate VHDL for diff table
		Table::newUniqueInstance(this, "X", diffOut,
															t.diffs, getName()+"_diff",
														 wIn, t.diffWordSize,
														 _logicTable);

		t.insertAdditionVHDL(this, "fullOut", subsamplingOut, diffOut);
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
