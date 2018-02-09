#include <iostream>
#include <sstream>
#include <list>
#include <vector>
//#include <tr1/memory>

#include "gmp.h"
#include "mpfr.h"

#include "BitHeap/BitHeap.hpp"
#include "GenericBinaryPolynomial.hpp"

using namespace flopoco;



static string vhdl_string_of_monomial_option
	(const Option<MonomialOfBits>& o)
{
	ostringstream vhdl;
	if (o.is_empty()) {
		vhdl << "'0'";
		return vhdl.str();
	}
	MonomialOfBits m = o.get_value();
	bool cont = false;
	size_t i;
	for (i = 0; i < m.data.size(); i++) {
		if (m.data[i]) {
			if (cont)
				vhdl << " and ";
			vhdl << "X" << of (m.data.size() - 1 - i);
			// because x_0 is the msb in ProductIR::identity(n)
			cont = true;
		}
	}
	if (!cont)
		vhdl << "'1'";
	return vhdl.str();
}

GenericBinaryPolynomial::GenericBinaryPolynomial(Target* target,
                                                 const Product& p,
						 std::map<std::string,double>
						 	inputDelays)
	:Operator(target,inputDelays), p(p) {

	ostringstream name;
//    name << "GenericBinaryPolynomial_" << p.mon_size << "_" << p.data.size() << "_uid" << Operator::getNewUId();
    name << "GenericBinaryPolynomial_" << p.mon_size;
    setName(name.str());
	setCopyrightString("Guillaume Sergent, Florent de Dinechin 2012");

	addInput ("X" , p.mon_size);
	addOutput("R" , p.data.size());

	if (p.data.size() == 0) {
		return;
	}

    nextCycle(); //ToDo: hotfix to produce results, remove me later!

    //Creating the Bit Heap, Guillaume's compressor trees were removed
	//	shared_ptr<BitHeap> bh(new BitHeap(this, ));
	// The bit heap
	BitHeap * bitHeap = new BitHeap(this, p.data.size());


	for (unsigned i = 0; i < p.data.size(); i++) { // i is a weight
		list<MonomialOfBits>::const_iterator it = p.data[i].data.begin();
		for (; it != p.data[i].data.end(); it++) {
			ostringstream rhs;
			rhs << vhdl_string_of_monomial_option (Option<MonomialOfBits>(*it));
			bitHeap -> addBit(i, rhs.str());
		}
	}

    bitHeap -> generateCompressorVHDL();
    vhdl << tab << "R" << " <= " << bitHeap-> getSumName() << range(p.data.size()-1, 0) << ";" << endl;

};


void GenericBinaryPolynomial::emulate(TestCase * tc) {
}

void GenericBinaryPolynomial::buildStandardTestCases(TestCaseList * tcl) {
}

void GenericBinaryPolynomial::buildRandomTestCases(TestCaseList *  tcl, int n) {
}

TestCase* GenericBinaryPolynomial::buildRandomTestCases(int i) {
  TestCase* tc = new TestCase(this);
  return tc;
}
