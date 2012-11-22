#ifndef flopoco_random_utils_simulate_hpp
#define flopoco_random_utils_simulate_hpp

#include "Operator.hpp"

#include <set>

namespace flopoco{
namespace random{
	
mpz_class simulate(Operator *op, const std::map<std::string,mpz_class> &m)
{
	TestCase tc(op);

	if(op->getNumberOfInputs()!=(int)m.size())
		throw std::string("simulate - Operator has different number of inputs than were passed.");
	if(op->getNumberOfOutputs()!=1)
		throw std::string("simulate - Operator must have one output.");
	
	std::map<std::string,mpz_class>::const_iterator it=m.begin();
	while(it!=m.end()){
		tc.addInput(it->first, it->second);	// Will throw if input does not exist
	}
	
	op->emulate(&tc);
	
	const std::vector<mpz_class> &r=tc.getExpectedOutputValues(op->getOutputSignal(0)->getName());
	if(r.size()!=0)
		throw std::string("simulate - Operator::emulate must produce one expected value.");
	return r[0];
}

//! Simulate a unary operator
/*! Will throw an exception if:
	- There is more than one input port
	- There is more than one output port
	- The output port produces more than one output value
*/
mpz_class simulate(Operator *op, const mpz_class &x)
{
	TestCase tc(op);
	if(op->getNumberOfInputs()!=1)
		throw std::string("simulate - Operator must have one input.");
	if(op->getNumberOfOutputs()!=1)
		throw std::string("simulate - Operator must have one output.");
	tc.addInput(op->getInputSignal(0)->getName(), x);
	
	op->emulate(&tc);
	
	const std::vector<mpz_class> &r=tc.getExpectedOutputValues(op->getOutputSignal(0)->getName());
	if(r.size()!=0)
		throw std::string("simulate - Operator::emulate must produce one expected value.");
	return r[0];
}

//! Simulate a binary operator
/*! Will throw an exception if:
	- There is more than one output port
	- There are more than two input ports
	- The named ports do not exist
	- The output port produces more than one output value
*/
mpz_class simulate(Operator *op, const std::string &na, const mpz_class &va, const std::string &nb, const mpz_class &vb)
{
	TestCase tc(op);
	if(op->getNumberOfInputs()!=2)
		throw std::string("simulate - Operator must have two inputs.");
	if(op->getNumberOfOutputs()!=1)
		throw std::string("simulate - Operator must have one input.");
	if(op->getSignalByName(na)==NULL)
		throw std::string("simulate - Named input a does not exist in operator.");
	if(op->getSignalByName(nb)==NULL)
		throw std::string("simulate - Named input b does not exist in operator.");
	
	tc.addInput(na, va);
	tc.addInput(nb, vb);
	
	op->emulate(&tc);
	
	const std::vector<mpz_class> &r=tc.getExpectedOutputValues(op->getOutputSignal(0)->getName());
	if(r.size()!=0)
		throw std::string("simulate - Operator::emulate must produce one expected value (use emulate instead).");
	return r[0];
}

//! Emulates a unary operator
/*! Will throw an exception if:
	- There is more than one input port
	- There is more than one output port
*/
std::vector<mpz_class> emulate(Operator *op, const mpz_class &x)
{	
	if(op->getNumberOfInputs()!=1)
		throw std::string("emulate - Operator must have one input.");
	if(op->getNumberOfOutputs()!=1)
		throw std::string("emulate - Operator must have one input.");
	
	TestCase tc(op);
	tc.addInput(op->getInputSignal(0)->getName(), x);
		
	op->emulate(&tc);
		
	const std::vector<mpz_class> &r=tc.getExpectedOutputValues(op->getOutputSignal(0)->getName());
	
	return std::vector<mpz_class>(r.begin(), r.end());
}

//! Emulates a unary operator
/*! Will throw an exception if:
	- There is more than one input port
	- There is more than one output port
*/
std::vector<mpz_class> emulate(Operator *op, const std::vector<mpz_class> &x)
{	
	if(op->getNumberOfInputs()!=1)
		throw std::string("emulate - Operator must have one input.");
	if(op->getNumberOfOutputs()!=1)
		throw std::string("emulate - Operator must have one input.");
	
	std::set<mpz_class> acc;
	for(unsigned i=0;i<x.size();i++){
		TestCase tc(op);
		tc.addInput(op->getInputSignal(0)->getName(), x[i]);
		
		op->emulate(&tc);
		
		const std::vector<mpz_class> &r=tc.getExpectedOutputValues(op->getOutputSignal(0)->getName());
		acc.insert(r.begin(), r.end());
	}
	
	return std::vector<mpz_class>(acc.begin(), acc.end());
}



//! Emulates a binary operator
/*! Will throw an exception if:
	- There is more than one output port
	- There are more than two input ports
	- The named ports do not exist
*/
std::vector<mpz_class> emulate(Operator *op, const char *na, const char *nb, const std::vector<std::pair<mpz_class,mpz_class> > &v)
{
	if(op->getNumberOfInputs()!=2)
		throw std::string("emulate - Operator must have two inputs.");
	if(op->getNumberOfOutputs()!=1)
		throw std::string("emulate - Operator must have one input.");
	
	if(op->getSignalByName(na)==NULL)
		throw std::string("simulate - Named input a does not exist in operator.");
	if(op->getSignalByName(nb)==NULL)
		throw std::string("simulate - Named input b does not exist in operator.");
	
	std::set<mpz_class> acc;
	for(unsigned i=0;i<v.size();i++){
		TestCase tc(op);
		tc.addInput(na, v[i].first);
		tc.addInput(nb, v[i].second);
		
		op->emulate(&tc);
		
		const std::vector<mpz_class> &r=tc.getExpectedOutputValues(op->getOutputSignal(0)->getName());
		acc.insert(r.begin(), r.end());
	}
	
	return std::vector<mpz_class>(acc.begin(), acc.end());
}


}; // random
}; // flopoco

#endif
