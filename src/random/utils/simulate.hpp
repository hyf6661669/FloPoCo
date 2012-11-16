


//! Simulate a unary operator
/*! Will throw an exception if:
	- There is more than one input port
	- There is more than one output port
	- The output port produces more than one output value
*/
mpz_class simulate(Operator *op, const mpz_class &x)
{
	TestCase tc(op);
	if(op->getInputCount()!=1)
		throw std::string("simulate - Operator must have one input.");
	if(op->getOutputCount()!=1)
		throw std::string("simulate - Operator must have one input.");
	tc->addInputValue(op->getInput(0)->name, x);
	
	op->emulate(tc);
	
	std::vector<mpz_class> &r=tc->getExpectedOutputs(op->getOutput(0)->name);
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
	if(op->getInputCount()!=2)
		throw std::string("simulate - Operator must have two inputs.");
	if(op->getOutputCount()!=1)
		throw std::string("simulate - Operator must have one input.");
	if(op->getSignalByName(na)==NULL)
		throw std::string("simulate - Named input a does not exist in operator.");
	if(op->getSignalByName(nb)==NULL)
		throw std::string("simulate - Named input b does not exist in operator.");
	
	tc->addInputValue(na, va);
	tc->addInputValue(nb, vb);
	
	op->emulate(tc);
	
	std::vector<mpz_class> &r=tc->getExpectedOutputs(op->getOutput(0)->name);
	if(r.size()!=0)
		throw std::string("simulate - Operator::emulate must produce one expected value.");
	return r[0];
}


//! Emulates a unary operator
/*! Will throw an exception if:
	- There is more than one input port
	- There is more than one output port
*/
std::vector<mpz_class> emulate(Operator *op, const std::vector<mpz_class> &x)
{	
	if(op->getInputCount()!=1)
		throw std::string("emulate - Operator must have one input.");
	if(op->getOutputCount()!=1)
		throw std::string("emulate - Operator must have one input.");
	
	std::set<mpz_class> acc;
	for(unsigned i=0;i<x.size();i++){
		TestCase tc(op);
		tc->addInputValue(op->getInput(0)->name);
		
		op->emulate(op);
		
		std::vector<mpz_class> &r=tc->getExpectedOutputs(op->getOutput(0)->name);
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
std::vector<mpz_class> emulate(Operator *op, const char *na, const char *nb, const std::pair<mpz_class,mpz_class> &v)
{
	if(op->getInputCount()!=2)
		throw std::string("emulate - Operator must have two inputs.");
	if(op->getOutputCount()!=1)
		throw std::string("emulate - Operator must have one input.");
	
	
	
	std::set<mpz_class> acc;
	for(unsigned i=0;i<x.size();i++){
		TestCase tc(op);
		tc->addInputValue(op->getInput(0)->name);
		
		op->emulate(op);
		
		std::vector<mpz_class> &r=tc->getExpectedOutputs(op->getOutput(0)->name);
		acc.insert(r.begin(), r.end());
	}
	
	return std::vector<mpz_class>(acc.begin(), acc.end());
}
