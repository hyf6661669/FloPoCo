#include <fstream>
#include <sstream>
#include "IntComplexAdder.hpp"


using namespace std;

namespace flopoco{

	extern vector<Operator *> oplist;



	IntComplexAdder::IntComplexAdder(Target* target, int wI_, int wF_, bool signedOperator_, map<string, double> inputDelays)
		: Operator(target), wI(wI_), wF(wF_), signedOperator(signedOperator_)
	{
		signedOperator ? w = 1 + wI + wF : w = wI + wF;
		
		ostringstream name;

		setCopyrightString ( "Istoan Matei, Florent de Dinechin (2008-2012)" );
		if(target->isPipelined())
			name << "IntComplexAdder_" << w << "_f"<< target->frequencyMHz() << "_uid" << getNewUId();
		else
			name << "IntComplexAdder_" << w << "_uid" << getNewUId();
		setName( name.str() );

		addInput( "Xi", 	w, true);
		addInput( "Xr", 	w, true);
		addInput( "Yi", 	w, true);
		addInput( "Yr", 	w, true);
		addInput( "Cinr", 	1		);
		addInput( "Cini", 	1		);
		addOutput("Zi", 	w, 2);
		addOutput("Zr", 	w, 2);
		
		setCriticalPath(getMaxInputDelays(inputDelays));

		IntAdder* addOperator =  new IntAdder(target, w, inDelayMap("X",getCriticalPath()));
		oplist.push_back(addOperator);
	
		inPortMap (addOperator, "X", "Xi");
		inPortMap (addOperator, "Y", "Yi");
		inPortMap (addOperator, "Cin", "Cini");
		outPortMap(addOperator, "R", "Zi", false);
		vhdl << instance(addOperator, "ADD_I");
		
		inPortMap (addOperator, "X", "Xr");
		inPortMap (addOperator, "Y", "Yr");
		inPortMap (addOperator, "Cin", "Cinr");
		outPortMap(addOperator, "R", "Zr", false);
		vhdl << instance(addOperator, "ADD_R");
		
		syncCycleFromSignal("Zr");
		setCriticalPath( addOperator->getOutputDelay("R") );
	
	}	


	IntComplexAdder::~IntComplexAdder()
	{
	}
	
	
	void IntComplexAdder::emulate(TestCase * tc)
	{
		mpz_class svXi = tc->getInputValue ( "Xi" );
		mpz_class svYi = tc->getInputValue ( "Yi" );
		mpz_class svCi = tc->getInputValue ( "Cini" );
		mpz_class svXr = tc->getInputValue ( "Xr" );
		mpz_class svYr = tc->getInputValue ( "Yr" );
		mpz_class svCr = tc->getInputValue ( "Cinr" );
		
		mpz_class svZi = svXi + svYi + svCi;
		mpz_class svZr = svXr + svYr + svCr;
		// Don't allow overflow
		mpz_clrbit ( svZi.get_mpz_t(), w );
		mpz_clrbit ( svZr.get_mpz_t(), w );
		
		tc->addExpectedOutput ( "Zi", svZi );
		tc->addExpectedOutput ( "Zr", svZr );
	}

}
