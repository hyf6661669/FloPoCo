
#include "Operator.hpp"
#include "utils.hpp"

#include "FPExp.hpp"
#include "FPAdderSinglePath.hpp"

#include "random/utils/operator_factory.hpp"

#include <cassert>

namespace flopoco
{
namespace random
{
namespace float_approx
{

class CompositeExpm1
	: public Operator
{
private:
	int m_wDomE, m_wDomF;	
	int m_wRanE, m_wRanF;	
public:	
	CompositeExpm1(
			Target *target,
			int wDomE,
			int wDomF,
			map<string, double> inputDelays = emptyDelayMap
		)
		: Operator(target)
		, m_wDomE(wDomE)
		, m_wDomF(wDomF)
		, m_wRanE(wDomE)
		, m_wRanF(wDomF)
	{
		if(wDomE<=0)	throw std::string("wDomE must be positive.");
		if(wDomF<=0)	throw std::string("wDomF must be positive.");
				
		std::stringstream name;
		name<<"Composite_expm1_e"<<wDomE<<"_f"<<wDomF<<"_"<<getNewUId();
		setName(name.str());
		
		addFPInput("iX", m_wDomE, m_wDomF);
		addFPOutput("oY", m_wRanE, m_wRanF);
		
		setCriticalPath( getMaxInputDelays(inputDelays) );
		
		assert(m_wDomE==m_wRanE);
		assert(m_wDomF==m_wRanF);
		
		FPExp *pExp = new FPExp(target, m_wDomE, m_wDomF, 0, 0, -1, false, 0.7, 
		                                         inDelayMap("X", target->localWireDelay() + getCriticalPath()));
		oplist.push_back(pExp); 
		
		inPortMap(pExp, "X", "iX");
		outPortMap(pExp, "R", "exp_R");
		vhdl<<instance(pExp, "theExp");
		syncCycleFromSignal("exp_R", pExp->getOutputDelay("R") );
		
		Operator *pAdd = new FPAdderSinglePath(target,
			m_wDomE, m_wDomF, // Input A
			m_wDomE, m_wDomF, // Input B (constant)
			m_wDomE, m_wDomF,
			inDelayMap("X", target->localWireDelay() + getCriticalPath()));
		oplist.push_back(pAdd);
		
		mpfr_t negOne;
		mpfr_init(negOne);
		mpfr_set_d(negOne, -1.0, MPFR_RNDN);
		std::string sNegOne=std::string("\"")+fp2bin(negOne, m_wRanE, m_wRanF)+"\"";
		mpfr_clear(negOne);
		
		inPortMap(pAdd, "X", "exp_R");
		inPortMapCst(pAdd, "Y", sNegOne);
		outPortMap(pAdd, "R", "add_R");
		vhdl<<instance(pAdd, "theAdd");
		syncCycleFromSignal("add_R", pAdd->getOutputDelay("R"));
		
		vhdl<<"oY <= add_R;\n";
	}

};

static void CompositeExpm1FactoryUsage(std::ostream &dst)
{
	OperatorFactory::classic_OP(dst, "CompositeExpm1", "wE wF", false);
	dst << "    Generates a naive (literal) implementation of exp(x)-1.\n";
	dst << "	      (wE,wF) - Floating-point format of input and output.\n";
	dst << "         NOTE : This is for testing purposes, and should never be used.";
}

static Operator *CompositeExpm1FactoryParser(Target *target ,const std::vector<std::string> &args,int &consumed)
{
	unsigned nargs = 2;
	if (args.size()<nargs)
		throw std::string("CompositeExpm1FactoryParser - Not enough arguments, check usage.");
	consumed += nargs;
	
	int wDomE = atoi(args[0].c_str());
	int wDomF = atoi(args[1].c_str());
	
	return new CompositeExpm1(target,      wDomE, wDomF);
}

void CompositeExpm1_registerFactory()
{
	DefaultOperatorFactory::Register(
		"CompositeExpm1",
		"private",
		CompositeExpm1FactoryUsage,
		CompositeExpm1FactoryParser
	);
}

}; // float_approx
}; // random
}; // flopoco
