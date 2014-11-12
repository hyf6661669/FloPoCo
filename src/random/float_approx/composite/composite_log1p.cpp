
#include "Operator.hpp"
#include "utils.hpp"

#include "FPLog.hpp"
#include "FPAdderSinglePath.hpp"

#include "random/utils/operator_factory.hpp"

#include <cassert>

namespace flopoco
{
namespace random
{
namespace float_approx
{

class CompositeLog1p
	: public Operator
{
private:
	int m_wDomE, m_wDomF;	
	int m_wRanE, m_wRanF;	
public:	
	CompositeLog1p(
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
		name<<"Composite_log1p_e"<<wDomE<<"_f"<<wDomF<<"_"<<getNewUId();
		setName(name.str());
		
		addFPInput("iX", m_wDomE, m_wDomF);
		addFPOutput("oY", m_wRanE, m_wRanF);
		
		setCriticalPath( getMaxInputDelays(inputDelays) );
		
		assert(m_wDomE==m_wRanE);
		assert(m_wDomF==m_wRanF);
		
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
		
		inPortMap(pAdd, "X", "iX");
		inPortMapCst(pAdd, "Y", sNegOne);
		outPortMap(pAdd, "R", "add_R");
		vhdl<<instance(pAdd, "theAdd");
		syncCycleFromSignal("add_R", pAdd->getOutputDelay("R"));
		
		FPLog *pLog = new FPLog(target, m_wDomE, m_wDomF, 0,
		                                         inDelayMap("X", target->localWireDelay() + getCriticalPath()));
		oplist.push_back(pLog); 
		
		inPortMap(pLog, "X", "add_R");
		outPortMap(pLog, "R", "log_R");
		vhdl<<instance(pLog, "theLog");
		syncCycleFromSignal("log_R", pLog->getOutputDelay("R") );		
		
		vhdl<<"oY <= log_R;\n";
	}

};

static void CompositeLog1pFactoryUsage(std::ostream &dst)
{
	OperatorFactory::classic_OP(dst, "CompositeLog1p", "wE wF", false);
	dst << "    Generates a naive (literal) implementation of log(x+1).\n";
	dst << "	      (wE,wF) - Floating-point format of input and output.\n";
	dst << "         NOTE : This is for testing purposes, and should never be used.";
}

static Operator *CompositeLog1pFactoryParser(Target *target ,const std::vector<std::string> &args,int &consumed)
{
	unsigned nargs = 2;
	if (args.size()<nargs)
		throw std::string("CompositeLog1pFactoryParser - Not enough arguments, check usage.");
	consumed += nargs;
	
	int wDomE = atoi(args[0].c_str());
	int wDomF = atoi(args[1].c_str());
	
	return new CompositeLog1p(target,      wDomE, wDomF);
}

void CompositeLog1p_registerFactory()
{
	DefaultOperatorFactory::Register(
		"CompositeLog1p",
		"private",
		CompositeLog1pFactoryUsage,
		CompositeLog1pFactoryParser
	);
}

}; // float_approx
}; // random
}; // flopoco
