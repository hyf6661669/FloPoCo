
#include "Operator.hpp"
#include "utils.hpp"


#include "FPExp.hpp"
#include "FPDiv.hpp"
#include "FPAdderSinglePath.hpp"

#include "random/utils/operator_factory.hpp"

#include <cassert>

namespace flopoco
{
namespace random
{
namespace float_approx
{

class CompositeSigmoid
	: public Operator
{
private:
	int m_wDomE, m_wDomF;	
	int m_wRanE, m_wRanF;	
public:	
	CompositeSigmoid(
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
		name<<"Composite_sigmoid_e"<<wDomE<<"_f"<<wDomF<<"_"<<getNewUId();
		setName(name.str());
		
		mpfr_t one;
		mpfr_init(one);
		mpfr_set_d(one, 1.0, MPFR_RNDN);
		std::string sOne=std::string("\"")+fp2bin(one, m_wRanE, m_wRanF)+"\"";
		mpfr_clear(one);
		
		addFPInput("iX", m_wDomE, m_wDomF);
		addFPOutput("oY", m_wRanE, m_wRanF);
		
		setCriticalPath( getMaxInputDelays(inputDelays) );
		
		assert(m_wDomE==m_wRanE);
		assert(m_wDomF==m_wRanF);
		
		vhdl << declare("negX",wDomE+wDomF+3) << " <= iX("<<(2+wDomE+wDomF)<<" downto "<<(1+wDomE+wDomF)<<")&(not iX("<<(wDomE+wDomF)<<"))&iX("<<(wDomE+wDomF-1)<<" downto 0);\n";
		
		
		FPExp *pExp = new FPExp(target, m_wDomE, m_wDomF, 0, 0, -1, false, 0.7,
		                                         inDelayMap("X", target->localWireDelay() + getCriticalPath()));
		oplist.push_back(pExp); 
		
		inPortMap(pExp, "X", "negX");
		outPortMap(pExp, "R", "exp_R");
		vhdl<<instance(pExp, "theExp");
		syncCycleFromSignal("exp_R", pExp->getOutputDelay("R") );
		
		
		Operator *pAdd = new FPAdderSinglePath(target,
			m_wDomE, m_wDomF, // Input A
			m_wDomE, m_wDomF, // Input B (constant)
			m_wDomE, m_wDomF,
			inDelayMap("X", target->localWireDelay() + getCriticalPath()));
		oplist.push_back(pAdd);
	
		inPortMap(pAdd, "X", "exp_R");
		inPortMapCst(pAdd, "Y", sOne);
		outPortMap(pAdd, "R", "add_R");
		vhdl<<instance(pAdd, "theAdd");
		syncCycleFromSignal("add_R", pAdd->getOutputDelay("R"));
		
		Operator *pDiv=new FPDiv(target, m_wDomE, m_wDomF);
		oplist.push_back(pDiv);
		
		inPortMapCst(pDiv, "X", sOne);
		inPortMap(pDiv, "Y", "add_R");
		outPortMap(pDiv, "R", "inv_R");
		vhdl<<instance(pDiv, "theInv");
		syncCycleFromSignal("inv_R", pDiv->getOutputDelay("R"));
		
		vhdl<<"oY <= inv_R;\n";
	}

};

static void CompositeSigmoidFactoryUsage(std::ostream &dst)
{
	OperatorFactory::classic_OP(dst, "CompositeSigmoid", "wE wF", false);
	dst << "    Generates a naive (literal) implementation of log(x+1).\n";
	dst << "	      (wE,wF) - Floating-point format of input and output.\n";
	dst << "         NOTE : This is for testing purposes, and should never be used.\n";
}

static Operator *CompositeSigmoidFactoryParser(Target *target ,const std::vector<std::string> &args,int &consumed)
{
	unsigned nargs = 2;
	if (args.size()<nargs)
		throw std::string("CompositeSigmoidFactoryParser - Not enough arguments, check usage.");
	consumed += nargs;
	
	int wDomE = atoi(args[0].c_str());
	int wDomF = atoi(args[1].c_str());
	
	return new CompositeSigmoid(target,      wDomE, wDomF);
}

void CompositeSigmoid_registerFactory()
{
	DefaultOperatorFactory::Register(
		"CompositeSigmoid",
		"private",
		CompositeSigmoidFactoryUsage,
		CompositeSigmoidFactoryParser
	);
}

}; // float_approx
}; // random
}; // flopoco
