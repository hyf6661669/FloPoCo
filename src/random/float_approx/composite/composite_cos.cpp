
#include "Operator.hpp"
#include "utils.hpp"
#include "sollya.h"

#include "random/utils/operator_factory.hpp"
#include "composite_base.hpp"

#include <cassert>

namespace flopoco
{
namespace random
{
namespace float_approx
{

class CompositeCos
	: public CompositeBase
{
private:
	int m_wDomE, m_wDomF;	
	int m_wRanE, m_wRanF;	
public:	


	CompositeCos(
			Target *target,
			int wDomE,
			int wDomF,
			map<string, double> inputDelays = emptyDelayMap
		)
		: CompositeBase(target)
		, m_wDomE(wDomE)
		, m_wDomF(wDomF)
		, m_wRanE(wDomE)
		, m_wRanF(wDomF)
	{
		if(wDomE<=0)	throw std::string("wDomE must be positive.");
		if(wDomF<=0)	throw std::string("wDomF must be positive.");
				
		std::stringstream name;
		name<<"Composite_cos_e"<<wDomE<<"_f"<<wDomF<<"_"<<getNewUId();
		setName(name.str());
		
		addFPInput("iX", m_wDomE, m_wDomF);
		addFPOutput("oY", m_wRanE, m_wRanF);
		
		assert(m_wDomE==m_wRanE);
		assert(m_wDomF==m_wRanF);

/* This is based on the approximation 4.3.99 from A&S:

cos7(x) := block([
    a2:-0.4999999963,
    a4:0.0416666418,
    a6:-0.0013888397,
    a8:0.0000247609,
    a10:-0.0000002605
    ],
    1+a2*x^2+a4*x^4+a6*x^6+a8*x^8+a10*x^10
);
*/
		setCriticalPath( getMaxInputDelays(inputDelays) );
		
		// TODO: this only does -pi/2..pi/2

		MakeFPSquarer("x2", "iX");
		MakeFPConstMult("mul10", "x2", "-0.0000002605");
		MakeFPConstAdd("add8", "mul10", "0.0000247609");
		MakeFPMultiplier("mul8", "x2", "add8");
		MakeFPConstAdd("add6", "mul8", "-0.0013888397");
		MakeFPMultiplier("mul6", "x2", "add6");
		MakeFPConstAdd("add4", "mul6", "0.0416666418");
		MakeFPMultiplier("mul4", "x2", "add4");
		MakeFPConstAdd("add2", "mul4", "-0.4999999963");
		MakeFPMultiplier("mul2", "x2", "add2");
		MakeFPConstAdd("add0", "mul2", "1.0");

		vhdl<<"oY <= add0;\n";
	}

};

static void CompositeCosFactoryUsage(std::ostream &dst)
{
	OperatorFactory::classic_OP(dst, "CompositeCos", "wE wF", false);
	dst << "    Generates an implementation of erf following A&S 4.3.97\n";
	dst << "	      (wE,wF) - Floating-point format of input and output.\n";
	dst << "         NOTE : This is for testing purposes, and should never be used.\n";
}

static Operator *CompositeCosFactoryParser(Target *target ,const std::vector<std::string> &args,int &consumed)
{
	unsigned nargs = 2;
	if (args.size()<nargs)
		throw std::string("CompositeSinFactoryParser - Not enough arguments, check usage.");
	consumed += nargs;
	
	int wDomE = atoi(args[0].c_str());
	int wDomF = atoi(args[1].c_str());
	
	return new CompositeCos(target,      wDomE, wDomF);
}

void CompositeCos_registerFactory()
{
	DefaultOperatorFactory::Register(
		"CompositeCos",
		"private",
		CompositeCosFactoryUsage,
		CompositeCosFactoryParser
	);
}

}; // float_approx
}; // random
}; // flopoco
