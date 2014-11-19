
#include "Operator.hpp"
#include "utils.hpp"
#include "sollya.h"

#include "FPMultiplier.hpp"
#include "FPSquarer.hpp"
#include "FPExp.hpp"
#include "FPDiv.hpp"
#include "ConstMult/FPConstMult.hpp"
#include "FPSquarer.hpp"
#include "FPAdderSinglePath.hpp"

#include "random/utils/operator_factory.hpp"
#include "composite_base.hpp"

#include <cassert>

namespace flopoco
{
namespace random
{
namespace float_approx
{

class CompositeNormPdf
	: public CompositeBase
{
private:
	int m_wDomE, m_wDomF;	
	int m_wRanE, m_wRanF;	
public:	
	CompositeNormPdf(
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
		name<<"Composite_normpdf_e"<<wDomE<<"_f"<<wDomF<<"_"<<getNewUId();
		setName(name.str());
		
		addFPInput("iX", m_wDomE, m_wDomF);
		addFPOutput("oY", m_wRanE, m_wRanF);
		
		assert(m_wDomE==m_wRanE);
		assert(m_wDomF==m_wRanF);

		setCriticalPath( getMaxInputDelays(inputDelays) );		

		MakeFPSquarer("x2", "iX");
		MakeFPNegate("neg_x2", "x2");
		MakeFPConstMult("half_neg_x2", "neg_x2", "1/2");
		MakeFPExp("exp_half_neg_x2", "half_neg_x2");
		MakeFPConstMult("res", "exp_half_neg_x2", "sqrt(2*pi)");

		vhdl<<"oY <= res;\n";
	}

};

static void CompositeNormPdfFactoryUsage(std::ostream &dst)
{
	OperatorFactory::classic_OP(dst, "CompositeNormPdf", "wE wF", false);
	dst << "    Generates an implementation of normpdf(x)=exp(-x^/2)/sqrt(2*pi)\n";
	dst << "	      (wE,wF) - Floating-point format of input and output.\n";
	dst << "         NOTE : This is for testing purposes, and should never be used.\n";
}

static Operator *CompositeNormPdfFactoryParser(Target *target ,const std::vector<std::string> &args,int &consumed)
{
	unsigned nargs = 2;
	if (args.size()<nargs)
		throw std::string("CompositeNormPdfFactoryParser - Not enough arguments, check usage.");
	consumed += nargs;
	
	int wDomE = atoi(args[0].c_str());
	int wDomF = atoi(args[1].c_str());
	
	return new CompositeNormPdf(target,      wDomE, wDomF);
}

void CompositeNormPdf_registerFactory()
{
	DefaultOperatorFactory::Register(
		"CompositeNormPdf",
		"private",
		CompositeNormPdfFactoryUsage,
		CompositeNormPdfFactoryParser
	);
}

}; // float_approx
}; // random
}; // flopoco
