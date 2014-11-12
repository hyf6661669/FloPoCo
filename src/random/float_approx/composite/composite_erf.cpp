
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

class CompositeErf
	: public CompositeBase
{
private:
	int m_wDomE, m_wDomF;	
	int m_wRanE, m_wRanF;	
public:	
/*
	void MakeFPNegate(std::string dst, std::string src)
	{
		Signal *pSrc=getSignalByName(src);
		if(pSrc==NULL)
			throw std::string("MakeFPNegate - Attempt to use non-existant signal as source.");
		if(!pSrc->isFP())
			throw std::string("MakeFPNegate - Attempt to use non-FP signal as source.");
		int wE=pSrc->wE();
		int wF=pSrc->wF();
		
		vhdl<<declareFP(dst,wE,wF)<<" <= ("<<src<<"("<<(2+wE+wF)<<" downto "<<(1+wE+wF)<<")&(not "<<src<<"("<<(wE+wF)<<"))&"<<src<<"("<<(wE+wF-1)<<" downto 0));\n";
	}
	
	// Output will be "name", instance will be "name"+"_inst"
	Operator *MakeFPExp(std::string name, std::string src)
	{
		Signal *pSrc=getSignalByName(src);
		if(pSrc==NULL)
			throw std::string("MakeFPExp - Attempt to use non-existant signal as source.");
		if(!pSrc->isFP())
			throw std::string("MakeFPExp - Attempt to use non-FP signal as source.");
		int wE=pSrc->wE();
		int wF=pSrc->wF();
		
		FPExp *pExp = new FPExp(getTarget(), wE, wF, 0, 0, -1, false, 0.7, inDelayMap("X", getTarget()->localWireDelay() + getCriticalPath()));
		oplist.push_back(pExp); 
		
		inPortMap(pExp, "X", src);
		outPortMap(pExp, "R", name+"_R");
		vhdl<<instance(pExp, name+"_inst");
		syncCycleFromSignal(name+"_R", pExp->getOutputDelay("R") );	

		vhdl<<declareFP(name,wE,wF)<<" <= "<<name<<"_R;\n";

		return pExp;
	}

	// Output will be "name", instance will be "name"+"_inst"
    Operator *MakeFPMultiplier(std::string name, std::string srcA, std::string srcB)
	{
		Signal *pSrcA=getSignalByName(srcA);
		if(pSrcA==NULL)
			throw std::string("MakeFPMultiplier - Attempt to use non-existant signal as source.");
		if(!pSrcA->isFP())
			throw std::string("MakeFPMultiplier - Attempt to use non-FP signal as source.");
		int wEA=pSrcA->wE();
		int wFA=pSrcA->wF();

		Signal *pSrcB=getSignalByName(srcB);
		if(pSrcB==NULL)
			throw std::string("MakeFPMultiplier - Attempt to use non-existant signal as source.");
		if(!pSrcB->isFP())
			throw std::string("MakeFPMultiplier - Attempt to use non-FP signal as source.");
		int wEB=pSrcB->wE();
		int wFB=pSrcB->wF();

		int wER=std::max(wEA,wEB);
		int wFR=std::max(wFA,wFB);
		
		std::map<std::string,double> delays;
		delays["X"]=getCriticalPath();
		delays["Y"]=getCriticalPath();
		
		FPMultiplier *pMul = new FPMultiplier(getTarget(), wEA, wFA, wEB, wFB, wER, wFR, true, true, 1, 1, delays);
		oplist.push_back(pMul); 
		
		inPortMap(pMul, "X", srcA);
		inPortMap(pMul, "Y", srcB);
		outPortMap(pMul, "R", name+"_R");
		vhdl<<instance(pMul, name+"_inst");
		syncCycleFromSignal(name+"_R", pMul->getOutputDelay("R") );	

		vhdl<<declareFP(name,wER,wFR)<<" <= "<<name<<"_R;\n";

		return pMul;
	}
	
	Operator *MakeFPConstMult(std::string name, std::string src, std::string constant)
	{
		Signal *pSrc=getSignalByName(src);
		if(pSrc==NULL)
			throw std::string("MakeFPNegate - Attempt to use non-existant signal as source.");
		if(!pSrc->isFP())
			throw std::string("MakeFPNegate - Attempt to use non-FP signal as source.");
		int wE=pSrc->wE();
		int wF=pSrc->wF();
		
		FPConstMult *pMul = new FPConstMult(getTarget(), wE, wF, wE, wF, wF, constant);
		oplist.push_back(pMul); 
		
		inPortMap(pMul, "X", src);
		outPortMap(pMul, "R", name+"_R");
		vhdl<<instance(pMul, name+"_inst");
		syncCycleFromSignal(name+"_R", pMul->getOutputDelay("R") );	

		vhdl<<declareFP(name,wE,wF)<<" <= "<<name<<"_R;\n";

		return pMul;
	}

    std::string EncodeFPConstant(std::string constant, int wE, int wF)
    {
		sollya_node_t node= parseString(constant.c_str());
		if (node == 0) {
		    ostringstream error;
		    error << srcFileName << ": Unable to parse string "<< constant << " as a numeric constant" <<endl;
		    throw error.str();
		}

		mpfr_t C;
		mpfr_init2(C, getToolPrecision());
		evaluateConstantExpression(C, node,  getToolPrecision());
		std::string sC=std::string("\"")+fp2bin(C, wE, wF)+"\"";
		mpfr_clear(C);

		return sC;
    }

    std::string EncodeFPConstant(double constant, int wE, int wF)
    {
	mpfr_t C;
	mpfr_init2(C, getToolPrecision());
	mpfr_set_d(C, constant, MPFR_RNDN);
	std::string sC=std::string("\"")+fp2bin(C, wE, wF)+"\"";
	mpfr_clear(C);
	return sC;
    }

	Operator *MakeFPInverse(std::string name, std::string src)
	{
		Signal *pSrc=getSignalByName(src);
		if(pSrc==NULL)
			throw std::string("MakeFPNegate - Attempt to use non-existant signal as source.");
		if(!pSrc->isFP())
			throw std::string("MakeFPNegate - Attempt to use non-FP signal as source.");
		int wE=pSrc->wE();
		int wF=pSrc->wF();

		
		FPDiv *pDiv = new FPDiv(getTarget(), wE, wF);
		oplist.push_back(pDiv); 
		
		inPortMapCst(pDiv, "X", EncodeFPConstant(1.0,wE,wF));
		inPortMap(pDiv, "Y", src);
		outPortMap(pDiv, "R", name+"_R");
		vhdl<<instance(pDiv, name+"_inst");
		syncCycleFromSignal(name+"_R", pDiv->getOutputDelay("R") );	

		vhdl<<declareFP(name,wE,wF)<<" <= "<<name<<"_R;\n";

		return pDiv;
	}

	Operator *MakeFPConstAdd(std::string name, std::string src, std::string constant)
	{
		Signal *pSrc=getSignalByName(src);
		if(pSrc==NULL)
			throw std::string("MakeFPConstAdd - Attempt to use non-existant signal as source.");
		if(!pSrc->isFP())
			throw std::string("MakeFPConstAdd - Attempt to use non-FP signal as source.");
		int wE=pSrc->wE();
		int wF=pSrc->wF();

		
		FPAdderSinglePath *pAdd = new FPAdderSinglePath(getTarget(), wE, wF, wE, wF, wE, wF,
						    inDelayMap("X", getTarget()->localWireDelay() + getCriticalPath()));
		oplist.push_back(pAdd); 
		
		inPortMap(pAdd, "X", src);
		inPortMapCst(pAdd, "Y", EncodeFPConstant(constant,wE,wF));
		outPortMap(pAdd, "R", name+"_R");
		vhdl<<instance(pAdd, name+"_inst");
		syncCycleFromSignal(name+"_R", pAdd->getOutputDelay("R") );	

		vhdl<<declareFP(name,wE,wF)<<" <= "<<name<<"_R;\n";

		return pAdd;
	}
	*/

	CompositeErf(
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
		name<<"Composite_erf_e"<<wDomE<<"_f"<<wDomF<<"_"<<getNewUId();
		setName(name.str());
		
		addFPInput("iX", m_wDomE, m_wDomF);
		addFPOutput("oY", m_wRanE, m_wRanF);
		
		assert(m_wDomE==m_wRanE);
		assert(m_wDomF==m_wRanF);

/* This is based on the approximation 7.1.26 from A&S:
p=0.3275911
a1=0.254829592
a2=-0.284496736
a3=1.421413741
a4=-1.453152027
a5=1.061405429

double t=1/(1+p*x);
erf = 1-(a1*t+a2*t^2+a3*t^3+a4*t^4+a5*t^5)*exp(-x^2)
*/

		setCriticalPath( getMaxInputDelays(inputDelays) );		

		MakeFPConstMult("px", "iX", "0.3275911");
		MakeFPConstAdd("px_plus_1", "px", "1");
		MakeFPInverse("t", "px_plus_1");

		// Now start doing the polynomial in t
		// t*(a1+t*(a2+t*(a3+t*(a4+t*a5))))
		MakeFPConstMult("mul5", "t", "1.061405429");    // t*a5
		MakeFPConstAdd("add4", "mul5", "-1.453152027"); // a4+(t*a5)
		MakeFPMultiplier("mul4", "t", "add4"); // t*(a4+(t*a5))
		MakeFPConstAdd("add3", "mul4", "1.421413741"); // a3+t*(a4+(t*a5))
		MakeFPMultiplier("mul3", "t", "add3"); // t*(a3+t*(a4+(t*a5)))
		MakeFPConstAdd("add2", "mul3", "-0.284496736"); // a2+t*(a3+t*(a4+(t*a5)))
		MakeFPMultiplier("mul2", "t", "add2"); // t*(a2+t*(a3+t*(a4+(t*a5))))
		MakeFPConstAdd("add1", "mul2", "0.254829592"); // a1+t*(a2+t*(a3+t*(a4+(t*a5))))
		MakeFPMultiplier("poly", "t", "add1"); // t*(a1+t*(a2+t*(a3+t*(a4+(t*a5)))))
		double polyDelay=getCriticalPath();

		// Ok, polynomial is done, now go back and do the exp(-x^2)

		setCycleFromSignal("iX");
		setCriticalPath( getMaxInputDelays(inputDelays) );		
		
		MakeFPNegate("neg_iX", "iX");
		MakeFPExp("exp_neg_iX", "neg_iX");
		double expDelay=getCriticalPath();

		syncCycleFromSignal("poly");// Join both back up
		setCriticalPath(std::max(polyDelay, expDelay));

		MakeFPMultiplier("poly_exp", "poly", "exp_neg_iX");
		MakeFPNegate("neg_poly_exp", "poly_exp");
		MakeFPConstAdd("res", "neg_poly_exp", "1.0");

		vhdl<<"oY <= res;\n";

	}

};

static void CompositeErfFactoryUsage(std::ostream &dst)
{
	OperatorFactory::classic_OP(dst, "CompositeErf", "wE wF", false);
	dst << "    Generates an implementation of erf following A&S 7.1.26\n";
	dst << "	      (wE,wF) - Floating-point format of input and output.\n";
	dst << "         NOTE : This is for testing purposes, and should never be used.\n";
}

static Operator *CompositeErfFactoryParser(Target *target ,const std::vector<std::string> &args,int &consumed)
{
	unsigned nargs = 2;
	if (args.size()<nargs)
		throw std::string("CompositeErfFactoryParser - Not enough arguments, check usage.");
	consumed += nargs;
	
	int wDomE = atoi(args[0].c_str());
	int wDomF = atoi(args[1].c_str());
	
	return new CompositeErf(target,      wDomE, wDomF);
}

void CompositeErf_registerFactory()
{
	DefaultOperatorFactory::Register(
		"CompositeErf",
		"private",
		CompositeErfFactoryUsage,
		CompositeErfFactoryParser
	);
}

}; // float_approx
}; // random
}; // flopoco
