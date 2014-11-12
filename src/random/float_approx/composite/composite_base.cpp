#include "composite_base.hpp"

#include "utils.hpp"
#include "sollya.h"

#include "FPMultiplier.hpp"
#include "FPSquarer.hpp"
#include "FPExp.hpp"
#include "FPDiv.hpp"
#include "ConstMult/FPConstMult.hpp"
#include "FPSquarer.hpp"
#include "FPAdderSinglePath.hpp"

#include <cassert>

namespace flopoco
{
namespace random
{
namespace float_approx
{

	void CompositeBase::MakeFPNegate(std::string dst, std::string src)
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
	
	Operator *CompositeBase::MakeFPExp(std::string name, std::string src)
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

	Operator *CompositeBase::MakeFPSquarer(std::string name, std::string src)
	{
		Signal *pSrc=getSignalByName(src);
		if(pSrc==NULL)
			throw std::string("MakeFPSquarer - Attempt to use non-existant signal as source.");
		if(!pSrc->isFP())
			throw std::string("MakeFPSquarer - Attempt to use non-FP signal as source.");
		int wE=pSrc->wE();
		int wF=pSrc->wF();
		
		FPSquarer *pSquarer = new FPSquarer(getTarget(), wE, wF, wF);
		oplist.push_back(pSquarer); 
		
		inPortMap(pSquarer, "X", src);
		outPortMap(pSquarer, "R", name+"_R");
		vhdl<<instance(pSquarer, name+"_inst");
		syncCycleFromSignal(name+"_R", pSquarer->getOutputDelay("R") );	

		vhdl<<declareFP(name,wE,wF)<<" <= "<<name<<"_R;\n";

		return pSquarer;
	}

	Operator *CompositeBase::MakeFPMultiplier(std::string name, std::string srcA, std::string srcB)
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
	
	Operator *CompositeBase::MakeFPConstMult(std::string name, std::string src, std::string constant)
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

    std::string CompositeBase::EncodeFPConstant(std::string constant, int wE, int wF)
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

    std::string CompositeBase::EncodeFPConstant(double constant, int wE, int wF)
    {
	mpfr_t C;
	mpfr_init2(C, getToolPrecision());
	mpfr_set_d(C, constant, MPFR_RNDN);
	std::string sC=std::string("\"")+fp2bin(C, wE, wF)+"\"";
	mpfr_clear(C);
	return sC;
    }

	Operator *CompositeBase::MakeFPInverse(std::string name, std::string src)
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

	Operator *CompositeBase::MakeFPConstAdd(std::string name, std::string src, std::string constant)
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


	CompositeBase::CompositeBase(
			Target *target
		)
		: Operator(target)
	{}

}; // float_approx
}; // random
}; // flopoco
