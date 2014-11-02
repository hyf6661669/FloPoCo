// general c++ library for manipulating streams
#include <iostream>
#include <sstream>
#include <math.h>	// for NaN

#include <random/utils/mpreal/boost_math_mpreal.hpp>



/* header of libraries to manipulate multiprecision numbers
  There will be used in the emulate function to manipulate arbitraly large
  entries */
#include "gmp.h"
#include "mpfr.h"
#include "FPNumber.hpp"

// include the header of the Operator
#include "OutputShuffle.hpp"
#include "OutputShuffleTransform.hpp"
#include "random/utils/operator_factory.hpp"


using namespace std;

namespace flopoco
{

extern vector<Operator *> oplist;	
	
namespace random
{

OutputShuffleTransform::OutputShuffleTransform(Target* target, unsigned k, RngTransformOperator *base)
	: RngTransformOperator(target)
	, m_k(k)
	, m_base(base)
{
	REPORT(DETAILED, "  OutputShuffleTransform() begin");	
	
	std::stringstream acc;
	acc<<"OutputShuffleTransform_"<<m_k<<"_"<<base->getName()<<"_uid"<<getNewUId();
	setName(acc.str());
	
	REPORT(DETAILED, "    Adding inputs.");
	
	addInput(uniformInputName(), uniformInputBits());
	vhdl<<tab<<declare("base_uniform", base->uniformInputBits())<<" <= "<<uniformInputName()<<"("<<base->uniformInputBits()-1<<" downto 0);\n";
	unsigned inOffset=base->uniformInputBits();
	for(unsigned i=0;i<nonUniformOutputCount();i++){
		addOutput(nonUniformOutputName(i), nonUniformOutputWidth(i));
		vhdl<<tab<<declare(join("shuffle_uniform_",i),m_k-1)<<" <= "<<uniformInputName()<<"("<<inOffset+m_k-2<<" downto "<<inOffset<<");\n";
		inOffset+=m_k-1;
	}	
	
	inPortMap(base, base->uniformInputName(), "base_uniform");
	for(unsigned i=0;i<nonUniformOutputCount();i++){
		outPortMap(base, base->nonUniformOutputName(i), join("base_output_",i));
	}
	
	vhdl<<tab<<instance(base, "transformBase");
	
	for(unsigned i=0;i<nonUniformOutputCount();i++){
		OutputShuffle *shuffle=new OutputShuffle(target, m_k, nonUniformOutputWidth(i));
		oplist.push_back(shuffle);
		
		// TODO : Note,there is a wasted opportunity to save regsters and SRLs here, as we
		// can short-cut the indices straight in. Don't want to mess around too much with
		// cycles right now though.
		
		setCycleFromSignal(join("base_output_",i));
		inPortMap(shuffle, "iData", join("base_output_",i));
		inPortMap(shuffle, "iIncrement", join("shuffle_uniform_",i));
		outPortMap(shuffle, "oData", join("shuffle_output_",i));
		vhdl<<instance(shuffle, join("shuffle_",i));
	}
	
	for(unsigned i=0;i<nonUniformOutputCount();i++){
		syncCycleFromSignal(join("shuffle_output_",i));
		vhdl<<tab<<nonUniformOutputName(i) << " <= "<<join("shuffle_output_",i)<<";\n";
	}
	
	REPORT(DETAILED, "  OutputShuffleTransform() complete.");
}

OutputShuffleTransform::~OutputShuffleTransform()
{}


void OutputShuffleTransform::emulate(TestCase * tc)
{
	throw std::string("OutputShuffle::emulate - Not implemented.");
}
	
TestCase* OutputShuffleTransform::buildRandomTestCase(int i)
{
	throw std::string("OutputShuffleTransform::buildRandomTestCase - Not implemented.");
}


static void OutputShuffleTransformFactoryUsage(std::ostream &dst)
{
	OperatorFactory::classic_OP(dst, "OutputShuffleTransform", "k", false);
	dst << "    Applies a shuffle according to the preceeding transform operator\n";
	dst << "	      k - Number of table address bits, iIndex will have k-1 bits.\n";
}

static Operator *OutputShuffleTransformFactoryParser(Target *target ,const std::vector<std::string> &args,int &consumed)
{
	unsigned nargs = 1;
	if (args.size()<nargs)
		throw std::string("OutputShuffleTransformFactory - Not enough arguments, check usage.");
	consumed += nargs;
	
	int k = atoi(args[0].c_str());
	
	vector<Operator*> *ops=target->getGlobalOpListRef();
	if(ops->size()==0)
		throw std::string("OutputShuffleTransformParser - No previous operator to work with.");
	
	Operator *op=ops->back();
	
	RngTransformOperator *base=dynamic_cast<RngTransformOperator*>(op);
	if(base==NULL)
		throw std::string("OutputShuffleTransformParser - Previous operator is not an RngTransformOperator.");
	
	
	if(::flopoco::verbose>=DETAILED)
		std::cerr<<"OutputShuffleParser : k="<<k<<", base="<<base->getName()<<"\n";

	if(k<1)
		throw std::string("OutputShuffleTransformParser - k must be at least 3.");
	
	return new OutputShuffleTransform(target, k ,base);
}

void OutputShuffleTransform::registerFactory()
{
	DefaultOperatorFactory::Register(
		"OutputShuffleTransform",
		"operator;rng_transform",
		flopoco::random::OutputShuffleTransformFactoryUsage,
		flopoco::random::OutputShuffleTransformFactoryParser
	);
}

};
};


