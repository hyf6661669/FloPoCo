// general c++ library for manipulating streams
#include <iostream>
#include <sstream>
#include <math.h>	// for NaN



/* header of libraries to manipulate multiprecision numbers
  There will be used in the emulate function to manipulate arbitraly large
  entries */
#include "gmp.h"
#include "mpfr.h"
#include "FPNumber.hpp"

// include the header of the Operator
#include "CLTTransform.hpp"

#include "random/utils/operator_factory.hpp"

using namespace std;

namespace flopoco
{

extern vector<Operator *> oplist;	
	
namespace random
{
	
CLTTransform::CLTTransform(Target* target, int baseWidth)
	: RngTransformOperator(target)
	, m_baseWidth(baseWidth)
{
	std::stringstream acc;
	acc<<"CLTTransform_w"<<baseWidth<<"_uid"<<getNewUId();
	setName(acc.str());
	
	addInput(uniformInputName(), uniformInputBits());
	addOutput(nonUniformOutputName(0), nonUniformOutputWidth(0));
	
	vhdl << declare("left",baseWidth) << "<= iUniformBits"<<range(baseWidth*2-1,baseWidth)<<";\n";
	vhdl << declare("right",baseWidth) << "<= iUniformBits"<<range(baseWidth-1,0)<<";\n";
	
	vhdl<<declare("res",baseWidth+1) << " <= "<<zeroExtend("left",baseWidth+1) << " - "<<zeroExtend("right",baseWidth+1)<<";\n";
	nextCycle();
	
	vhdl<<nonUniformOutputName(0)<<" <= res;\n";
	
}

CLTTransform::~CLTTransform()
{}
	
mpz_class toTwosComplement(const mpz_class &x, unsigned bits)
{
	if(x>=0)
		return x;
	
	// TODO : Lamest conversion ever!
	// mpz bitwise stuff is in twos-complement already
	mpz_class res;
	for(int i=bits-1;i>=0;i--){
		res=(res<<1) + (mpz_tstbit(x.get_mpz_t(), i));
	}
	return res;
}

mpz_class fromTwosComplement(const mpz_class &x, unsigned bits)
{
	if(mpz_tstbit(x.get_mpz_t(),bits-1)==0)
		return x;
	
	// TODO : second lamest conversion ever!
	mpz_class res=x;
	mpz_tdiv_r_2exp(res.get_mpz_t(), res.get_mpz_t(), bits-1);	// Chop off the positive part
	return res + (mpz_class(-1)<<(bits-1));	// and add the sign...
}

void CLTTransform::emulate(TestCase * tc)
{
	mpz_class bits=tc->getInputValue(uniformInputName());
	
	mpz_class left=bits>>m_baseWidth;
	mpz_class right=bits;
	mpz_tdiv_r_2exp(right.get_mpz_t(), right.get_mpz_t(), m_baseWidth);
	mpz_class res=left-right;
	
	mpz_class enc=toTwosComplement(res, m_baseWidth+1);
	assert(enc>=0);
	tc->addExpectedOutput(nonUniformOutputName(0), enc);
	
	assert(res==fromTwosComplement(enc, m_baseWidth+1));
}

	
TestCase* CLTTransform::buildRandomTestCase(int i)
{
	TestCase *tc=new TestCase(this);
	
	tc->addInput(uniformInputName(), getLargeRandom(uniformInputBits()));
	emulate(tc);
	
  	return tc;
}

static void CLTFactoryUsage(std::ostream &dst)
{
	OperatorFactory::classic_OP(dst, "clt_transform", "baseWidth k", false);
	dst << "       Generates a CLT transform with k uniform inputs\n";
	dst << "	baseWidth - How many bits per base uniform generator.\n";
	dst << "	k - Number of input uniforms (must be even, and greater than zero)\n";
}

static Operator *CLTFactoryParser(Target *target ,const std::vector<std::string> &args,int &consumed)
{
	unsigned nargs = 2;
	if (args.size()<nargs)
		throw std::string("Not enough arguments.");
	
	if(::flopoco::verbose >= DETAILED)
		std::cerr<<"CLTFactoryParser, wBase="<<args[0]<<", k="<<args[1]<<"\n";
	
	int wBase = atoi(args[0].c_str());
	int k = atoi(args[1].c_str());
	
	if(k!=2)
		throw std::string("clt_transform - Only k==2 is supported at the moment.");
	
	consumed = nargs;
	return new flopoco::random::CLTTransform(target, wBase);
}

void CLTTransform::registerFactory()
{
	DefaultOperatorFactory::Register(
		"clt_transform",
		"operator;rng_transform",
		CLTFactoryUsage,
		CLTFactoryParser,
		DefaultOperatorFactory::ParameterList(
			DefaultOperatorFactory::Parameters("8", "2"),
			DefaultOperatorFactory::Parameters("16", "2")
		)
	);
}

};
};
