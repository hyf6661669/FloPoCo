#include <iostream>
#include <sstream>
#include <math.h>

#include "RngTransformOperator.hpp"

#include "random/utils/operator_factory.hpp"

using namespace std;

namespace flopoco
{

extern vector<Operator *> oplist;	
	
namespace random
{
	
template<class T>
void DumpStats(std::ostream &dest_file, moptl::DiscreteDistribution<T>::TypePtr got, moptl::DiscreteDistribution<T>::TypePtr target, T startx)
{
	Accumulator<T> cdfCdf, targetCdf;
	
	
}

static void CLTFactoryUsage(std::ostream &dst)
{
	OperatorFactory::classic_OP(dst, "TransformStats", "dest_file", false);
	dst << "       Calculates \n";
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
