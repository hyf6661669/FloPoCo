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
#include "TableTransform.hpp"
#include "random/utils/operator_factory.hpp"

// For twos complement stuff
#include "CLTTransform.hpp"

#include "FixedPointFunctions/Function.hpp"
#include "Table.hpp"
#include "UtilSollya.hh"

#include "random/distributions/table_distribution.hpp"
#include "random/distributions/gaussian_distribution.hpp"
#include "random/moment_correction/correct_distribution.hpp"

#include "random/utils/fft/convolve_mpreal.hpp"

using namespace std;

namespace flopoco
{

extern vector<Operator *> oplist;	
	
namespace random
{

class CLTCorrectedTransform :
	public RngTransformOperator
{
	Distribution<mpfr::mpreal>::TypePtr m_distribution;
	
	CLTCorrectedTransform(Target *target)
		: RngTransformOperator(target)
		, IRngTransformDistributions
	{};
	
	virtual void outputVHDL(std::ostream& o, std::string name)
	{ throw std::string("CLTCorrectedTransform::outputVHDL - Not implemented (distribution only at the moment)."); }

	virtual unsigned uniformInputBits() const
	{ return 128; }
	
	virtual std::string uniformInputName() const
	{ return "urng"; }

	virtual unsigned nonUniformOutputCount() const
	{ return 1; }

	virtual bool nonUniformOutputsAreHomogenous() const
	{ return true; }

	virtual unsigned nonUniformOutputWidth(int i) const
	{ return 16; }

	//! The name of the output port associated with this signal
	virtual std::string nonUniformOutputName(int i) const
	{ return "X"; }
	
	virtual typename Distribution<mpfr::mpreal>::TypePtr nonUniformOutputDistribution(int i, unsigned prec) const
	{
		if(!m_distribution){
			// First, build up the CLT distribution. The stated technique is to start
			// from 16 bit samples, then to add, but every time they add one LSB is
			// dropped. Worryingly, this seems to imply quite a lot of bias, as the
			// thing won't be zero mean to start with, and is going to get dragged
			// down with every lost LSB.
			
			mpfr::mpreal one(1.0, prec);
			
			std::vector<mpfr::mpreal> ones(1<<16, ldexp(one, -16));
			
			// First distribution in Q(16,15)
			HistogramDistribution<mpfr::mpreal>::TypePtr curr=boost::make_shared<HistogramDistribution<mpfr::mrpeal> >(-one, ldexp(one, -15), ones);
			
			// First adder and truncation
			curr=SelfAddDistribution(curr, 2);
			curr=
		}
		return m_distribution;
	}
};	



static Operator *GRNGTableFactoryParser(Target *target ,const std::vector<std::string> &args,int &consumed)
{
	unsigned nargs = 5;
	if (args.size()<nargs)
		throw std::string("GRNGTableFactory - Not enough arguments, check usage.");
	consumed += nargs;
	
	int k = atoi(args[0].c_str());
	int wF = atoi(args[1].c_str());
	double stddev=parseSollyaConstant(args[2]);
	std::string correction=args[3];
	std::string quantisation=args[4];
	
	return MakeGRNGTable(target, k, wF, stddev, correction, quantisation);
}

void GRNGTableTransform_registerFactory()
{
	DefaultOperatorFactory::Register(
		"grng_table_transform",
		"operator;rng_transform",
		GRNGTableFactoryUsage,
		GRNGTableFactoryParser,
		DefaultOperatorFactory::ParameterList(
			DefaultOperatorFactory::Parameters("6", "12", "1.0", "auto", "auto")
		)
	);
}

};
};


