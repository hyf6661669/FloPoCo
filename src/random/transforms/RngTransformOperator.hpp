#ifndef random_transform_rng_transform_operator_hpp
#define random_transform_rng_transform_operator_hpp

#include "Operator.hpp"

#include "random/distributions/distribution.hpp"

#include "mpreal.h"
//#include "random/utils/mpreal/boost_math_mpreal.hpp"

namespace flopoco
{
namespace random
{

class RngTransformOperator
	: public flopoco::Operator
{
protected:
	RngTransformOperator(Target* target, map<string, double> inputDelays = emptyDelayMap)
		: Operator(target, inputDelays)
	{}
public:
	virtual ~RngTransformOperator()
	{}

	virtual unsigned uniformInputBits() const=0;
	virtual std::string uniformInputName() const =0;

	virtual unsigned nonUniformOutputCount() const=0;

	//! Return true if all outputs have the same types and dstributions
	virtual bool nonUniformOutputsAreHomogenous() const=0;

	//! The width of the output port (no matter what the representation is)
	virtual unsigned nonUniformOutputWidth(int i) const=0;

	//! The name of the output port associated with this signal
	virtual std::string nonUniformOutputName(int i) const=0;
		
	virtual void Simulate(
		unsigned nSamples,
		const mpz_class *pUniformInputs,	// Array of nSamples uniform input bits, with at least uniformInputBits() content
		mpz_class *pNonUniformOutputs		// Array of nSamples*nonUniformOutputCount() gmp instances
	) const {
		// There should be a default version done in terms of emulate, but I'm out of time
		throw std::runtime_error("Simulate not implemented on this RngTransformOperator.");
	}
};

/*! Used for transforms which are able to precisely specify their output distribution. Not all
	transforms can do this.
*/
class IRngTransformDistributions
{
public:
	virtual ~IRngTransformDistributions()
	{}
		
	//! Return the distribution of the output, with precision of prec or higher
	virtual Distribution<mpfr::mpreal>::TypePtr nonUniformOutputDistribution(int i, unsigned prec) const=0;
};

};
};

#endif
