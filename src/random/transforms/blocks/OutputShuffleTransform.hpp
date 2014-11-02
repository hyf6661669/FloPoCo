/* Each Operator declared within the flopoco framework has 
to inherit the class Operator and overload some functions listed below*/
#ifndef output_shuffle_transform_hpp
#define output_shuffle_transform_hpp
#include <vector>
#include <sstream>
#include <string>
#include <fstream>
#include "gmp.h"
#include "mpfr.h"
#include <math.h>
#include <assert.h>
#include <iomanip>
#include <time.h>

#include "../RngTransformOperator.hpp"

/* This file contains a lot of useful functions to manipulate vhdl */
#include "utils.hpp"

#include "random/distributions/table_distribution.hpp"

/*  All flopoco operators and utility functions are declared within
the flopoco namespace.
*/
namespace flopoco{
namespace random{

/* Given a base transform, applies a random shuffle to each output stream */
class OutputShuffleTransform
	: public RngTransformOperator
	// Need to think about what it means to support this
	//, public IRngTransformDistributions
	{
private:
	unsigned m_k;
	RngTransformOperator *m_base;
	
public:
	/*! 
		\param k Table size is 2^k, and will advance by 2^(k-1) bits each time
	*/
	OutputShuffleTransform(Target* target, unsigned k, RngTransformOperator *base);
	~OutputShuffleTransform();

	void emulate(TestCase * tc);

	void buildRandomTestCases(TestCaseList* tcl, int n);

	TestCase* buildRandomTestCase(int i);

	virtual unsigned uniformInputBits() const
	{ return m_base->uniformInputBits() + m_base->nonUniformOutputCount()*(m_k-1); }
	
	virtual std::string uniformInputName() const
	{ return m_base->uniformInputName(); }

	virtual unsigned nonUniformOutputCount() const
	{ return m_base->nonUniformOutputCount(); }
	
	virtual bool nonUniformOutputsAreHomogenous() const
	{ return m_base->nonUniformOutputsAreHomogenous(); }
	
	virtual unsigned nonUniformOutputWidth(int i) const
	{ return m_base->nonUniformOutputWidth(i); }
	
	virtual std::string nonUniformOutputName(int i) const
	{ return m_base->nonUniformOutputName(i); }
	
	static void registerFactory();
};


}; // random
}; // flopoco
#endif
