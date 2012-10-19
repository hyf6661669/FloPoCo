/* Each Operator declared within the flopoco framework has 
to inherit the class Operator and overload some functions listed below*/
#ifndef random_clt_transform_hpp
#define random_clt_transform_hpp
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



/*  All flopoco operators and utility functions are declared within
the flopoco namespace.
*/
namespace flopoco{
namespace random{

/* Given an operator, it will create an array of them and apply a hadamard transform */
class TableTransform
	: public RngTransformOperator
	{
private:
	unsigned m_wElts;
	bool m_addRandomSign;
	std::vector<mpz_class> m_elements;
	unsigned m_log2n;
public:
	/*! Construct a table that will be uniformly sampled.
		\param wElts How wide each element of the table is
		\param elements A set of elements, which must be a binary power size
		\param addRandomSize If true, then elements are treated as unsigned numbers, and output will be signed(eElts+1) with a random sign
	*/
	TableTransform(Target* target, int wElts, const std::vector<mpz_class> &elements, bool addRandomSign);
	~TableTransform();

	void emulate(TestCase * tc);

	void buildRandomTestCases(TestCaseList* tcl, int n);

	TestCase* buildRandomTestCase(int i);

	virtual unsigned uniformInputBits() const
	{ return m_log2n + (m_addRandomSign?1:0); }
	
	virtual std::string uniformInputName() const
	{ return "iUniformBits"; }

	virtual unsigned nonUniformOutputCount() const
	{ return 1; }
	
	virtual unsigned nonUniformOutputWidth(int) const
	{ return m_wElts+(m_addRandomSign?1:0); }
	
	virtual std::string nonUniformOutputName(int i) const
	{ return "oRng"; }
};


}; // random
}; // flopoco
#endif
