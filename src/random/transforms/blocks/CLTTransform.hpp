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
class CLTTransform : public RngTransformOperator {
private:
	unsigned m_baseWidth;
public:
	CLTTransform(Target* target, int baseW);
	~CLTTransform();

	void emulate(TestCase * tc);

	void buildRandomTestCases(TestCaseList* tcl, int n);

	TestCase* buildRandomTestCase(int i);

	virtual unsigned uniformInputBits() const
	{ return m_baseWidth*2; }
	
	virtual std::string uniformInputName() const
	{ return "iUniformBits"; }

	virtual unsigned nonUniformOutputCount() const
	{ return 1; }
	
	virtual unsigned nonUniformOutputWidth(int) const
	{ return m_baseWidth+1; }
	
	virtual std::string nonUniformOutputName(int i) const
	{ return "oRng"; }
	
	static void registerFactory();
};

mpz_class toTwosComplement(const mpz_class &x, unsigned bits);
mpz_class fromTwosComplement(const mpz_class &x, unsigned bits);

}; // random
}; // flopoco
#endif
