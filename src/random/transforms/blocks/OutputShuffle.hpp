/* Each Operator declared within the flopoco framework has 
to inherit the class Operator and overload some functions listed below*/
#ifndef output_shuffle_hpp
#define output_shuffle_hpp
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

/* Randomly delay outputs */
class OutputShuffle
	: public Operator
	{
private:
	unsigned m_k, m_w;
public:
	/*! 
		\param k Table size is 2^k, and will advance by 2^(k-1) bits each time
	*/
	OutputShuffle(Target* target, unsigned k, unsigned w);
	~OutputShuffle();

	void emulate(TestCase * tc);

	void buildRandomTestCases(TestCaseList* tcl, int n);

	TestCase* buildRandomTestCase(int i);

	static void registerFactory();

	void outputVHDL(std::ostream&, std::string);
};


}; // random
}; // flopoco
#endif
