#include <iostream>
#include <sstream>

#include "gmp.h"
#include "mpfr.h"

#include "IntPower.hpp"
#include "FormalBinaryProduct.hpp"

using namespace flopoco;

IntPower::IntPower(Target* target,
                   size_t wIn, size_t n,
		   std::map<std::string,double> inputDelays)
	:GenericBinaryPolynomial(target,
			ProductIR::identity(wIn).toPow(n).simplifyInPlace(),
			inputDelays), wIn(wIn), n(n) {
};

	
void IntPower::emulate(TestCase * tc) {
	mpz_class x = tc->getInputValue ("X"), res(1);
	for (size_t i = 0; i < n; i++) {
		res *= x;
	}
	res &= ((mpz_class(1) << p.data.size()) - 1);
	tc->addExpectedOutput ("R", res);
}

void IntPower::buildStandardTestCases(TestCaseList * tcl) {
}

void IntPower::buildRandomTestCases(TestCaseList *  tcl, int n) {
}

TestCase* IntPower::buildRandomTestCases(int i) {
  TestCase* tc = new TestCase(this);
  return tc;
}

OperatorPtr IntPower::parseArguments( Target *target, vector<string> &args ) {
    int wIn, n;
    UserInterface::parseInt( args, "wIn", &wIn );
    UserInterface::parseInt( args, "n", &n );

    return new IntPower( target, wIn, n);
}

void IntPower::registerFactory() {
    UserInterface::add( "IntPower", // name
                        "Integer power of n, unsigned, with precision wIn (NPY)", // description, string
                        "BasicInteger", // category, from the list defined in UserInterface.cpp
                        "", //seeAlso
                        // where parameterDescription is parameterName (parameterType)[=defaultValue]: parameterDescriptionString
                        "wIn (int): Input Wordsize; \
                        n (int): power.",
                        "",
                        IntPower::parseArguments
                      ) ;
}
