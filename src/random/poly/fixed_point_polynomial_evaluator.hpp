#ifndef flopoco_random_poly_fixed_point_polynomial_evaluator_hpp
#define flopoco_random_poly_fixed_point_polynomial_evaluator_hpp

#include "Operator.hpp"
#include "Target.hpp"

#include "mpreal.h"

#include "fixed_format_t.hpp"

namespace flopoco
{
namespace random
{

class FixedPointPolynomialEvaluator
	: public Operator
{
public:

protected:
	FixedPointPolynomialEvaluator(Target *target);

	fixed_format_t MultiplyType(const fixed_format_t &aType, const fixed_format_t &bType) const;
	fixed_format_t AddType(const fixed_format_t &a, const fixed_format_t &b) const;

	std::string ExtendExpr(const fixed_format_t &outType, const std::string srcExpr, const fixed_format_t &srcType) const;
	
	// If the src type has more lsbs than the output type, then round them off. The msbs should be the same
	std::string RoundExpr(const fixed_format_t &outType, const std::string srcExpr, const fixed_format_t &srcType);

	// Not const, as they modify vhdl
	fixed_format_t MultiplyStatement(std::string resName, std::string aName, const fixed_format_t &aType, std::string bName, const fixed_format_t &bType);
	fixed_format_t AddStatement(std::string resName, std::string aName, const fixed_format_t &aType, std::string bName, const fixed_format_t &bType);
public:
	virtual int getPolynomialDegree() const=0;
	virtual fixed_format_t getInputFormat() const=0;
	virtual fixed_format_t getCoefficientFormat(unsigned i) const=0;
	
	virtual fixed_format_t getOutputFormat() const=0;

	void emulate(TestCase *tc);	
	void buildStandardTestCases(TestCaseList* tcl);
};

FixedPointPolynomialEvaluator *CreateFixedPointPolynomialEvaluator(
	const fixed_format_t &inputFormat,
	const std::vector<fixed_format_t> &coefficientFormats,
	int outputLsb,
	const mpfr::mpreal &errorBudget,
	Target *target
);

FixedPointPolynomialEvaluator *CreateFixedPointPolynomialEvaluator(
	std::string selector,
	const fixed_format_t &inputFormat,
	const std::vector<fixed_format_t> &coefficientFormats,
	int outputLsb,
	const mpfr::mpreal &errorBudget,
	Target *target
);

/*
struct polynomial_spec_t
{
	mpfr::mpreal minX, maxX;
	std::vector<mpfr::mpreal> coefficients;
};

FixedPointPolynomialEvaluator *CreateFixedPointPolynomialEvaluator(
	const fixed_format_t &inputFormat,
	const std::vector<polynomial_spec_t> &polynomials,
	int outputLsb,
	const mpfr::mpreal &errorBudget,
	Target *target
);
*/

}; // random
}; // flopoco

#endif
