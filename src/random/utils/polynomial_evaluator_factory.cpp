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

#include "FixedPointFunctions/PolynomialEvaluator.hpp"

using namespace std;

namespace flopoco
{

extern vector<Operator *> oplist;	
	
namespace random
{

static void PolynomialEvaluatorTableFactoryUsage(std::ostream &dst)
{
	OperatorFactory::classic_OP(dst, "PolynomialEvaluator", "degree targetPrec inputFormat coeffFormat+", false);
	dst << "    Create a polynomial evaluator with given input and coefficient formats\n";
	dst << "	      degree - Degree of polynomial\n";
	dst << "	      targetPrec - Binary weight of the faithful rounding, i.e. target error of 2^(targetPrec-1).\n";
	dst << "           inputFormat - Format for the polynomial argument (currently must be unsigned)\n";
	dst << "           coeffFormat - Format for each coefficient from a[0] up to a[degree] (currently must be signed).\n";
	dst << "      Formats are specified as \"U|S;msb;lsb\", where u or s says whether it is two's complement,\n";
	dst << "      msb gives the weight of the msb, and lsb gives the weight of the lsb. Some examples are:,\n";
	dst << "         \"U;15;0\"  : 16 bit unsigned integer, bits represent 2^15,2^14,...,2^1,2^0\n";
	dst << "         \"S;32;0\"  : 32 bit signed integer, bits represent -2^31,2^30,...,2^1,2^0\n";
	dst << "         \"U;-1;-16\"  : 16 bit unsigned number with 16 fractional bits in range [0,1)\n";
	dst << "         \"S;0;-15\"  : 16 bit signed number with 15 fractional bits in range [-1,1)\n";
	dst << "         \"S;7;-8\"  : 16 bit signed number with 8 fraction bits in range [-128,+128)\n";
	dst << "\n";
}



PolynomialEvaluator::format_t ParseFormat(const std::string &x)
{
	std::string left=x;
	if(left.size()<5)
		throw std::string("ParseFixedFormat('"+x+"') - Need at least five characters.");
	
	PolynomialEvaluator::format_t res;
	
	if(left[0]=='U'){
		res.isSigned=false;
	}else if(left[0]=='S'){
		res.isSigned=true;
	}else{
		throw std::string("ParseFixedFormat('"+x+"') - Format should start with 'U' or 'S'.");
	}
	
	if(left[1]!=';')
		throw std::string("ParseFixedFormat('"+x+"') - 'U' or 'S' must be followed by ';'");
	
	left=left.substr(2,-1);
	
	int split=left.find(';');
	if(split==-1)
		throw std::string("ParseFixedFormat('"+x+"') - Couldn't find second ';'");
	
	std::string msbStr=left.substr(0,split);
	std::string lsbStr=left.substr(split+1, -1);
	
	if(msbStr.size()==0)
		throw std::string("ParseFixedFormat('"+x+"') - Couldn't convert msb to a number");
	res.msb=boost::lexical_cast<int>(msbStr);	
	
	if(lsbStr.size()==0)
		throw std::string("ParseFixedFormat('"+x+"') - Couldn't convert lsb to a number");
	res.lsb=boost::lexical_cast<int>(lsbStr);	
	
	return res;
}

static Operator *PolynomialEvaluatorFactoryParser(Target *target ,const std::vector<std::string> &args,int &consumed)
{
	typedef PolynomialEvaluator::format_t format_t;
	
	if(args.size()<4)
		throw std::string("PolynomialEvaluatorFactory - Not enough arguments, check usage.");
	
	int degree=boost::lexical_cast<int>(args[0]);
	if(args.size()<4+degree)
		throw std::string("PolynomialEvaluatorFactory - Not enough arguments, check usage.");
	
	int targetPrec=boost::lexical_cast<int>(args[1]);
	format_t inputFormat=ParseFormat(args[2]);
	std::vector<format_t> coeffFormats;
	for(int i=0;i<=degree;i++){
		coeffFormats.push_back(ParseFormat(args.at(3+i)));
	}
	consumed=4+degree;
	
	return PolynomialEvaluator::Create(Target *target,
		coeffFormats,	//! Format for each of the coefficients, with coeffFormats[0]=a0, etc.
		format_t inputFormat,	//! The polynomial input format
		targetPrec,				//! The LSB to which the evaluator should be accurate
		0.0					//! The maximum approximation error which already occurred. Must have approxError < 2^(outputLsb-1)
	);
}

void PolynomialEvaluator_registerFactory()
{
	DefaultOperatorFactory::Register(
		"PolynomialEvaluator",
		"operator",
		PolynomialEvaluatorFactoryUsage,
		PolynomialEvaluatorFactoryParser,
		DefaultOperatorFactory::ParameterList(
			DefaultOperatorFactory::Parameters("1", "-8", "U;-1;-8", "S;-1;-8", "S;-1;-8")
		)
	);
}

};
};


