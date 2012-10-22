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

using namespace std;

namespace flopoco
{

extern vector<Operator *> oplist;	
	
namespace random
{

std::vector<double> CorrectTable(std::string correction, const std::vector<double> &elements)
{
	// These are empirical limits on how many elements are needed before the solver becomes stable
	// using double precision. For higher-precision solvers more might be possible, but the value is debatable.
	static const unsigned stability={
		1,					// linear is always fine
		8,					// cubic is very stable, and goes way down, but the curve is hugely distorted
		64,				// quintic is very stable, and works down to 32, but the distribution looks very distorted
		(1<<16),		// heptic, stable, but not convincing
		(1<<17)		// nonic, seems more stable than heptic (actually works at 2^14), but is very worrying
	}
	
	if(correction=="auto"){
		// These are empirically derived, and by no means optimal
		if(elements.size()>=2*stability[4]){
			correction="poly9";
		}else if(elements.size()>=2*stability[3]){
			correction="poly7";
		}else if(elements.size()>=2*stability[2]){
			correction="poly5";
		}else if(elements.size()>=2*stability[1]){
			correction="poly3";
		}else{
			correction="poly1";
		}
		if(::flopoco::verbose>=INFO)
			std::cerr<<"CorrectTable : chose correction="<<correction<<" for table with "<<elements.size()<<" elements.";
	}
	
	if(correction.substr(0,4)=="poly"){
		std::string index=correction.substr(4,-1);
		int index=atoi(index.c_str());
		
		#error "Here"
		if(::flopoco::verbose>=DETAILED)
			std::cerr<<"CorrectTable : got correction="<<correction<<" for table with "<<elements.size()<<" elements.";
	}else{
		throw std::string("CorrectTable : Didn't understand correction method '"+correction+"'.");
	}
}

std::vector<mpz_class> BuildTable(
	int k, int w, double stddev, std::string correction, std::string fixation
){
	ContinuousDistribution<double>::TypePtr norm=boost::make_shared<GaussianDistribution<double> >(0, stddev);
	
	int n=1<<k;
	std::vector<double> x(n);
	for(int i=0;i<n/2;i++){
		double u=(i+0.5)/n;
		x[i]=target->InvCdf(u);
		x[n-i-1]=-x[i];
	}
	
	TableDistribution<double>::TypePtr current=boost::make_shared<TableDistribution<double> >(&x[0], &x[n]);
	
	
}
	

static void GRNGTableFactoryUsage(std::ostream &dst)
{
	OperatorFactory::classic_OP(dst, "grng_table_transform", "k w stddev correctMethod fixMethod", false);
	dst << "    Generates a table transform which approximates the Gaussian distribution\n";
	dst << "	      k - Number of input address bits, table will have 2^k elements.\n";
	dst << "	      wF - Fractional bits of each table element\n";
	dst << "          stddev - Target standard deviation for the table (sollya expression)\n";
	dst << "          correctMethod - How to correct the table in the continuous domain.\n";
	dst << "            auto - Select some sensible correction method.\n";
	dst << "            poly1,poly3,...,poly9 - Apply odd polynomial correction of given degree.\n";
	dst << "          fixMethod - How to convert the continous table to fixed.\n";
	dst << "            auto - Select some sensible fixation method.\n";
	dst << "            round - Simply round the continuous entries.\n";
	dst << "      The transform will take a (k+1) bit uniform input, and produce a (w+1) bit signed output.\n";
}



static Operator *GRNGTableFactoryParser(Target *target ,const std::vector<std::string> &args,int &consumed)
{
	unsigned nargs = 4;
	if (args.size()<nargs)
		throw std::string("TableFactory - Not enough arguments, check usage.");
	consumed += nargs;
	
	int k = atoi(args[0].c_str());
	int w = atoi(args[1].c_str());
	std::string funcStr = args[2];
	bool addSign=atoi(args[3].c_str())!=0;

	if(k<1)
		throw std::string("TableFactory - k must be a positive integer.");
	if(w<1)
		throw std::string("TableFactory - w must be a positive integer.");
	if(w>=48)
		throw std::string("TableFactory - w must be less than 48 (currently table is built in double-precision).");
	
	if(DETAILED<=::flopoco::verbose)
		std::cerr<<"  Parsing sollya string \""<<funcStr<<"\" ... ";
	Function func(funcStr);
	if(DETAILED<=::flopoco::verbose)
		std::cerr<<"done\n";
	
	std::vector<mpz_class> contents(1<<k);
	for(int i=0;i<(1<<k);i++){
		double x=(i+0.5)/(1<<k);

		double y=func.eval(x);
		double ry=round(ldexp(y,w));
		
		if(DEBUG<=::flopoco::verbose)
			std::cerr<<"    "<<i<<", x="<<x<<", y="<<y<<", ry="<<ry<<"\n";
		
		if((ry<0.0) || (ry>=ldexp(1.0,w))){
			std::stringstream acc;
			acc<<"TableFactoryParser : For index "<<i<<", the value "<<y<<"=f("<<x<<") leads to out of range value "<<ry;
			throw std::string(acc.str());
		}
		contents[i]=ry;
	}
	
	return new TableTransform(target, w, contents, addSign);
}

void GRNGTableTransform_registerFactory()
{
	DefaultOperatorFactory::Register(
		"grng_table_transform",
		"operator;rng_transform",
		flopoco::random::TableFactoryUsage,
		flopoco::random::TableFactoryParser,
		DefaultOperatorFactory::ParameterList(
			DefaultOperatorFactory::Parameters("8", "8", "sin(x)", "0"),
			DefaultOperatorFactory::Parameters("8", "8", "sin(x)", "1")
		)
	);
}

};
};


