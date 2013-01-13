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

template<class T>
typename TableDistribution<T>::TypePtr CorrectTable(
	std::string correction,
	typename ContinuousDistribution<T>::TypePtr targetDist,
	typename TableDistribution<T>::TypePtr current
){
	unsigned n=current->ElementCount();
	
	// These are empirical limits on how many elements are needed before the solver becomes stable
	// using double precision. For higher-precision solvers more might be possible, but the value of doing
	// so is debatable, as it would require very high precision tables.
	static const unsigned stability[]={
		1,					// linear is always fine
		8,					// cubic is very stable, and goes way down, but the curve is hugely distorted
		64,				// quintic is very stable, and works down to 32, but the distribution looks very distorted
		(1<<16),		// heptic, stable, but not convincing
		(1<<17)		// nonic, seems more stable than heptic (actually works at 2^14), but is very worrying
	};
	
	if(correction=="auto"){
		// These are empirically derived, and by no means optimal
		if(n>=2*stability[4]){
			correction="poly9";
		}else if(n>=2*stability[3]){
			correction="poly7";
		}else if(n>=2*stability[2]){
			correction="poly5";
		}else if(n>=2*stability[1]){
			correction="poly3";
		}else{
			correction="poly1";
		}
		if(::flopoco::verbose>=INFO)
			std::cerr<<"CorrectTable : chose correction="<<correction<<" for table with "<<n<<" elements.\n";
	}
	
	if(correction=="none"){
		std::cerr<<"WARNING : CorrectTable is using correction=none, this should only be used if you intend to create broken GRNGs.\n";
		return current;
	}else if(correction.substr(0,4)=="poly"){
		std::string index=correction.substr(4,-1);
		int degree=atoi(index.c_str());
		
		if(::flopoco::verbose>=DETAILED)
			std::cerr<<"CorrectTable : got correction="<<correction<<" for table with "<<n<<" elements.\n";
		
		if((degree<1) || (degree>9) || ((degree%2)==0))
			throw std::string("CorrectTable - Degree (N) for polyN must be odd and between 1 and 9.");
		
		assert(current->IsSymmetric() && targetDist->IsSymmetric());
		
		// Now we know the degree, look for the correction
		std::vector<T> poly=FindPolynomialCorrection<T>(current, targetDist, degree);
		if(::flopoco::verbose>=INFO){
			std::cerr<<"CorrectTable : correction poly = "<<poly[0];
			for(unsigned i=1;i<poly.size();i++){
				std::cerr<<"+x^"<<i<<"*"<<poly[i];
			}
			std::cerr<<"\n";
		}
		
		if(::flopoco::verbose>=DETAILED)
			std::cerr<<"CorrectTable :applying correction.\n";
		return current->ApplyPolynomial(poly);
	}else{
		throw std::string("CorrectTable : Didn't understand correction method '"+correction+"'.");
	}
}

template<class T>
T square(const T &x)
{ return x*x; }

template<class T>
typename TableDistribution<T>::TypePtr QuantiseTable(
	std::string correction,
	typename ContinuousDistribution<T>::TypePtr targetDistrib,	//! Target distribution
	typename TableDistribution<T>::TypePtr current,	//! discrete distribution with continuous elements
	int wF	//! Number of fractional bits to quantise to
){
	int n=current->ElementCount();
	
	if(correction=="auto"){
		correction="stddev_greedy";
	}
	
	std::vector<T> contents(current->ElementCount());
	for(int i=0;i<n;i++){
		contents[i]=current->RangeFromIndex(i);
	}
	
	if(::flopoco::verbose>=INFO)
			std::cerr<<"QuantiseTable : applying quantisation '"+correction+"'.\n";
	if(correction=="round"){
		for(int i=0;i<n/2;i++){
			contents.at(i)=ldexp(round(ldexp(contents.at(i), wF)),-wF);
			contents.at(n-i-1)=-contents[i];
		}
	}else if(correction=="stddev_greedy"){
		assert((n%2)==0);
		T delta=pow(2.0, -wF);
		T target=targetDistrib->StandardMoment(2);
		T acc=0;
		for(int i=n/2;i<n;i++){
			contents[i]=ldexp(round(ldexp(contents.at(i), wF)),-wF);
			contents[n-i-1]=-contents[i];
			acc += contents[i] * contents[i];
		}
		T curr=sqrt(acc/(n/2));
		
		for(int i=n-1;i>=n/2;i--){
			//std::cerr<<"  curr="<<curr<<", target="<<target<<", err="<<(curr-target)/target<<"\n";
			while(contents[i]>0){
				//T accDown=acc-square(contents[i])+square(contents[i]-eps);
				//T accDown=acc-contents[i]*contents+(contents[i]-eps)*(contents[i]-eps);
				//T accDown=acc-contents[i]^2+contents[i]^2-2*eps*contents[i]+eps^2;
				T accDown=acc-2*delta*contents[i]+square(delta);
				T gotDown=sqrt(accDown/(n/2));
				if(abs(gotDown-target) < abs(curr-target)){
					contents[i] -= delta;
					acc=accDown;
					curr=gotDown;
				}else{
					break;
				}
			}
			while(true){
				T accUp=acc+2*delta*contents[i]-square(delta);
				T gotUp=sqrt(accUp/(n/2));
				if(abs(gotUp-target) < abs(curr-target)){
					contents[i] += delta;
					acc=accUp;
					curr=gotUp;
				}else{
					break;
				}
			}
		}
		for(int i=0;i<n/2;i++){
			contents[i]=-contents[n-i-1];
		}
		
		
	}else{
		throw std::string("QuantiseTable - Unknown quantisation method '"+correction+"'");
	}
	
	if(::flopoco::verbose>=DETAILED)
		std::cerr<<"CorrectTable : rebuilding distribution.\n";
	return boost::make_shared<TableDistribution<T> >(&contents[0], &contents[n]);
}

template<class T>
std::vector<mpz_class> BuildTable(
	int k, T stddev, int wF,
	std::string correction, std::string quantisation
){
	ContinuousDistribution<double>::TypePtr targetDist=boost::make_shared<GaussianDistribution<T> >(0, stddev);
	
	int n=1<<k;
	std::vector<T> x(n);
	T half=0.5;
	for(int i=0;i<n/2;i++){
		T u=(i+half)/n;
		x[i]=targetDist->InvCdf(u);
		x[n-i-1]=-x[i];
	}
	
	typename TableDistribution<T>::TypePtr original=boost::make_shared<TableDistribution<double> >(&x[0], &x[n]);
	typename TableDistribution<T>::TypePtr corrected=CorrectTable<T>(correction, targetDist, original);
	typename TableDistribution<T>::TypePtr quantised=QuantiseTable<T>(quantisation, targetDist, corrected, wF);
	
	if(::flopoco::verbose>=DETAILED){
		std::cerr<<"  metric, raw, corrected, quantised, target\n";
		for(int i=2;i<=12;i++){
			T want=targetDist->StandardMoment(i), got=quantised->StandardMoment(i);
			std::cerr<<" mu_"<<i<<", "<<original->StandardMoment(i)<<", "<<corrected->StandardMoment(i)<<", "<<got<<", "<<want<<", "<<(got-want)/want<<"\n";
		}
	}
	
	if(::flopoco::verbose>=DETAILED)
		std::cerr<<"BuildTable : converting table to fixed-point.\n";
	std::vector<mpz_class> res(n/2);
	T scale=pow(2.0, wF);
	
	std::vector<std::pair<T,T> > elts=quantised->GetElements();
	for(unsigned i=0;i<res.size();i++){
		res[i]=(mpz_class)round(scale*elts[i+n/2].first);
	}
	
	return res;
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
	dst << "            none - Don't correct at all (not advised; at least do poly1 to fix the stddev).\n";
	dst << "            poly1,poly3,...,poly9 - Apply odd polynomial correction of given degree.\n";
	dst << "          fixMethod - How to convert the continous table to fixed.\n";
	dst << "            auto - Select some sensible fixation method.\n";
	dst << "            round - Simply round the continuous entries.\n";
	dst << "      The transform will take a (k+1) bit uniform input, and produce a (w+1) bit signed output.\n";
}

TableTransform *MakeGRNGTable(Target *target, int k, int wF, mpfr::mpreal stddev, std::string correction, std::string quantisation)
{
	if(k<1)
		throw std::string("TableFactory - k must be a positive integer.");
	if(wF<1)
		throw std::string("TableFactory - w must be a positive integer.");
	
	if(wF>32)
		throw std::string("TableFactory - w must be less than 32 (currently table is built in double-precision).");
	
	std::vector<mpz_class> contents=BuildTable<double>(k, stddev.toDouble(), wF, correction, quantisation);
	unsigned wO=0;
	for(unsigned i=0;i<contents.size();i++){
		unsigned ww=mpz_sizeinbase(contents[i].get_mpz_t(), 2);
		if(::flopoco::verbose>=DEBUG){
			std::cerr<<"table["<<i<<"]="<<contents[i]<<", w="<<ww<<"\n";
		}
		wO=std::max(wO, ww);
	}
	
	TableTransform *result=new TableTransform(target, wO, contents, true, wF);
	
	return result;
}

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


