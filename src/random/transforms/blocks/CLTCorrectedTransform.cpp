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
#include "random/distributions/histogram_distribution.hpp"
#include "random/distributions/sum_distribution.hpp"

#include "random/utils/fft/convolve_mpreal.hpp"

using namespace std;

namespace flopoco
{

extern vector<Operator *> oplist;	
	
namespace random
{

class CLTCorrectedTransform
	: public RngTransformOperator
	, public  IRngTransformDistributions
{
	mutable Distribution<mpfr::mpreal>::TypePtr m_distribution;
	
public:
	CLTCorrectedTransform(Target *target)
		: RngTransformOperator(target)
		
	{};
	
	virtual void outputVHDL(std::ostream& o, std::string name)
	{ throw std::string("CLTCorrectedTransform::outputVHDL - Not implemented (distribution only at the moment)."); }

	virtual unsigned uniformInputBits() const
	{ return 128; }
	
	virtual std::string uniformInputName() const
	{ return "urng"; }

	virtual unsigned nonUniformOutputCount() const
	{ return 1; }

	virtual bool nonUniformOutputsAreHomogenous() const
	{ return true; }

	virtual unsigned nonUniformOutputWidth(int i) const
	{ return 16; }

	//! The name of the output port associated with this signal
	virtual std::string nonUniformOutputName(int i) const
	{ return "X"; }
	
	struct half_op
	{
		int fb;
		
		half_op(int _fb)
			: fb(_fb)
		{}
			
		mpfr::mpreal operator()(mpfr::mpreal x) const
		{
			return ldexp(floor(ldexp(x,fb)),-fb);
		}
	};
	
	/*struct corrector_op
	{
		int fb;
		std::vector<std::pair<double,double> > table;
		
		corrector_op(int _fb)
			: fb(_fb)
		{
			
		}
		
		mpfr::mpreal operator()(mpfr::mpreal x) const
		{
			
		}
	};*/
	/*
	struct pseudo_minimax
	{
		std::vector<double> x;	// Values in the input 
		std::vector<double> y;	// Target values in the output
		
		// We have k segments, giving k+1 endpoints, defined in terms of the indices in the domain, and the endpoints in the range
		std::vector<int> indices;	// k+1 indices
		std::vector<double> endpoints; // k indices
		
		std::vector<double> errors;
		
		std::pair<double,double> GetTransform(int i)
		{
			int start_index=indices[i], finish_index=indices[i+1]-1;
			
			// We want:
			//	x[start_index]*scale+offset=endpoints[i]
			//	x[finish_index]*scale+offset=endpoints[i+1];
			//So:
			double scale=(endpoints[i]-endpoints[i+1]) / (x[start_index]-x[finish_index]);
			double offset=(endpoints[i+1]*x[start_index] - x[finish_index]*endpoints[i]) / / (x[start_index]-x[finish_index]);
			return std::make_pair(scale,offset);
		}
		
		double SegmentError(int i, std::pair<double,double> transform) const
		{
			int start_index=indices[i], finish_index=indices[i+1]-1;
			
			double worst=0;
			for(int i=start_index;i<end_index;i++){
				double err=std::fabs(y[i] - x[i]*transform.first+transform.second);
				worst=std::max(worst, err);
			}
			return worst;
		}
		
		void PerturbEndpoint(int boundary)
		{
			
		}
	};
	*/
	
	virtual Distribution<mpfr::mpreal>::TypePtr nonUniformOutputDistribution(int i, unsigned prec) const
	{
		if(!m_distribution){
			// First, build up the CLT distribution. The stated technique is to start
			// from 16 bit samples, then to add, but every time they add one LSB is
			// dropped. Worryingly, this seems to imply quite a lot of bias, as the
			// thing won't be zero mean to start with, and is going to get dragged
			// down with every lost LSB.
			
			mpfr::mpreal one(1.0, prec);
			
			int fb=12; // 16 for full version
			
			std::vector<mpfr::mpreal> ones(1<<fb, ldexp(one, -fb));
			
			// First distribution in Q(16,15)
			HistogramDistribution<mpfr::mpreal>::TypePtr curr=boost::make_shared<HistogramDistribution<mpfr::mpreal> >(-one, ldexp(one, -(fb-1)), ones);
			
			assert(curr->Pmf(0).get_prec()>=(int)prec);
			
			// First adder and truncation to Q(16,14) for CLT-2
			REPORT(INFO,"Building CLT-2");
			curr=SelfAddHistogramDistributions<mpfr::mpreal>(curr, 2);
			curr=HistogramDistribution<mpfr::mpreal>::BuildFromTransform(curr, half_op(fb-2), mpfr::mpreal(pow(2.0,-(fb-2)),prec));
			
			// Second adder and truncation to Q(16,13) for CLT-4
			REPORT(INFO,"Building CLT-4");
			curr=SelfAddHistogramDistributions<mpfr::mpreal>(curr, 2);
			curr=HistogramDistribution<mpfr::mpreal>::BuildFromTransform(curr, half_op(fb-3), mpfr::mpreal(pow(2.0,-(fb-3)),prec));
			
			// Final adder and truncation to Q(16,12) for CLT-8
			REPORT(INFO,"Building CLT-8");
			curr=SelfAddHistogramDistributions<mpfr::mpreal>(curr, 2);
			curr=HistogramDistribution<mpfr::mpreal>::BuildFromTransform(curr, half_op(fb-4), mpfr::mpreal(pow(2.0,-(fb-4)),prec));
			
			for(int i=0;i<=10;i++){
				std::cerr<<"Mom["<<i<<"] = "<<curr->StandardMoment(i)<<"\n";
			}
			
			/* From the paper we can't work out exactly what the non-uniform segmentation
				is, as they only list the coefficients and not the boundaries. The whole curvature
				method for segmentation is not really explained, so I'm going to find a psuedo-minimax,
				as that seems to be their goal.
			
				Their error metric seems to be to look at the PMF for a particular point, then to
				invert it back to the place where it should be. I have no idea what that is supposed
				to mean, as it is mixing PMF and PDF, so I'm going to look at CDF.
			*/
			std::pair<int64_t,int64_t> support=curr->IndexSupport();
			support.second=curr->IndexFromRange(0);
			int64_t points=support.second-support.first+1;
			std::vector<mpfr::mpreal> tmp(points);
			
			curr->RangeByIndex(support.first, support.second, &tmp[0]);
			std::vector<double> x_clt(points);	// The points of the untransformed CLT
			for(int i=0;i<points;i++){
				x_clt[i]=tmp[i].toDouble();
			}
			
			curr->CdfByIndex(support.first, support.second, &tmp[0]);
			std::vector<double> p_clt(points);	// The CDF of the untransformed CLT
			for(int i=0;i<points;i++){
				p_clt[i]=tmp[i].toDouble();
			}
			
			boost::math::normal_distribution<double> n(0.0,1.0);
			
			std::vector<double> x_gauss(points);
			for(int i=0;i<points;i++){
				std::cerr<<x_clt[i]<<" -> ";
				if(p_clt[i]==0)
					x_gauss[i]=0;
				else
					x_gauss[i]=boost::math::quantile(n, p_clt[i]);
				std::cerr<<p_clt[i]<<" -> "<<x_gauss[i]<<"\n";
			}
			
			// Ok, now we have x_clt, which is where we are, and x_gauss, which is where we want to be
			
			m_distribution=curr;
		}
		return m_distribution;
	}
};	

static void CLTCorrectedFactoryUsage(std::ostream &dst)
{
	OperatorFactory::classic_OP(dst, "clt_corrected_transform", "", false);
	dst<<"    Builds transform from 'Generating High Tail Accuracy\n";
	dst<<"    Gaussian Random Numbers in Hardware Using Central Limit Theorem'.\n";
	dst<<"    Currently no actual VHDL operator is built, just distribution.\n";
}

static Operator *CLTCorrectedFactoryParser(Target *target ,const std::vector<std::string> &args,int &consumed)
{
	consumed =0;
	
	return new CLTCorrectedTransform(target);
}

void CLTCorrectedTransform_registerFactory()
{
	DefaultOperatorFactory::Register(
		"clt_corrected_transform",
		"operator;rng_transform",
		CLTCorrectedFactoryUsage,
		CLTCorrectedFactoryParser
	);
}

};
};


