#include <iostream>
#include <sstream>
#include <math.h>

#include "RngTransformOperator.hpp"

#include "random/distributions/discrete_gaussian_distribution.hpp"
#include "random/distributions/quantised_gaussian_distribution.hpp"

#include "random/utils/operator_factory.hpp"

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/make_shared.hpp>

#include "UtilSollya.hh"

using namespace std;

namespace flopoco
{

extern vector<Operator *> oplist;	
	
namespace random
{
	
template<class T>
void DumpCdfStats(std::ostream &dest_file, typename DiscreteDistribution<T>::TypePtr got, typename DiscreteDistribution<T>::TypePtr target, T startx)
{
	if(!(target->IsSymmetric()&&got->IsSymmetric()))
		throw std::string("DumpStats - Currently only strictly symmetric distributions are supported.");
	
	if(!(target->StandardMoment(1)==0 &&got->StandardMoment(1)==0))
		throw std::string("DumpStats - Currently only zero mean distributions are supported.");
	
	int64_t begin_index_target=target->ClosestIndexFromRange(startx), end_index_target=target->ClosestIndexFromRange(0)+1;
	int64_t n=end_index_target-begin_index_target;
	
	std::vector<T> cdf_target(n), range_target(n);
	
	target->RangeByIndex(begin_index_target, end_index_target, &range_target[0]);
	target->CdfByIndex(begin_index_target, end_index_target, &cdf_target[0]);
	
	T worstAbsCdfErr=0, worstRelCdfErr=0;
	T sumAbsCdfErr=0, sumRelCdfErr=0;
	
	for(int i=0;i<n;i++){
		T p=got->Cdf(range_target[i]);
		
		T absErr=abs(p-cdf_target[i]);
		T relErr=absErr/cdf_target[i];
		
		//if(absErr>worstAbsCdfErr || relErr>worstRelCdfErr)
		//std::cerr<<i<<", "<<range_target[i]<<", "<<cdf_target[i]<<", "<<p<<", "<<absErr<<", "<<relErr<<"\n";
		
		worstAbsCdfErr=std::max(worstAbsCdfErr, absErr);
		worstRelCdfErr=std::max(worstRelCdfErr, relErr);
	}
	
	dest_file<<"MaxCdfErr, "<<startx<<"..0, rel="<<worstRelCdfErr<<", abs="<<worstAbsCdfErr<<"\n";
}

template<class T>
void DumpMomentStats(std::ostream &dest_file, typename DiscreteDistribution<T>::TypePtr got, typename DiscreteDistribution<T>::TypePtr target)
{
	for(int i=2;i<=8;i+=2){
		T e=target->StandardMoment(i), o=got->StandardMoment(i);
		std::cerr<<"m["<<i<<"], target="<<e<<", got="<<o<<", rel="<<(o-e)/e<<"\n";
	}
}

static void TransformStatsUsage(std::ostream &dst)
{
	OperatorFactory::classic_OP(dst, "TransformStats", "distribution", false);
	dst << "       Calculates the accuracy of the previous statistical transform operator.\n";
	dst <<"          distribution - Specifies the distribution to check against.\n";
	dst<<"    distributions are:\n";
	dst<<"      QuantisedGaussian[stddev,fracbits]\n";
	dst<<"      DiscreteGaussian[stddev,fracbits]\n";
}

DiscreteDistribution<mpfr::mpreal>::TypePtr ParseDistribution(std::string spec, int prec)
{
	std::vector<std::string> parts;
	boost::split(parts, spec, boost::is_any_of("[,]"));
	
	if(parts.size()==0)
		throw std::string("ParseDistribution('")+spec+"') - Doesn't contain any parts.";
	
	for(int i=0;i<parts.size();i++){
		std::cerr<<i<<" : '"<<parts[i]<<"', len="<<parts[i].size()<<"\n";
	}
	
	if(parts[0]=="QuantisedGaussian" || parts[0]=="DiscreteGaussian"){
		if(parts.size()<3)
			throw std::string("ParseDistribution('")+spec+"') - Need exactly three parts for this distribution.";
		mpfr::mpreal stddev(0, prec);
		parseSollyaConstant(get_mpfr_ptr(stddev), parts[1], MPFR_RNDN);
		int fracbits=boost::lexical_cast<int>(parts[2]);
		
		if(parts[0]=="QuantisedGaussian"){
			return boost::make_shared<QuantisedGaussianDistribution<mpfr::mpreal> >(stddev, fracbits);
		}else{
			return boost::make_shared<DiscreteGaussianDistribution<mpfr::mpreal> >(stddev, fracbits);
		}
	}else{
		throw std::string("ParseDistribution('")+spec+"') - Didn't understand first part '"+parts[0]+"'.";
	}
}

static Operator *TransformStatsParser(Target *target ,const std::vector<std::string> &args,int &consumed)
{
	if(args.size()<1)
		throw std::string("TransformStatsParser - No distribution specification.");
	consumed = 1;
	
	int prec=128;
	
	DiscreteDistribution<mpfr::mpreal>::TypePtr targetDist=ParseDistribution(args[0], prec);
	
	vector<Operator*> *ops=target->getGlobalOpListRef();
	if(ops->size()==0)
		throw std::string("TransformStatsParser - No previous operator to work with.");
	
	Operator *op=ops->back();
	
	RngTransformOperator *transform=dynamic_cast<RngTransformOperator*>(op);
	if(transform==NULL)
		throw std::string("TransformStatsParser - Previous operator is not an RngTransformOperator.");
	
	IRngTransformDistributions *opDists=dynamic_cast<IRngTransformDistributions*>(op);
	if(opDists==NULL)
		throw std::string("TransformStatsParser - Previous operator does not support IRngTransformDistributions.");
	
	if(!transform->nonUniformOutputsAreHomogenous() && transform->nonUniformOutputCount() > 1)
		throw std::string("TransformStatsParser - Can't deal with multiple non-homogeneous output distributions.");
	
	Distribution<mpfr::mpreal>::TypePtr got=opDists->nonUniformOutputDistribution(0, prec);
	DiscreteDistribution<mpfr::mpreal>::TypePtr gotDist=boost::dynamic_pointer_cast<DiscreteDistribution<mpfr::mpreal> >(got);
	
	DumpCdfStats(std::cout, gotDist, targetDist, -8 * targetDist->StandardMoment(2));
	DumpMomentStats<mpfr::mpreal>(std::cout, gotDist, targetDist);
	
	return NULL;
}

void TransformStats_registerFactory()
{
	DefaultOperatorFactory::Register(
		"transform_stats",
		"utility",
		TransformStatsUsage,
		TransformStatsParser
	);
}

};
};
