#ifndef flopoco_random_distributions_sum_distribution_hpp
#define flopoco_random_distributions_sum_distribution_hpp

#include "random/distributions/table_distribution.hpp"
#include "random/distributions/histogram_distribution.hpp"

#include <map>

namespace flopoco
{
namespace random
{
	
	template<class T>
	std::vector<T> self_convolve(const std::vector<T> &values, int n);
	
	template<class T>
	std::vector<T> convolve(const std::vector<T> &a, const std::vector<T> &b);
	
	template<class T,class TF>
	typename EnumerableDistribution<T>::TypePtr ConvolveDistributions(
		typename EnumerableDistribution<T>::TypePtr a,
		typename EnumerableDistribution<T>::TypePtr b,
		const TF &f
	){
		if(a->ElementCount() * b->ElementCount() > pow(2.0,40))
			throw std::runtime_error("Attempt to perform additive convolution with more than 2^40 combinations.");
		
		std::vector<std::pair<T,T> > va, vb;
		a->GetElements(0, a->ElementCount(), &va[0]);
		b->GetElements(0, b->ElementCount(), &vb[0]);
				
		std::map<T,T> acc;	// should really be a hashtable
		
		for(size_t ia=0;ia<va.size();ia++){
			std::pair<T,T> ca=va[ia];
			if(ca.second==0)
				continue;
			for(size_t ib=0;ib<vb.size();ib++){
				std::pair<T,T> cb;
				acc[ f(ca.first,cb.first) ] += ca.second * cb.second;
			}
		}
		
		return boost::make_shared<TableDistribution<T> >(
			acc
		);
	}
	
	template<class T>
	typename HistogramDistribution<T>::TypePtr AddHistogramDistributions(
		typename HistogramDistribution<T>::TypePtr a,
		typename HistogramDistribution<T>::TypePtr b
	){
		if(a->RangeStep()!=b->RangeStep())
			throw std::runtime_error("AddDistributions - histograms with unequal ranges are not supported yet.");
		
		std::vector<T> va(a->ElementCount()), vb(b->ElementCount());
		a->GetProbabilities(0, a->ElementCount(), &va[0]);
		b->GetProbabilities(0, b->ElementCount(), &vb[0]);
		
		std::vector<T> vr=convolve(va,vb);
		
		return boost::make_shared<HistogramDistribution<T> >(
			a->RangeBase()+b->RangeBase(),
			a->RangeStep(),
			vr
		);
	}
	
	template<class T>
	typename HistogramDistribution<T>::TypePtr SelfAddHistogramDistributions(
		typename HistogramDistribution<T>::TypePtr a,
		int k
	){
		std::vector<T> va(a->ElementCount());
		a->GetProbabilities(0, a->ElementCount(), &va[0]);
		
		std::vector<T> vr=self_convolve(va, k);
		
		return boost::make_shared<HistogramDistribution<T> >(
			k*a->RangeBase(),
			a->RangeStep(),
			vr
		);
	}
	
	template<class T>
	typename EnumerableDistribution<T>::TypePtr AddDistributions(
		typename EnumerableDistribution<T>::TypePtr a,
		typename EnumerableDistribution<T>::TypePtr b
	){
		typename HistogramDistribution<T>::TypePtr ha=boost::dynamic_pointer_cast<HistogramDistribution<T> >(a);
		typename HistogramDistribution<T>::TypePtr hb=boost::dynamic_pointer_cast<HistogramDistribution<T> >(b);
		
		if(ha&&hb){
			if(ha->RangeStep()==hb->RangeStep())
				return AddHistogramDistributions<T>(ha,hb);
		}
		
		if(a->ElementCount() < b->ElementCount())
			return AddDistributions<T>(b,a);
		return ConvolveDistributions<T>(a,b,std::plus<T>());
	}
	
	template<class T>
	typename EnumerableDistribution<T>::TypePtr SelfAddDistributions(
		typename EnumerableDistribution<T>::TypePtr a,
		int k
	){
		if(k<=0)
			throw std::invalid_argument("SelfAddDistributions - Can't convolve less than once.");
		if(k==1)
			return a;
		
		typename HistogramDistribution<T>::TypePtr ha=boost::dynamic_pointer_cast<HistogramDistribution<T> >(a);
		
		if(ha){
			return SelfAddHistogramDistributions<T>(ha,k);
		}
		
		typename EnumerableDistribution<T>::TypePtr res;
		typename EnumerableDistribution<T>::TypePtr ss=a;
		while(k){
			if(k&1){
				if(!res)
					res=ss;
				else
					res=AddDistributions<T>(res, ss);
			}
			k=k/2;
			if(k==0)
				break;
			ss=AddDistributions<T>(ss, ss);
		}
		
		return res;
	}
	
	template<class T>
	typename EnumerableDistribution<T>::TypePtr SubtractDistributions(
		typename EnumerableDistribution<T>::TypePtr a,
		typename EnumerableDistribution<T>::TypePtr b
	){
		if(a->ElementCount() < b->ElementCount)
			return AddDistributions(b,a);
		return ConvolveDistributions(a,b,std::minus<T>());
	}
	
};
};

#endif
