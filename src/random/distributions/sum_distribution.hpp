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
	std::vector<T> self_convolve(const std::vector<T> &values, unsigned n);
	
	template<class T>
	std::vector<T> convolve(const std::vector<T> &a, const std::vector<T> &b);
	
	template<class T,class TF>
	typename DiscreteDistribution<T>::TypePtr ConvolveDistributions(
		typename DiscreteDistribution<T>::TypePtr a,
		typename DiscreteDistribution<T>::TypePtr b,
		const TF &f
	){
		if(a->ElementCount() * b->ElementCount() > pow(2.0,40))
			throw std::runtime_error("Attempt to perform additive convolution with more than 2^40 combinations.");
		
		std::pair<int64_t,int64_t> aIndices=a->IndexSupport();
		std::pair<int64_t,int64_t> bIndices=b->IndexSupport();
		
		std::vector<std::pair<T,T> > va(aIndices.second-aIndices.first), vb(bIndices.second-bIndices.first);
		a->ElementsByIndex(aIndices.first, aIndices.second, &va[0]);
		b->ElementsByIndex(bIndices.first, bIndices.second, &vb[0]);
				
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
		a->PmfByIndex(0, a->ElementCount(), &va[0]);
		b->PmfByIndex(0, b->ElementCount(), &vb[0]);
		
		std::vector<T> vr=convolve(va,vb);
		
		if(a->IsSymmetric() && b->IsSymmetric()){
			for(int i=0;i<(int)vr.size()/2;i++){
				vr[vr.size()-i-1]=vr[i];
			}
		}
		
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
		a->PmfByIndex(0, a->ElementCount(), &va[0]);
		
		if(a->IsSymmetric()){
			for(int i=0;i<(int)va.size()/2;i++){
				assert(va[i]==va[va.size()-i-1]);
			}
		}
		
		std::vector<T> vr=self_convolve(va, k);
		
		if(a->IsSymmetric()){
			T origSum=sum(vr.begin(), vr.end());
			
			for(int i=0;i<(int)vr.size()/2;i++){
				std::cerr<<i<<" / "<<(vr.size()-i-1)<<" = "<<vr[i]<<" / "<<vr[vr.size()-i-1]<<"\n";
				T avg=(vr[vr.size()-i-1]+vr[i])>>1;
				vr[vr.size()-i-1]=avg;
				vr[i]=avg;
			}
			assert(vr.front()==vr.back());
			
			T currSum=sum(vr.begin(), vr.end());
			
			std::cerr<<"OrigSum="<<origSum<<", currSum="<<currSum<<"\n";
		}
		
		for(int i=0;i<(int)vr.size();i++){
			if(vr[i]<0)
				vr[i]=0;
		}
		
		return boost::make_shared<HistogramDistribution<T> >(
			k*a->RangeBase(),
			a->RangeStep(),
			vr
		);
	}
	
	template<class T>
	typename DiscreteDistribution<T>::TypePtr AddDistributions(
		typename DiscreteDistribution<T>::TypePtr a,
		typename DiscreteDistribution<T>::TypePtr b
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
	typename DiscreteDistribution<T>::TypePtr SelfAddDistributions(
		typename DiscreteDistribution<T>::TypePtr a,
		int k
	){
		if(k<=0)
			throw std::invalid_argument("SelfAddDistributions - Can't convolve less than once.");
		if(k==1)
			return a;
		
		typename HistogramDistribution<T>::TypePtr ha=boost::dynamic_pointer_cast<HistogramDistribution<T> >(a);
		typename TableDistribution<T>::TypePtr ta=boost::dynamic_pointer_cast<TableDistribution<T> >(a);
		
		if( !!ta && (pow((double)a->ElementCount(), k) > 100000)){
			int fb=ta->FixedPointResolution();
			if(fb != INT_MAX){
				std::vector<std::pair<T,T> > elts=ta->GetElements();
				
				int iFirst=boost::math::tools::real_cast<long>(ldexp(elts.front().first,fb)), iLast=boost::math::tools::real_cast<long>(ldexp(elts.back().first,fb));
				std::vector<T> pdf(iLast-iFirst+1);
				
				for(int i=0;i<(int)elts.size();i++){
					pdf[boost::math::tools::real_cast<long>(ldexp(elts[i].first,fb)-iFirst)] += elts[i].second;
				}
				
				ha=boost::make_shared<HistogramDistribution<T> >(elts.front().first, T(pow(2.0,-fb)), pdf);
			}
		}
		
		if(ha){
			return SelfAddHistogramDistributions<T>(ha,k);
		}
		
		typename DiscreteDistribution<T>::TypePtr res;
		typename DiscreteDistribution<T>::TypePtr ss=a;
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
	
	/*
	template<class T>
	typename DiscreteDistribution<T>::TypePtr SubtractDistributions(
		typename DiscreteDistribution<T>::TypePtr a,
		typename DiscreteDistribution<T>::TypePtr b
	){
		return ConvolveDistributions(a,b,std::minus<T>());
	}
	*/
	
};
};

#endif
