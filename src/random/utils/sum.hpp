#ifndef flopoco_random_utils_sum_hpp
#define flopoco_random_utils_sum_hpp

#include "random/utils/ladder_accumulator.hpp"

#include <iterator>

#include "mpreal.h"

namespace flopoco
{
namespace random
{
	template<class T>
	struct SelectAccumulator;
	
	template<>
	struct SelectAccumulator<float>
	{ typedef double type; };
	
	template<>
	struct SelectAccumulator<double>
	{ typedef LadderAccumulator type; };
	
	template<>
	struct SelectAccumulator<mpfr::mpreal>
	{ typedef mpfr::mpreal type; };
	
	template<class TIt>
	typename std::iterator_traits<TIt>::value_type sum(TIt begin, TIt end)
	{
		if(begin==end)
			return 0;
		if(begin+1==end)
			return *begin;
		typename SelectAccumulator<typename std::iterator_traits<TIt>::value_type>::type acc;
		acc=0;
		while(begin!=end){
			acc += *begin;
			++begin;
		}
		return acc;
	}
	
	template<class TIt, class TF>
	typename TF::result_type sum(TIt begin, TIt end, const TF &f)
	{
		if(begin==end)
			return 0;
		if(begin+1==end)
			return f(*begin);
		typename SelectAccumulator<typename TF::return_type>::type acc;
		acc=0;
		while(begin!=end){
			acc += f(*begin);
			++begin;
		}
		return acc;
	}
	
	template<class TIt>
	typename std::iterator_traits<TIt>::value_type sum_powers(TIt begin, TIt end, int k)
	{
		if(begin==end)
			return 0;
		if(begin+1==end)
			return pow(*begin, k);
		typename SelectAccumulator<typename std::iterator_traits<TIt>::value_type>::type acc;
		acc=0;
		while(begin!=end){
			acc += pow(*begin, k);
			++begin;
		}
		return acc;
	}
	
	template<class TIt>
	typename std::iterator_traits<TIt>::value_type::second_type
		sum_weighted_powers(TIt begin, TIt end, int k)
	{
		if(begin==end)
			return 0.0;
		if(begin+1==end)
			return pow(begin->first, k) * begin->second;
		typename SelectAccumulator<typename std::iterator_traits<TIt>::value_type::second_type>::type acc;
		acc=0;
		if(k==0){
			while(begin!=end){
				acc += begin->second;
				++begin;
			}
		}else if(k==1){
			while(begin!=end){
				acc += begin->first * begin->second;
				++begin;
			}
		}else{
			while(begin!=end){
				acc += pow(begin->first, k) * begin->second;
				++begin;
			}
		}
		return acc;
	}

}; // random
}; // flopoco

#endif
