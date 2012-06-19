#ifndef flopoco_random_fixed_point_exp_find_close_value_hpp
#define flopoco_random_fixed_point_exp_find_close_value_hpp

#include "residual_type.hpp"

#include <map>

namespace flopoco{
namespace random{

template<class T>	
struct close_value_t
{
	T input;
	T offset;
	T exact;	// == exp(input-offset)
	T approx;	// == round(exact,fracBits)
	unsigned fracBits;
};

template<class T>
close_value_t<T> FindCloseValue(
	T mu, T sigma,
	T input,
	const residual_type<T> &offsetSpace,
	T errorMin,
	T errorMax,
	unsigned fracBits // Current number of fracBits
){
	close_value_t<T> res;
	res.input=input;
	
	while(true){
		typename residual_type<T>::iterator curr=offsetSpace.begin();
		typename residual_type<T>::iterator end=offsetSpace.end();
	
		if(fracBits+1 >= (unsigned)std::numeric_limits<T>::digits){
			std::stringstream acc;
			acc<<"FindCloseValue - fracBits has risen to level of type fraction width, fracBits="<<fracBits<<", digits="<<std::numeric_limits<T>::digits<<"\n";
			throw std::invalid_argument(acc.str());
		}
		
		while(curr!=end){
			
			T offset=*curr;
			
			T ref=exp(mu+sigma*(input-offset));
			int e;
			T refF=frexp(ref,&e);
			T gotF=ldexp(ceil(ldexp(refF, fracBits+1)), -(int)(fracBits+1));
			T got=ldexp(gotF, e);
			T err=(got-ref)/ref;
			
			if((errorMin<=err) && (err<=errorMax)){
				res.fracBits=fracBits;
				res.offset=offset;
				res.exact=ref;
				res.approx=got;
				return res;
			}
			
			gotF=ldexp(floor(ldexp(refF, fracBits+1)), -(int)(fracBits+1));
			got=ldexp(gotF, e);
			err=(got-ref)/ref;
			
			if((errorMin<=err) && (err<=errorMax)){
				res.fracBits=fracBits;
				res.offset=offset;
				res.exact=ref;
				res.approx=got;
				return res;
			}
		
			++curr;
		}
		++fracBits;
	}
}

template<class T>
struct close_value_table_t
{
	T mu, sigma;
	residual_type<T> inputSpace;
	residual_type<T> targetOffsetSpace;
	T targetErrorMin, targetErrorMax;
	residual_type<T> offsetSpace;
	result_type<T> resultSpace;
	std::map<T,close_value_t<T> > table;
};

/* Run through the entire input space, and find the smallest number of
	fracBits that will support the given relative error, and a mapping
	from values in the inputSpace to an offset*/
template<class T>
close_value_table_t<T> FindCloseValues(
	T mu, T sigma,
	const residual_type<T> &inputSpace,
	const residual_type<T> &offsetSpace,
	T errorMin,
	T errorMax
){
	close_value_table_t<T> res={
		mu, sigma,
		inputSpace, offsetSpace,
		errorMin, errorMax,
		offsetSpace,
		result_type<T>(0),
		std::map<T,close_value_t<T> >()
	};
	
	unsigned fracBits=0;
	
	typename residual_type<T>::iterator curr=inputSpace.begin();
	typename residual_type<T>::iterator end=inputSpace.end();
	
	//std::cerr<<"errMin="<<errorMin<<", errMax="<<errorMax<<"\n";
	
	while(curr!=end){
		close_value_t<T> v=FindCloseValue(mu, sigma, *curr, offsetSpace, errorMin, errorMax, fracBits);
		fracBits=std::max(fracBits, v.fracBits);
		//std::cerr<<"  in="<<*curr<<", fracBits="<<v.fracBits<<", relErr="<<(v.approx-v.exact)/v.exact<<"\n";
		res.table[*curr]=v;
		res.offsetSpace.Add(v.offset);
		++curr;
	}
	
	res.resultSpace=result_type<T>(fracBits);
	
	typename std::map<T,close_value_t<T> >::iterator it=res.table.begin();
	while(it!=res.table.end()){
		res.resultSpace.Add(it->second.approx, it->second.exact);
		++it;
	}
	
	return res;
}

}; // random
}; // flopoco

#endif
