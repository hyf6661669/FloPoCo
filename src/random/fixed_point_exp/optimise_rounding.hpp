#ifndef flopoco_random_optimise_rounding_hpp
#define flopoco_random_optimise_rounding_hpp

#include "random/fixed_point_exp/result_type.hpp"
#include "random/fixed_point_exp/fixed_point_exp_stage.hpp"

namespace flopoco{
namespace random{

/*! Describes how to multiply two floating-point fractions together */
template<class T>
struct multiplier_rounding_spec
{
	multiplier_rounding_spec()
		: resultType(0)
	{}
	
	bool roundHalfUp;			// If true then round half up, otherwise round half down
	bool postNormUp;		// True if the result is ever normalised up
	bool postNormDown;	// True if the result is ever normalised down
	result_type<T> resultType;
};

namespace detail{
	template<class T>
	int getExpnt(T x)
	{
		int e;
		frexp(x, &e);
		return e;
	}
};

template<class TItA, class TItB, class T>
multiplier_rounding_spec<T> MultiplyWithResultBits(
	const TItA &a,
	const TItB &b,
	int outputFracBits,
	bool roundHalfUp,
	T relErrorMin,
	T relErrorMax
){
	using detail::getExpnt;
	
	if(a.size() < b.size())
		return MultiplyWithResultBits(b,a,outputFracBits, roundHalfUp, relErrorMin, relErrorMax);
	
	std::string srcFileName="<none>";
	REPORT(DEBUG,"  MultiplyWithResultBits : outputFracBits = "<<outputFracBits);
	
	multiplier_rounding_spec<T> res;
	res.roundHalfUp=roundHalfUp;
	res.postNormUp=false;
	res.postNormDown=false;
	res.resultType=result_type<T>(outputFracBits);
	
	typename TItA::const_iterator curr_a=a.begin(), end_a=a.end();
	while(curr_a!=end_a){
		T exact_a=curr_a->first, approx_a=curr_a->second;
		
		typename TItB::const_iterator curr_b=b.begin(), end_b=b.end();
		while(curr_b!=end_b){
			T exact_b=curr_b->first, approx_b=curr_b->second;
			
			T exact=exact_a*exact_b;
			T approx=approx_a*approx_b;
			
			T hd=roundHalfUp ? res.resultType.RoundHalfUp(approx) : res.resultType.RoundHalfDown(approx);
			res.resultType.Add(hd, exact);
			
			REPORT(FULL, "expnt(a)="<<getExpnt(approx_a)<<", expnt(b)="<<getExpnt(approx_b)<<", expnt(hd)="<<getExpnt(hd));
			
			if(getExpnt(hd)==getExpnt(approx_a)+getExpnt(approx_b)){
				res.postNormUp=true;
			}else{
				res.postNormDown=true;
				assert(getExpnt(hd)+1==getExpnt(approx_a)+getExpnt(approx_b));
			}
			
			curr_b++;
		}
		
		if( (res.resultType.relErrorMin < relErrorMin) || (relErrorMax < res.resultType.relErrorMax)){
			throw unsatisfiable_constraints("MultiplyWithGuardBits - constraint not met.");
		}
		
		curr_a++;
	}
	
	REPORT(DEBUG, "  MultiplyWithResultBits : relErrorMin = "<<relErrorMin<<" <= ["<<res.resultType.relErrorMin<<","<<res.resultType.relErrorMax<<"] <= "<<relErrorMax<<" = relErrorMax" );
	
	return res;
}

template<class TItA, class TItB,class T>
multiplier_rounding_spec<T> FindResultBitsForMultiply(
	int wA,
	const TItA &a,
	int wB,
	const TItB &b,
	T relErrorMin,
	T relErrorMax
){
	std::string srcFileName="<none>";
	REPORT(DEBUG,"FindResultBitsForMultiply\n");
	
	for(int i=std::max(wA,wB);i<=wA+wB+1;i++){
		try{
			multiplier_rounding_spec<T> res=MultiplyWithResultBits<TItA,TItB,T>(a, b, i, true, relErrorMin, relErrorMax);
			return res;
		}catch(unsatisfiable_constraints){
			// try the next one
		}
		
		try{
			multiplier_rounding_spec<T> res=MultiplyWithResultBits<TItA,TItB,T>(a, b, i, false, relErrorMin, relErrorMax);
			return res;
		}catch(unsatisfiable_constraints){
			// try the next one
		}
	}	
	
	throw unsatisfiable_constraints("FindResultBitsForMultiply - constraints cannot be met.");
}

}; // random
}; // flopoco

#endif
