#ifndef random_fixed_point_exp_result_type_hpp
#define random_fixed_point_exp_result_type_hpp

#include <limits>
#include <cmath>
#include <stdint.h>
#include <cmath>
#include <climits>
#include <cassert>
#include <cfloat>
#include <stdexcept>
#include <stdio.h>

namespace flopoco
{
namespace random
{

template<class T>
struct result_type
{
	int expMin;
	int expMax;
	int fracWidth;
	T fracMin;
	T fracMax;
	T valueMin;
	T valueMax;
	T absErrorMin;
	T absErrorMax;
	T relErrorMin;
	T relErrorMax;
	
	// Default constructor is for first element in chain, where there is no existing answer
	result_type(int _fracWidth)
		: expMax(INT_MIN)
		, expMin(INT_MAX)
		, fracWidth(_fracWidth)	// doesn't include the implicit bit
		, fracMin(DBL_MAX)
		, fracMax(DBL_MIN)
		, valueMin(DBL_MAX)
		, valueMax(DBL_MIN)
		, absErrorMin(DBL_MAX)
		, absErrorMax(DBL_MIN)
		, relErrorMin(DBL_MAX)
		, relErrorMax(DBL_MIN)
	{
		if(_fracWidth >= std::numeric_limits<T>::digits){
			throw std::invalid_argument("result_type - Type can't represent this many fractional bits.");
		}
	}
		
	// Width as bits, does not include implicit bit
	unsigned Width() const
	{
		if(expMax <= expMin)
			return fracWidth;
		unsigned expWidth=(unsigned)ceil(log(expMax-expMin+1)/log(2.0));
		return expWidth+fracWidth;
	}
	
	T FromBits(uint64_t bits) const
	{
		assert(bits < (1ull<<Width()));
		uint64_t f_bits=bits&((1ull<<fracWidth)-1);
		int e=(bits>>fracWidth)+expMin;
		T f=ldexp(T(f_bits), -fracWidth-1)+0.5;
		return ldexp(f, e);
	}
	
	uint64_t ToBits(const T &x) const
	{
		assert(x>0);
		int e;
		T of=frexp(x, &e);
		T f=of;
		assert((0<=f)&&(f<1));
		f=ldexp(f,fracWidth+1);
		assert(round(f)==f);
		assert((expMin <=e) && (e<=expMax));
		uint64_t f_bits=(uint64_t)f;
		assert(f_bits & (1ull<<fracWidth));
		f_bits &= ~(1ull<<fracWidth);
		assert(f_bits<(1ull<<fracWidth));
		uint64_t e_bits=e-expMin;
		fprintf(stderr, "x=%lg, f=%lg, e=%d, f_bits=0x%llx, e_bits=0x%llx\n", x, of, e, f_bits, e_bits);
		uint64_t res= (e_bits<<fracWidth) | f_bits;
		assert(res<(1ull<<Width()));
		return res;
	}
	
	T Round(const T &x) const
	{
		if(x<=0)
			throw std::invalid_argument("x must be strictly positive.");
		int e;
		T f=frexp(x, &e);
		return ldexp(round(ldexp(f,fracWidth+1)), e-fracWidth-1);
	}
	
	void Add(const T &x, const T &correct)
	{
		assert(x==Round(x));
		
		int e;
		T f=frexp(x, &e);
		
		expMin=std::min(expMin, e);
		expMax=std::max(expMax, e);
		fracMin=std::min(fracMin, f);
		fracMax=std::max(fracMax, f);
		valueMin=std::min(valueMin, x);
		valueMax=std::max(valueMax, x);
		
		if(x!=correct){
			T absError=x-correct;
			T relError=absError/correct;
			absErrorMin=std::min(absErrorMin, absError);
			absErrorMax=std::max(absErrorMax, absError);
			relErrorMin=std::min(relErrorMin, relError);
			relErrorMax=std::max(relErrorMax, relError);
		}
	}
	
	void Add(const T &x)
	{
		Add(x,x);
	}
	
	T RoundAndAdd(const T &x)
	{
		T got=Round(x);
		Add(got, x);
		return got;
	}
		
	
};

}; // random
}; // flopoco

#endif
