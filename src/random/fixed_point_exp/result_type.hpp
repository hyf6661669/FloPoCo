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
#include <iostream>
#include <sstream>

#include <gmpxx.h>

namespace flopoco
{
namespace random
{
	
mpz_class to_mpz_class(uint64_t x)
{
	uint32_t high=x>>32;
	if(high==0){
		return mpz_class((uint32_t)(x&0xFFFFFFFFu));
	}else{
		mpz_class tmp=mpz_class(high);
		tmp=tmp<<32;
		tmp=tmp+uint32_t(x&0xFFFFFFFF);
		return tmp;
	}
}

uint64_t to_uint64_t(const mpz_class &x)
{
	if(x<0)
		throw std::invalid_argument("to_uint64_t(mpz_class) - Value is negative.");
	if(mpz_sizeinbase (x.get_mpz_t(),2)>64)
		throw std::invalid_argument("to_uint64_t(mpz_class) - Value is larger than 64-bits.");
	uint64_t dst=0;
	mpz_export (&dst, 0 /*size_t *COUNTP*/, -1 /*int ORDER*/,
          sizeof(uint64_t) /*size_t SIZE*/, 0 /*int ENDIAN*/, 0 /*size_t NAILS*/, x.get_mpz_t() /*mpz_t OP*/);
	return dst;
}

template<class T>
struct result_type
{
protected:
	T random() const
	{ return drand48(); }
	
public:
	
	int expMax;
	int expMin;
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
	
	unsigned ExpWidth() const
	{
		if(expMax<expMin)
			return fracWidth;
		return (unsigned)ceil(log(expMax-expMin+1)/log(2.0));
	}
	
	unsigned FracWidth() const
	{
		if(expMax < expMin)
			return 0;
		return fracWidth;
	}
	
	// Width as bits, does not include implicit bit
	unsigned Width() const
	{ return FracWidth()+ExpWidth(); }
	
	T TypeMin() const
	{
		if(expMax<expMin)
			return 1;
		return ldexp(0.5, expMin);
	}
	
	T TypeMax() const
	{
		if(expMax<expMin)
			return 1;
		return ldexp(1-pow(2.0,-(fracWidth+1)), expMax);
	}
	
	T RangeMin() const
	{
		if(valueMax<valueMin)
			return TypeMin();
		return valueMin;
	}
	
	T RangeMax() const
	{
		if(valueMax<valueMin)
			return TypeMax();
		return valueMax;
	}
	
	T FracMin() const
	{
		if(fracMax<fracMin)
			return 0.5;
		return fracMin;
	}
	
	T FracMax() const
	{
		if(fracMax < fracMin)
			return 1-pow(2.0,-(fracWidth+1));
		return fracMax;
	}
	
	//! Returns the number of bits which are non-zero, based on FracMin()..FracMax()
	/*! This is only able to take into account values which end up as 0.1000000XXXXX,
		doesn't do anything with things like  0.10000011XXX, so it would count the (non-leading) constant ones*/
	unsigned FracWidthNonZero() const
	{
		if(fracMax < fracMin)
			return FracWidth();
		if(fracMax==0.5)
			return 0;
		double tmp=fracMax-0.5;
		int res=FracWidth();
		while(tmp<0.25){
			tmp=tmp*2;
			res--;
			assert(res <(int)FracWidth());
		}
		return res;
	}
	
	bool IsInRange(const T &x) const
	{
		if(expMax < expMin){
			if(x==1.0)
				return true;
		}
		int e;
		if(valueMin <= valueMax){
			if((x<valueMin) || (valueMax <x))
				return false;
		}
		T f=frexp(x, &e);
		if((e<expMin) || (expMax<e))
			return false;
		if(fracMin<=fracMax){
			if((f<fracMin) || (fracMax<f))
				return false;
		}
		return true;
	}
	
	T RandomElement() const
	{
		if((expMax<expMin))
			return 1.0;
		while(true){
			int e=(int)floor(expMin+random()*(expMax+1-expMin));
			T frac=0.5+random()/2;
			T v=ldexp(frac, e);
			fprintf(stderr, "e=%d, frac=%lf, v=%lg\n", e, frac, v);
			v=Round(v);
			if(IsInRange(v)){
				return v;
			}
		}
	}
	
	T FromBits(uint64_t bits) const
	{
		if(expMax<expMin){
			assert(bits==0);
			return 1.0;
		}
		assert(bits < (1ull<<Width()));
		uint64_t f_bits=bits&((1ull<<fracWidth)-1);
		int e=(bits>>fracWidth)+expMin;
		T f=ldexp(T(f_bits), -fracWidth-1)+0.5;
		return ldexp(f, e);
	}
	
	T FromMpz(const mpz_class &x) const
	{ return FromBits(to_uint64_t(x)); }
	
	uint64_t ToBits(const T &x) const
	{
		if(expMax < expMin){
			std::cerr<<"   ToBits, type="<<*this<<"\n";
			assert(x==1.0);
			return 0;
		}
		
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
		//fprintf(stderr, "x=%lg, f=%lg, e=%d, f_bits=0x%llx, e_bits=0x%llx\n", x, of, e, f_bits, e_bits);
		uint64_t res= (e_bits<<fracWidth) | f_bits;
		assert(res<(1ull<<Width()));
		return res;
	}
	
	mpz_class ToMpz(const T &x) const
	{ return to_mpz_class(ToBits(x)); }
	
	T Next(const T &x) const
	{
		// This is a hack that works well for strictly positive floating-point
		return FromBits(ToBits(x)+1);
	}
	
	T Round(const T &x) const
	{
		if(x<=0)
			throw std::invalid_argument("x must be strictly positive.");
		int e;
		T f=frexp(x, &e);
		return ldexp(round(ldexp(f,fracWidth+1)), e-fracWidth-1);
	}
	
	T RoundHalfDown(const T &x) const
	{
		if(x<=0)
			throw std::invalid_argument("x must be strictly positive.");
		int e;
		T f=frexp(x, &e);
		return ldexp(ceil(ldexp(f,fracWidth+1)-0.5), e-fracWidth-1);
	}
	
	T Floor(const T &x) const
	{
		if(x<=0)
			throw std::invalid_argument("x must be strictly positive.");
		int e;
		T f=frexp(x, &e);
		return ldexp(floor(ldexp(f,fracWidth+1)), e-fracWidth-1);
	}
	
	T Ceil(const T &x) const
	{
		if(x<=0)
			throw std::invalid_argument("x must be strictly positive.");
		int e;
		T f=frexp(x, &e);
		return ldexp(ceil(ldexp(f,fracWidth+1)), e-fracWidth-1);
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

template<class T>
std::ostream &operator<<(std::ostream &dst, const result_type<T> &rt)
{
	if(rt.expMax<rt.expMin){
		dst<<"PosFloat[e=null;w="<<rt.fracWidth<<"]";
	}else{
		if(rt.valueMax < rt.valueMin){
			dst<<"PosFloat[e="<<rt.expMin<<".."<<rt.expMax<<";w="<<rt.fracWidth<<"]";
		}else{
			dst<<"PosFloat[e="<<rt.expMin<<".."<<rt.expMax<<";w="<<rt.fracWidth<<";value=("<<rt.valueMin<<"..."<<rt.valueMax<<")]";
		}
	}
	return dst;
}

}; // random
}; // flopoco

#endif
