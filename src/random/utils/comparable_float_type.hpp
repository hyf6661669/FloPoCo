#ifndef flopoco_random_utils_comparable_float_type_hpp
#define flopoco_random_utils_comparable_float_type_hpp

#include <mpfr.h>
#include <gmpxx.h>
#include <assert.h>
#include <math.h>

#include "Operator.hpp"

namespace flopoco
{
namespace random
{
	
class ComparableFloatEncoder;

/* An ordered float is a format where floating-point numbers can
	be directly compared. The format is a prefix|exponent|fraction format
	
	10|100|111	0.875*2^0
	10|100|000	0.5*2^0
	10|011|111	0.875*2^-1
	10|000|000	0.5*2^-4
	01|000|000	0
	00|111|111	-0.5*2^-4
	00|100|000	-0.875*2^-1
	00|011|111	0.5*2^0
	00|011|000	-0.875*2^0
	
	Prefixes are:
	11 : positive regular
	10 : positive zero
	01 : negative zero
	00 : negative regular
	
	exponent field is  "e+2^(wE-1)" if positive, or "2^(wE-1)-e-1" if negative
	fraction field is "f" if positive, or "2^wF-f" if negative
*/
class ComparableFloatType
{
public:
	const int wE, wF;
private:
	double DecodePrefix(unsigned prefix) const
	{
		switch(prefix){
		case 3: return +1.0;
		case 2: return +0.0;
		case 1: return -0.0;
		case 0: return -1.0;
		default: throw std::invalid_argument("ComparableFloatType - invalid prefix.");
		}
	}
	
	unsigned EncodePrefix(mpfr_t x) const
	{
		if(mpfr_zero_p(x))
			return mpfr_signbit(x)==0 ? 2 : 1;		
		if(mpfr_regular_p(x))
			return mpfr_signbit(x)==0 ? 3 : 0;
		mpfr_fprintf(stderr, "  x=%Rg\n", x);
		throw std::invalid_argument("ComparableFloatType - cannot handle non-regular numbers.");
	}
	
	unsigned EncodeExponent(mpfr_t x) const
	{
		if(mpfr_zero_p(x))
			return 0;
		if(!mpfr_regular_p(x)){
			mpfr_fprintf(stderr, "  x=%Rg\n", x);
			throw std::invalid_argument("ComparableFloatType - cannot handle non-regular numbers.");
		}
		int e=mpfr_get_exp(x)-1;
		int bias=(1<<(wE-1));
		
		if(e < -bias){
		    std::stringstream acc;
		    acc<<"ComparableFloatType - exponent too small, got="<<e<<", min="<<-bias;
		    throw std::invalid_argument(acc.str());
		}
		if(e >= bias)
			throw std::invalid_argument("ComparableFloatType - exponent too large.");
		if(mpfr_sgn(x)>=0)
			e=e+bias;
		else
			e=bias-e-1;
		assert(e>=0);
		assert(e<(1<<wE));
		return e;
	}
	
	mpz_class EncodeFraction(mpfr_t x) const
	{		
		if(mpfr_zero_p(x))
			return mpz_class(0);
		
		mpfr_t tmp;
		mpfr_init2(tmp, wF+1);
		
		if(mpfr_set(tmp, x, MPFR_RNDN)!=0){
			mpfr_clear(tmp);
			throw std::invalid_argument("ComparableFloatType - mpfr source cannot be represented in target type.");
		}
		
		mpz_class wFpow=(mpz_class(1)<<wF);
		
		mpfr_mul_2si(tmp, tmp, -mpfr_get_exp(tmp)+wF+1, MPFR_RNDN);
	
		mpz_class res;
		if(mpfr_get_z(res.get_mpz_t(), tmp, MPFR_RNDN)!=0){
			//mpfr_fprintf(stderr, "tmp=%Rg*2^%d, res=%Zd,  wFpow=%Zu\n", tmp, mpfr_get_exp(x)-wF-1, res.get_mpz_t(), wFpow.get_mpz_t());
			mpfr_clear(tmp);
			throw std::logic_error("ComparableFloatType - problem while rounding (internal error).");
		}
		if(res<0)
			res=-res;
		
		res=res-wFpow;
		
		assert(res >= 0);
		assert(res < wFpow);
		
		if(mpfr_sgn(x)<0)
			res=wFpow-res-1;
		
		//mpfr_fprintf(stderr, "res=%Zd\n", res.get_mpz_t());
		assert(res >= 0);
		assert(res < wFpow);
		
		return res;
	}
	
public:
	ComparableFloatType(int _wE, int _wF)
		: wE(_wE)
		, wF(_wF)
	{}
		
	int Width() const
	{ return 2+wE+wF; }
	
	void Round(mpfr_t y, mpfr_t x)
	{
		mpfr_t tmp;
		mpfr_init2(tmp, wF+1);
		
		mpfr_set(tmp, x, MPFR_RNDN);
		mpfr_set(y, tmp, MPFR_RNDN);
		
		mpfr_clear(tmp);
	}

	mpz_class ToBits(mpfr_t x, bool doRound=false)
	{
		mpfr_t tmp;
		mpfr_init2(tmp, wF+1);
		try{
			//mpfr_fprintf(stderr, "Encode(%Rg)\n", x);
			
			if(doRound){
				mpfr_set(tmp, x, MPFR_RNDN);
			}else{
				if(mpfr_set(tmp, x, MPFR_RNDN)!=0)
					throw std::invalid_argument("ComparableFloatType::Encode - Given value cannot be exactly represented in this format.");
			}
			
			mpz_class res=EncodePrefix(tmp);
			res=(res<<wE)+EncodeExponent(tmp);
			res=(res<<wF)+EncodeFraction(tmp);
	
			mpfr_clear(tmp);
			
			return res;
		}catch(...){
			mpfr_clear(tmp);
			throw;
		}
	}
	
	void FromBits(mpfr_t res, mpz_class x)
	{
		mpfr_t tmp;
		mpfr_init2(tmp, wF+1);
		try{
			mpz_class tmp=x;
			
			mpz_class fraction;
			mpz_tdiv_r_2exp(fraction.get_mpz_t(), tmp.get_mpz_t(), wF);
			tmp=tmp>>wF;
			
			mpz_class exponent;
			mpz_tdiv_r_2exp(exponent.get_mpz_t(), tmp.get_mpz_t(), wE);
			tmp=tmp>>wE;
			
			double prefix=DecodePrefix(tmp.get_ui());
			
			//mpfr_fprintf(stderr, "Decode: prefix=%g, exponent=%Zd, fraction=%Zd\n", prefix, exponent.get_mpz_t(), fraction.get_mpz_t());
			
			int e=exponent.get_ui();
			
			int bias=1<<(wE-1);
			mpz_class wFpow=(mpz_class(1)<<wF);
			
			if(prefix==0){
				mpfr_set_zero(res, prefix==-0.0 ? -1 : +1);
			}else{
				if(prefix<0){
					// e'=bias-e-1;
					// e'-bias+1=-e
					// bias-e'-1=e
					e=bias-e-1;
					fraction=wFpow-fraction-1;
				}else{
					e=e-bias;
				}
				assert(-bias<=e);
				assert(e<bias);
				assert(fraction>=0);
				assert(fraction<wFpow);
				
				fraction=fraction+wFpow;
				
				if(0!=mpfr_set_z(res, fraction.get_mpz_t(), MPFR_RNDN))
					throw std::invalid_argument("ComparableFloatType::Decode - destination is not precise enough to hold type.");
				
				mpfr_mul_2si(res, res, e-wF, MPFR_RNDN);
				if(prefix<0){
					mpfr_neg(res, res, MPFR_RNDN);
				}
			}
		}catch(...){
			mpfr_clear(tmp);
			throw;
		}
	}
	
	//! Make an encoded which maps from a flopoco float to a comparable float
	ComparableFloatEncoder *MakeEncoder(Target *target, map<string, double> inputDelays = emptyDelayMap);
};

class ComparableFloatEncoder
	: public Operator
{
private:
	friend class ComparableFloatType;

	ComparableFloatType type;
	const int &wE, &wF;

	ComparableFloatEncoder(Target *target, ComparableFloatType type, map<string, double> inputDelays = emptyDelayMap);

	void addTestCase(TestCaseList* tcl, mpfr_t val);
public:
	~ComparableFloatEncoder()
	{}
		
	void emulate(TestCase * tc);
	void buildStandardTestCases(TestCaseList* tcl);
	TestCase* buildRandomTestCase(int i);
};

}; // random
}; // flopoco

#endif
