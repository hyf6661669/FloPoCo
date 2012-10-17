#include "comparable_float_type.hpp"

#include "mpfr_vec.hpp"

#define BOOST_TEST_MODULE ComparableFloatType
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(TestCompile)
{
}

BOOST_AUTO_TEST_CASE(TestCross)
{
	gmp_randstate_t r;
	gmp_randinit_default(r);
	
	srand48(1);
	
	for(int wE=3;wE<=8;wE++){
		for(int wF=3;wF<=16;wF++){
			fprintf(stderr, "wE=%d,wF=%d\n", wE, wF);
			
			flopoco::random::ComparableFloatType codec(wE, wF);
	
			int n=1000;
			
			flopoco::random::MPFRVec orig(n, wF+1), tmp(1, wF+1);
			std::vector<mpz_class> encoded(n);
			
			mpz_class wFpow(mpz_class(1)<<wF);
			
			mpz_class fraction;
			
			mpfr_set_zero(orig[0], +1);
			mpfr_set_zero(orig[1], -1);
			for(int i=2;i<n;i++){
				int e=(int)floor(drand48()*(1<<wE))-(1<<(wE-1));
				if((i%7)==0)
					e=-(1<<(wE-1));
				if((i%13)==0)
					e=(1<<(wE-1))-1;
				
				mpz_urandomb(fraction.get_mpz_t(), r, wF);
				fraction=fraction+wFpow;
				if((i%17)==0)
					fraction=wFpow;
				if((i%19)==0)
					fraction=2*wFpow-1;
				
				mpfr_set_z(orig[i], fraction.get_mpz_t(), MPFR_RNDN);
				mpfr_mul_2si(orig[i], orig[i], e-wF, MPFR_RNDN);
				
				if(drand48()>0.5)
					mpfr_setsign(orig[i], orig[i], -1, MPFR_RNDN);
			}	
			
			for(int i=0;i<n;i++){
				encoded[i]=codec.ToBits(orig[i]);
				codec.FromBits(tmp[0], encoded[i]);
				
				assert(mpfr_equal_p(tmp[0], orig[i]));
			}
			
			for(int i=0;i<n;i++){
				for(int j=0;j<n;j++){
					//mpfr_fprintf(stderr, "Left=%Rg (%Zx), Right=%Rg (%Zx)\n", orig[i], encoded[i].get_mpz_t(), orig[j], encoded[j].get_mpz_t());
					
					if(mpfr_zero_p(orig[i]) && mpfr_zero_p(orig[j])){
						assert((mpfr_signbit(orig[i]) > mpfr_signbit(orig[j])) == (encoded[i] < encoded[j]));
					}else{
						assert(mpfr_less_p(orig[i], orig[j]) == (encoded[i]<encoded[j]));
						assert(mpfr_greater_p(orig[i], orig[j]) == (encoded[i]>encoded[j]));
						assert(mpfr_equal_p(orig[i], orig[j]) == (encoded[i]==encoded[j]));
					}
				}
			}
		}
	}
	
	gmp_randclear(r);
}
