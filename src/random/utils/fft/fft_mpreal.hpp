#ifndef flopoco_random_utils_fft_mpreal_hpp
#define flopoco_random_utils_fft_mpreal_hpp

#include <cassert>

#include "mpreal.h"

#include "random/utils/fft/fft.hpp"

namespace flopoco
{
namespace random
{
	namespace detail{
		
		// FFT based on four1 from numerical recipes. Input is as a zero-based array (i.e. normal C, not NumRec style)
		// TODO : What is the license on numerical recipes...
		/*! \param wPrec is the working precision we want for calculations
			\note The trig recurrence from numerical recipes is supposed to be accurate to O(sqrt(n)). Putting
			a constant on that lets say the error is about (8*eps)*sqrt(n), given the number of operations
			in the recurrence. So if we do the recurrence with wt bits, then eps=2^-wt, and we have
			error ~ sqrt(n)*8*2^-wt. To meet a given working precision wp, then we want
			sqrt(n)*2^(3-wt) < 2^-wp ->
			sqrt(n) < 2^-wp / 2^(2-wt) = 2^(-wp+wt-3)
			n < 2^[2(wt-wp-3)]
			Or:
			2^-wt < 2^(-wp-3) / sqrt(n)
			-wt < log2( 2^(-wp-3) / sqrt(n) )
			wt > log2( 2^(-wp-3) / sqrt(n) )
		
			Technically we could derive a balance between the cost for the recurrence vs sin/cos
			at a particular level, but it probably makes sense to just allow n to grow up to some maximum,
			then always do a fresh trig call at that point. This error analysis is a bit hokey
			anyway, as the target precision (wPrec) is quite poorly defined :)
		*/
		inline void fft_radix2(mpfr::mpreal *data, unsigned nn, bool inverse, unsigned wPrec)
		{
			if(!IsBinaryPower(nn))
				throw std::string("four1 - FFT Input data size must be power of 2.");
			
			// We've been told to use wPrec bits, and that's what we're going to do, damnit!
			for(unsigned i=0;i<2*nn;i++){
				mpfr_prec_round(data[i].mpfr_ptr(), wPrec, MPFR_RNDN);
			}
			
			const unsigned tMaxLen=4096;	// Maximum distance before refreshing sincos. Presumably 4096 recurrences cost less than one sin/cos?
			const unsigned tSafety=8;	// Safety bits to tack on, as we have some precalculations
			// -log2( pow(2,-wPrec-3) / sqrt(2^log2(tMaxLen)) ) + tSafety
			// -log2( pow(2,-wPrec-3) / (2^log2(tMaxLen)/2) ) + tSafety
			// -log2( pow(2,-wPrec-3 - log2(tMaxLen)/2 ) + tSafety
			//  wPrec+3+log2(tMaxLen)/2 + tSafety
			unsigned tPrec=wPrec +log2(tMaxLen)/2 + tSafety;
			// tMaxLen = 4096 -> wPrec+3+5+8 = wPrec+ 16
			// So it's not costing us a lot in extra precision
			
			mpfr::mpreal pi=mpfr::const_pi(tPrec);
			
			unsigned n, mmax, m, j, istep, i;
			
			n=nn<<1;
			j=1;
			for(i=1;i<n;i+=2){
				if(j>i){
					assert(j < 2*nn);
					assert(i < 2*nn);
					std::swap(data[j-1], data[i-1]);
					std::swap(data[j], data[i]);
				}
				m=n>>1;
				while(m>=2 && j>m){
					j-=m;
					m>>=1;
				}
				j+=m;
			}
			
			mpfr::mpreal twiddle_r(0,tPrec), twiddle_i(0,tPrec), twiddle_r_old(0,tPrec), delta_theta(0,tPrec);
			mpfr::mpreal alpha(0,tPrec), beta(0,tPrec);
			
			mpfr::mpreal twiddle_r_wPrec(0,wPrec), twiddle_i_wPrec(0,wPrec);
			
			bool extraWorking=false;	// This doesn't really help at all
			mpfr::mpreal temp_i(0,extraWorking?2*wPrec:wPrec), temp_r(0, extraWorking?2*wPrec:wPrec);
			
			mmax=2;
			while(n>mmax){
				istep=mmax<<1;
				
				delta_theta=(pi<<1)/mmax;
				if(inverse)
					delta_theta=-delta_theta;
				alpha= sqr(sin(delta_theta>>1))<<1;
				beta =sin(delta_theta);
				unsigned trig_steps=0;
				
				for(m=0;m<mmax;m+=2){
					if(trig_steps==0){
						twiddle_r=1.0;	// cos(0)==cos(trig_steps*delta_theta)
						twiddle_i=0.0;	// sin(0)==sin(trig_steps*delta_theta)
					}else if((trig_steps%tMaxLen)==0){
						twiddle_r=cos(delta_theta*trig_steps);
						twiddle_i=sin(delta_theta*trig_steps);
					}else{
						twiddle_r_old=twiddle_r;
						twiddle_r = twiddle_r - (alpha*twiddle_r + beta*twiddle_i);
						twiddle_i = twiddle_i - (alpha*twiddle_i - beta*twiddle_r_old);
					}
					trig_steps++;
					
					// Lose the extra precision for the actual calculations
					mpfr_set(twiddle_r_wPrec.mpfr_ptr(), twiddle_r.mpfr_ptr(), MPFR_RNDN);
					mpfr_set(twiddle_i_wPrec.mpfr_ptr(), twiddle_i.mpfr_ptr(), MPFR_RNDN);
					
					for(i=m;i<n;i+=istep){
						j=i+mmax;
						assert(j+1 < 2*nn);
						assert(i+1 < 2*nn);
						
						if(extraWorking){
							mpfr_mul(temp_r.mpfr_ptr(), twiddle_i_wPrec.mpfr_ptr(), data[j+1].mpfr_ptr(), MPFR_RNDN);
							mpfr_fms(temp_r.mpfr_ptr(), twiddle_r_wPrec.mpfr_ptr(), data[j].mpfr_ptr(), temp_r.mpfr_ptr(), MPFR_RNDN);
							
							mpfr_mul(temp_i.mpfr_ptr(), twiddle_r_wPrec.mpfr_ptr(), data[j+1].mpfr_ptr(), MPFR_RNDN);
							mpfr_fma(temp_i.mpfr_ptr(), twiddle_i_wPrec.mpfr_ptr(), data[j].mpfr_ptr(), temp_i.mpfr_ptr(), MPFR_RNDN);
							
							mpfr_sub(data[j].mpfr_ptr(), data[i].mpfr_ptr(), temp_r.mpfr_ptr(), MPFR_RNDN);
							mpfr_sub(data[j+1].mpfr_ptr(), data[i+1].mpfr_ptr(), temp_i.mpfr_ptr(), MPFR_RNDN);
							
							mpfr_add(data[i].mpfr_ptr(), data[i].mpfr_ptr(), temp_r.mpfr_ptr(), MPFR_RNDN);
							mpfr_add(data[i+1].mpfr_ptr(), data[i+1].mpfr_ptr(), temp_i.mpfr_ptr(), MPFR_RNDN);
						}else{
							temp_r=twiddle_r_wPrec*data[j]    	- 	twiddle_i_wPrec*data[j+1];
							temp_i=twiddle_r_wPrec*data[j+1]	+	twiddle_i_wPrec*data[j];
							
							data[j]=data[i]-temp_r;
							data[j+1]=data[i+1]-temp_i;
							data[i]+=temp_r;
							data[i+1]+=temp_i;
						}
					}
					
					
				}
				mmax=istep;
			}
			
			if(inverse){
				mpfr::mpreal scale(1.0, wPrec);
				scale=1.0;
				scale=scale/nn;
				// Divide everything by nn. Given nn is a binary power, this is just shifting
				int places=(int)round(log2(nn));
				for(unsigned i=0;i<2*nn;i++){
					data[i] >>=places;
				}
			}
		}
	};
	
	inline void fft_radix2_complex(mpfr::mpreal *data, unsigned nn, unsigned wPrec)
	{ detail::fft_radix2(data, nn, false, wPrec); }
	
	inline void ifft_radix2_complex(mpfr::mpreal *data, unsigned nn, unsigned wPrec)
	{ detail::fft_radix2(data, nn, true, wPrec); }
};
};

#endif
