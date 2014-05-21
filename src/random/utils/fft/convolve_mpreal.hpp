#ifndef flopoco_random_utils_convolve_mpreal_hpp
#define flopoco_random_utils_convolve_mpreal_hpp

#include <cmath>
#include <complex>
#include <vector>

#include "random/utils/fft/fft.hpp"
#include "random/utils/fft/fft_mpreal.hpp"

#include "random/utils/fft/convolve.hpp"

namespace flopoco
{
namespace random
{
	inline std::vector<mpfr::mpreal> convolve_mpreal(const std::vector<mpfr::mpreal> &a, const std::vector<mpfr::mpreal> &b, unsigned wPrec)
	{
		unsigned n=a.size()+b.size()-1;
		unsigned nn=detail::NextBinaryPower(n);
		
		if(wPrec==0){
			wPrec=mpfr_get_default_prec();
			for(unsigned i=0;i<a.size();i++){
				wPrec=std::max(wPrec, (unsigned)a[i].get_prec());
			}
			for(unsigned i=0;i<b.size();i++){
				wPrec=std::max(wPrec, (unsigned)b[i].get_prec());
			}
		}
		
		std::vector<mpfr::mpreal> aa(2*nn, mpfr::mpreal(0,wPrec)), bb(2*nn, mpfr::mpreal(0,wPrec));
		for(unsigned i=0;i<a.size();i++){
			mpfr_set(aa[2*i].mpfr_ptr(), a[i].mpfr_srcptr(), MPFR_RNDN);
		}
		for(unsigned i=0;i<b.size();i++){
			mpfr_set(bb[2*i].mpfr_ptr(), b[i].mpfr_srcptr(), MPFR_RNDN);
		}
		
		fft_radix2_complex(&aa[0], nn, wPrec);
		fft_radix2_complex(&bb[0], nn, wPrec);
		
		for(unsigned i=0;i<nn;i++){
			std::complex<mpfr::mpreal> av(aa[2*i], aa[2*i+1]);
			std::complex<mpfr::mpreal> bv(bb[2*i], bb[2*i+1]);
			av=av*bv;
			aa[2*i]=real(av);
			aa[2*i+1]=imag(av);
		}
		
		ifft_radix2_complex(&aa[0], nn, wPrec);
		
		std::vector<mpfr::mpreal> res(n);
		for(unsigned i=0;i<n;i++){
			std::swap(res[i], aa[2*i]);
		}
		
		return res;
	}
	
	inline std::complex<mpfr::mpreal> mini_pow(std::complex<mpfr::mpreal> x, int k, unsigned wPrec)
	{
		assert(0);	// This needs to be checked, never been tested or run yet
		
		std::complex<mpfr::mpreal> acc(mpfr::mpreal(1, wPrec));
		while(k>0){
			if(k&1){
				acc=acc*x;
			}
			k=k/2;
			if(k>0){
				x=x*x;
			}
		}
		return acc;
	}

	inline std::vector<mpfr::mpreal> self_convolve_mpreal(const std::vector<mpfr::mpreal> &a, unsigned k, unsigned wPrec)
	{
		if(k==0)
			throw std::string("self_convolve_mpreal - Self convolution of degree 0 is probably not intended.");
		if(k==1)
			return a;
		
		if(wPrec==0){
			wPrec=mpfr_get_default_prec();
			for(unsigned i=0;i<a.size();i++){
				wPrec=std::max(wPrec, (unsigned)a[i].get_prec());
			}
		}
		
		unsigned n=a.size()*k-(k-1);
		unsigned nn=detail::NextBinaryPower(n);
		
		std::vector<mpfr::mpreal> aa(2*nn, mpfr::mpreal(0,wPrec));
		for(unsigned i=0;i<a.size();i++){
			mpfr_set(aa[2*i].mpfr_ptr(), a[i].mpfr_srcptr(), MPFR_RNDN);
		}
		
		fft_radix2_complex(&aa[0], nn, wPrec);
		
		std::complex<mpfr::mpreal> one(mpfr::mpreal(1, wPrec));
		for(unsigned i=0;i<nn;i++){
			std::complex<mpfr::mpreal> av(aa[2*i], aa[2*i+1]);
			//av=std::pow(av, k);	// No longer works in C++11...
			av=mini_pow(av, k, wPrec);
			aa[2*i]=real(av);
			aa[2*i+1]=imag(av);
		}
		
		ifft_radix2_complex(&aa[0], nn, wPrec);
		
		std::vector<mpfr::mpreal> res(n);
		for(unsigned i=0;i<n;i++){
			std::swap(res[i], aa[2*i]);
		}
		
		return res;
	}
	
	template<>
	inline std::vector<mpfr::mpreal> convolve(const std::vector<mpfr::mpreal> &a, const std::vector<mpfr::mpreal> &b)
	{
		return convolve_mpreal(a,b,0);
	}
	
	template<>
	inline std::vector<mpfr::mpreal> self_convolve(const std::vector<mpfr::mpreal> &a, unsigned k)
	{
		return self_convolve_mpreal(a,k,0);
	}
	
	
};
};

#endif
