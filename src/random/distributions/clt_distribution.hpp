#ifndef random_base_clt_distribution_hpp
#define random_base_clt_distribution_hpp

#include "table_distribution.hpp"

#include <vector>
#include <algorithm>
#include <numeric>

#include "gsl/gsl_fft_real.h"
#include "gsl/gsl_fft_halfcomplex.h"

namespace flopoco
{
namespace random
{
	template<class T>
	std::vector<T> convolve(const std::vector<T> &a, const std::vector<T> &b);
	
	template<>
	std::vector<double> convolve(const std::vector<double> &a, const std::vector<double> &b)
	{
		unsigned n=a.size()+b.size()-1;
		n=pow(2.0, ceil(log(n) / log(2.0)));
		
		// This is a bit of a hack, it should really defer to the arbitrary precision version in moptl2
		std::vector<double> aa(n, 0.0);
		std::copy(a.begin(), a.end(), aa.begin());
		if(0!=gsl_fft_real_radix2_transform(&aa[0], 1, n))
			throw std::runtime_error("Error from gsl_fft_real_radix2_transform");
		
		std::vector<double> bb(n, 0.0);
		std::copy(b.begin(), b.end(), bb.begin());
		if(0!=gsl_fft_real_radix2_transform(&bb[0], 1, n))
			throw std::runtime_error("Error from gsl_fft_real_radix2_transform");
		
		aa[0]=aa[0]*bb[0];
		for(unsigned i=1;i<n/2;i++){
			std::complex<double> va(aa[i], aa[n-i]);
			std::complex<double> vb(bb[i], bb[n-i]);
			std::complex<double> vr=va*vb;
			aa[i]=real(vr);
			aa[n-i]=imag(vr);
		}
		aa[n/2]=aa[n/2]*bb[n/2];
		
		bb.clear();
		
		if(0!=gsl_fft_halfcomplex_radix2_inverse(&aa[0], 1, n))
			throw std::runtime_error("Error from gsl_fft_halfcomplex_radix2_inverse");
		
		aa.resize(a.size()+b.size()-1);
		return aa;
	}
	
	template<class T>
	std::vector<T> self_convolve(const std::vector<T> &a, unsigned k);
	

	template<>
	std::vector<double> self_convolve(const std::vector<double> &a, unsigned k)
	{
		unsigned n=k*a.size()-1;
		n=pow(2.0, ceil(log(n) / log(2.0)));
		
		// This is a bit of a hack, it should really defer to the arbitrary precision version in moptl2
		std::vector<double> res(n, 0.0);
		std::copy(a.begin(), a.end(), res.begin());
		
		if(0!=gsl_fft_real_radix2_transform(&res[0], 1, n))
			throw std::runtime_error("Error from gsl_fft_real_radix2_transform");
		
		res[0]=pow(res[0], (int)k);
		for(unsigned i=1;i<n/2;i++){
			std::complex<double> v(res[i], res[n-i]);
			v=pow(v, (int)k);
			res[i]=real(v);
			res[n-i]=imag(v);
		}
		res[n/2]=pow(res[n/2], (int)k);
		
		if(0!=gsl_fft_halfcomplex_radix2_inverse(&res[0], 1, n))
			throw std::runtime_error("Error from gsl_fft_halfcomplex_radix2_inverse");
		
		res.resize(a.size()*k-(k-1));
		return res;
	}
	
	template<class T>
	inline typename TableDistribution<T>::TypePtr MakeCLTDistribution(
		unsigned w,	//! Number of bits used in each base component
		unsigned k,	//! Number of accumulation steps, must be even and >=2
		unsigned shift=0	//! Post-shift to apply: shift=0 -> [0,1),  shift=w -> [0,2^w)
	){
		T scale;
		scale=1;
		scale=scale/pow(2.0,w);
		std::vector<T> probs(1<<w, scale);
		
		probs=self_convolve(probs, k);
		T acc=0.0;
		for(unsigned i=0;i<probs.size()/2;i++){
			if(probs[i]<0){
				if(probs[i] < -10e-13)
					throw std::runtime_error("Received negative probability less than 10^-13.");
				probs[i]=0;
			}
			acc += 2*probs[i];
		}
		if(probs.size()%2)
			acc += probs[probs.size()/2];
				
		if(fabs(acc-1) > pow(2.0,-10)){
			fprintf(stderr, "acc=%lg\n",acc);
			throw std::runtime_error("Sum of PDF after convolution differs from 1.0 by more than 10^-10.");
		}
		
		scale=1.0/acc;
		for(unsigned i=0;i<=probs.size()/2;i++){
			T p=probs[i] * scale;
			probs[i]=p;
			probs[probs.size()-i-1]=p;
		}
	
		std::vector<std::pair<T,T> > table(probs.size());
		scale=pow(2.0, int(shift)-int(w));
		T base=-scale*(probs.size()/2);
		
		for(unsigned i=0;i<probs.size();i++){
			table[i]=std::make_pair(base + scale*i, probs[i]);
		}
		
		return boost::make_shared<TableDistribution<T> >(table);
	}

}; // random
}; // flopoco

#endif
