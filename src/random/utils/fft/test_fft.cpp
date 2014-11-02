#include "fft_mpreal.hpp"
#include "fft_double.hpp"

#include "convolve_mpreal.hpp"

#include <vector>
#include <stdlib.h>

void compare(unsigned n, unsigned wPrec)
{
	std::vector<double> raw(2*n);
	
	for(unsigned i=0;i<2*n;i++){
		raw[i]=drand48();
	}
	
	std::vector<mpfr::mpreal> alt(2*n);
	for(unsigned i=0;i<2*n;i++){
		alt[i]=mpfr::mpreal(raw[i], wPrec);
	}
	
	flopoco::random::fft_radix2_complex(&raw[0], n);
	flopoco::random::fft_radix2_complex(&alt[0], n, wPrec);
	
	//flopoco::random::ifft_radix2_complex(&raw[0], n);
//	flopoco::random::ifft_radix2_complex(&alt[0], n, wPrec);
	
	for(unsigned i=0;i<n;i++){
		std::cerr<<i<<" : dbl("<<raw[2*i]<<","<<raw[2*i+1]<<"), mpfr("<<alt[2*i]<<","<<alt[2*i+1]<<")\n";
		std::cerr<<i<<" : "<< (alt[2*i]-raw[2*i])/raw[2*i]<<"\n";
	}
}

void compare_two(unsigned n, unsigned wPrec)
{
	gmp_randstate_t state;
	
	gmp_randinit_default (state);
	
	std::vector<mpfr::mpreal> low(2*n, mpfr::mpreal(0, wPrec)), high(2*n, mpfr::mpreal(0, 2*wPrec));
	
	for(unsigned i=0;i<n;i++){
		//mpfr_grandom (get_mpfr_ptr(low[2*i]), get_mpfr_ptr(high[2*i+1]), state, MPFR_RNDN);
		mpfr_urandomb(get_mpfr_ptr(low[2*i]), state);
		mpfr_urandomb(get_mpfr_ptr(low[2*i+1]), state);
		mpfr_set(get_mpfr_ptr(high[2*i]), get_mpfr_ptr(low[2*i]), MPFR_RNDN);
		mpfr_set(get_mpfr_ptr(high[2*i+1]), get_mpfr_ptr(low[2*i+1]), MPFR_RNDN);
	}
	
	flopoco::random::fft_radix2_complex(&low[0], n, wPrec);
	flopoco::random::fft_radix2_complex(&high[0], n, 2*wPrec);
	
	mpfr::mpreal errMax(0), errSumSqr(0);
	
	for(unsigned i=0;i<2*n;i++){
		mpfr::mpreal err=mpfr::abs( (low[i]-high[i]) / high[i] );
		errMax=std::max(errMax, err);
		errSumSqr += err*err;
	}
	
	mpfr::mpreal errRms=sqrt(errSumSqr/(2*n));
	
	std::cout<<n<<", "<<wPrec<<", 2^"<<log2(errMax)<<", 2^"<<log2(errRms)<<"\n";
	
	gmp_randclear(state);
}

void compare_three(unsigned n, unsigned wPrec)
{
	gmp_randstate_t state;
	
	gmp_randinit_default (state);
	
	std::vector<mpfr::mpreal> low(2*n, mpfr::mpreal(0, wPrec));
	std::vector<double> high(2*n);
	
	for(unsigned i=0;i<n;i++){
		//mpfr_grandom (get_mpfr_ptr(low[2*i]), get_mpfr_ptr(high[2*i+1]), state, MPFR_RNDN);
		mpfr_urandomb(get_mpfr_ptr(low[2*i]), state);
		mpfr_urandomb(get_mpfr_ptr(low[2*i+1]), state);
		mpfr_log(get_mpfr_ptr(low[2*i]), get_mpfr_ptr(low[2*i]), MPFR_RNDN);
		mpfr_log(get_mpfr_ptr(low[2*i+1]), get_mpfr_ptr(low[2*i+1]), MPFR_RNDN);
		high[2*i]=low[2*i].toDouble();
		high[2*i+1]=low[2*i+1].toDouble();
		//high[2*i]=drand48();
		//high[2*i+1]=drand48();
		//low[2*i]=high[2*i];
		//low[2*i+1]=high[2*i+1];
	}
	
	std::vector<mpfr::mpreal> orig(low);
	
	flopoco::random::fft_radix2_complex(&low[0], n, wPrec+1+(unsigned)log2(n));
	flopoco::random::ifft_radix2_complex(&low[0], n, wPrec+1+(unsigned)log2(n));
	flopoco::random::fft_radix2_complex(&high[0], n);
	flopoco::random::ifft_radix2_complex(&high[0], n);
	
	double errMax(0), errSumSqr(0);
	double errAltMax(0), errAltSumSqr(0);
	
	for(unsigned i=0;i<2*n;i++){
		//std::cerr<<low[i]<<", "<<high[i]<<"\n";
		mpfr::mpreal correct(orig[i]);
		correct.setPrecision(128);
		double err=std::fabs( ((low[i]-correct) / correct).toDouble() );
		double errAlt=std::fabs( ((high[i]-correct) / correct).toDouble());
		if(err>=1)
			std::cerr<<"correct="<<correct<<", got="<<low[i]<<"\n";
		errMax=std::max(errMax, err);
		errSumSqr += err*err;
		errAltMax=std::max(errAltMax, errAlt);
		errAltSumSqr += errAlt*errAlt;
	}
	
	mpfr::mpreal errRms=sqrt(errSumSqr/(2*n));
	mpfr::mpreal errAltRms=sqrt(errAltSumSqr/(2*n));
	
	std::cout<<n<<", "<<wPrec<<", 2^"<<log2(errMax)<<", 2^"<<log2(errRms)<<" | 2^"<<log2(errAltMax)<<", 2^"<<log2(errAltRms)<<"\n";
	
	gmp_randclear(state);
}

void compare_convolve(int n, int k, int prec)
{
	std::vector<mpfr::mpreal> ones(n, mpfr::mpreal(1,prec));
	
	std::vector<mpfr::mpreal>  res(flopoco::random::self_convolve(ones, k));
	
	for(int i=0;i<res.size()/2;i++){
		mpfr::mpreal fwd=res[i], rev=res[res.size()-i-1];
		
		std::cerr<<i<<", "<<fwd<<", "<<rev<<", "<<fwd-rev<<"\n";
	}
}

int main(int argc, char *argg[])
{
	srand48(1);
	
	unsigned n=16;
	unsigned wPrec=1024;
	
	std::vector<double> raw(2*n);
	std::vector<mpfr::mpreal> alt(2*n, mpfr::mpreal(0,wPrec));
	
	raw[0]=1;
	alt[0]=1;
	raw[2*n-2]=1;
	alt[2*n-2]=1;
	
	flopoco::random::fft_radix2_complex(&alt[0], n, wPrec);
	flopoco::random::fft_radix2_complex(&raw[0], n);
	
	flopoco::random::ifft_radix2_complex(&alt[0], n, wPrec);
	flopoco::random::ifft_radix2_complex(&raw[0], n);
	
	for(unsigned i=0;i<n;i++){
		//if(i>n-128){
			std::cerr<<i<<" : dbl("<<raw[2*i]<<","<<raw[2*i+1]<<"), mpfr("<<alt[2*i]<<","<<alt[2*i+1]<<")\n";
		//}
	}
	
	for(unsigned i=8;i<=1024;i*=2){
		for(unsigned w=16;w<=128;w+=4){
			//compare_three(i,w);
		}
	}
	
	for(unsigned i=8;i<=1024;i*=4){
		for(unsigned w=16;w<=64;w+=16){
			//compare_two(i,w);
		}
	}
	
	compare_convolve(256, 16, 128);
	
	return 0;
}
