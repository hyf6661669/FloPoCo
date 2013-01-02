#ifndef flopoco_random_table_approx_quantised_gaussian_distribution_hpp
#define flopoco_random_table_approx_quantised_gaussian_distribution_hpp

#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/erf.hpp>

#include "distribution.hpp"

namespace flopoco
{
namespace random
{
	
/*
	
This is the "obvious" discrete distribution, where quantiles are correct, but the moments
are completely wrong.
	
Maxima:
load(distrib);
kill(Phi);
Phi[x] := bfloat(cdf_normal(bfloat(x),0,1));
mom_k_true(sig,k) :=
    if k=0 then 1
    elseif k=1 then 0
    elseif k=2 then sig*sig
    elseif mod(k,2)=1 then 0
    else sig*sig*(k-1)*mom_k_true(sig,k-2);
mom_k_q(sig,k):=block([x_scale:1/sig,acc:0b0],
    for i:-round(32*sig) thru -1 do(
        acc:acc+(Phi[(i+0.5)*x_scale]-Phi[(i-0.5)*x_scale])*bfloat(i)^k
    ),
    2*acc
);
fpprec:32;

log2(x):=log(x)/log(bfloat(2));
errbits(a,e) := if a=e then -inf else ceiling(log2(abs((a-e)/e)));
pmf_errbits(sig,k) := errbits(mom_k_q(sig,k),mom_k_true(sig,k));
for log2_sig:1.0b0 thru 8.0b0 step 1.0 do(
    sig:2^log2_sig,
    print(
        "sig=",sig,
        " m2=2^",pmf_errbits(sig,2),
        " m4=2^",pmf_errbits(sig,4),
        " m8=2^",pmf_errbits(sig,8),
        " m16=2^",pmf_errbits(sig,16),
        " m32=2^",pmf_errbits(sig,32),
        " m64=2^",pmf_errbits(sig,64),
        " m128=2^",pmf_errbits(sig,128)
    )
  );
 
 Output:
 
 "sig="2.0b0" m2=2^"-5" m4=2^"-4" m8=2^"-3" m16=2^"-2" m32=2^"-1" m64=2^"0" m128=2^"1
"sig="4.0b0" m2=2^"-7" m4=2^"-6" m8=2^"-5" m16=2^"-4" m32=2^"-3" m64=2^"-2" m128=2^"-1
"sig="8.0b0" m2=2^"-9" m4=2^"-8" m8=2^"-7" m16=2^"-6" m32=2^"-5" m64=2^"-4" m128=2^"-7
"sig="1.6b1" m2=2^"-11" m4=2^"-10" m8=2^"-9" m16=2^"-8" m32=2^"-7" m64=2^"-6" m128=2^"-3
"sig="3.2b1" m2=2^"-13" m4=2^"-12" m8=2^"-11" m16=2^"-10" m32=2^"-9" m64=2^"-8" m128=2^"-3
"sig="6.4b1" m2=2^"-15" m4=2^"-14" m8=2^"-13" m16=2^"-12" m32=2^"-11" m64=2^"-10" m128=2^"-3
"sig="1.28b2" m2=2^"-17" m4=2^"-16" m8=2^"-15" m16=2^"-14" m32=2^"-13" m64=2^"-12" m128=2^"-3
"sig="2.56b2" m2=2^"-19" m4=2^"-18" m8=2^"-17" m16=2^"-16" m32=2^"-15" m64=2^"-14" m128=2^"-3

So the standard-deviation is generally awful, and as a consequence all the higher order moments
rapidly fall apart.

It would be tempted to try to correct such that m2 is correct, in which case the higher order moments
would also become more accurate, but it is very difficult to achieve. For example, we could choose a
correction of the form:
pmf(sig,alpha,x) = Phi( (x+0.5)*alpha/sig ) - Phi( (x-0.5)*alpha/sig )
We could then state it as an optimisation of the form:
  sig^2 = sum(x^2*pmf(sig,alpha,x),x=-inf..+inf)
and search for alpha. This is equivalent to searching for a different
standard deviation sig', such that:
  sig^2 = mom_k_q(sig',2)
which could be done using bisection search. However, it's not very satisfying, and is
going to be pretty damn slow as we have to do ~32*sig evaluations of Phi for
every search step.

Maxima:
std_quick(sig):=block([
        acc:0.0,
        ss:float(sig)
    ],
    for i:-float(round(32*sig)) thru -1 do(
        acc:acc+(float(cdf_normal(i+0.5,0,ss))-float(cdf_normal(i-0.5,0,ss)))*i^2
    ),
    2*acc
);

find_sig_alt(sig):=block([hi_sig:sig*1.01,lo_sig:sig*0.99,hi_std,lo_std,mid_sig, mid_std],
    do(
        lo_std:sqrt(std_quick(lo_sig)),
        if lo_std < sig then return (),
        lo_sig:lo_sig*0.99
    ),
    do(
        hi_std:sqrt(std_quick(hi_sig)),
        if hi_std > sig then return (),
        hi_sig:hi_sig*1.01
    ),
    do(
        mid_sig:(lo_sig+hi_sig)/2,
        mid_std:sqrt(std_quick(mid_sig)),
        if abs((mid_std-sig)/sig)<1e-12 then
            return([mid_sig,mid_std]),
        if mid_std > sig then (
            hi_sig:mid_sig,
            hi_std:mid_std
        )else(
            lo_sig:mid_sig,
            lo_std:mid_std
        )
    )
);

for log2_sig:1.0b0 thru 6.0b0 step 1.0 do(
    sig:2^log2_sig,
    sigalt:find_sig_alt(sig)[1],
    print(
        "sig=",sig,
        "sigalt=",sigalt,
        " m2=2^",errbits(mom_k_q(sigalt,2),mom_k_true(sig,2)),
        " m4=2^",errbits(mom_k_q(sigalt,4),mom_k_true(sig,4)),
        " m8=2^",errbits(mom_k_q(sigalt,8),mom_k_true(sig,8)),
        " m16=2^",errbits(mom_k_q(sigalt,16),mom_k_true(sig,16)),
        " m32=2^",errbits(mom_k_q(sigalt,32),mom_k_true(sig,32)),
        " m64=2^",errbits(mom_k_q(sigalt,64),mom_k_true(sig,64)),
        " m128=2^",errbits(mom_k_q(sigalt,128),mom_k_true(sig,128))
    )
  );
  
  Output:
  
 "sig="2.0b0"sigalt="1.979057014505937677b0" m2=2^"-41" m4=2^"-12" m8=2^"-9" m16=2^"-7" m32=2^"-5" m64=2^"-3
" m128=2^"0
"sig="4.0b0"sigalt="3.9895697345305234101b0" m2=2^"-39" m4=2^"-16" m8=2^"-13" m16=2^"-11" m32=2^"-9" m64=2^"-
6" m128=2^"0
"sig="8.0b0"sigalt="7.9947899701446294738b0" m2=2^"-39" m4=2^"-20" m8=2^"-17" m16=2^"-15" m32=2^"-13" m64=2^"-
7" m128=2^"0
"sig="1.6b1"sigalt="1.5997395621370524166b1" m2=2^"-43" m4=2^"-24" m8=2^"-21" m16=2^"-19" m32=2^"-17" m64=2^"-
6" m128=2^"0
"sig="3.2b1"sigalt="3.1998697890192270278b1" m2=2^"-39" m4=2^"-28" m8=2^"-25" m16=2^"-23" m32=2^"-21" m64=2^"-
7" m128=2^"0
"sig="6.4b1"sigalt="6.3999348955005407333b1" m2=2^"-40" m4=2^"-32" m8=2^"-29" m16=2^"-27" m32=2^"-25" m64=2^"-
7" m128=2^"0

So that technique sucks as well, and is much worse than the uncorrected one for uncorrected moments.


*/

template<class T>
class QuantisedGaussianDistribution
	: public DiscreteDistribution<T>
{
private:	
	T m_one, m_zero, m_stddev;
	int m_fb;

	boost::math::normal_distribution<T> m_dist;

	T m_step;
	mutable std::vector<T> m_cdf;

	mutable std::vector<T> m_moments;

	double ERFC(double x) const
	{
		return boost::math::erfc(x);
	}

	mpfr::mpreal ERFC(mpfr::mpreal x) const
	{
		mpfr::mpreal res;
		mpfr_erfc(res.mpfr_ptr(), x.mpfr_ptr(), MPFR_RNDN);
		return res;
	}

	void ExtendCdf(int target) const
	{
		int curr=((int)m_cdf.size())-1;
		if(curr>target)
			return;
		
		T x_step=m_step/(m_stddev*sqrt(m_one*2));
		
		T x=x_step*curr - (x_step/2);
		m_cdf.resize(target+1);
		while(curr<target){
			curr++;
			x += x_step;
			T tmp=ERFC(x)*0.5;
			m_cdf.at(curr)=tmp;
			assert(tmp>0);
		}
	}
	
	T PmfAt(int64_t index) const
	{
		index=std::abs(index);
		if(index+1>=(int64_t)m_cdf.size())
			ExtendCdf(index+1);
		return m_cdf[index]-m_cdf[index+1];
	}
	
	T CdfAt(int64_t index) const
	{
		if(index>0)
			return 1-CdfAt(-index-1);
		index=-index;
		if(index+1>=(int64_t)m_cdf.size())
			ExtendCdf(index+1);
		return m_cdf.at(index);
		
		// This works, but only for doubles. erfc seems broken for mpreal
		//return boost::math::cdf(m_dist, T(index)*m_step+m_step/2);
		
		// This works for doubles and mpreals
		//T x=(T(index)*m_step+m_step/2) / (m_stddev * sqrt(2*m_one));
		//return 0.5*ERFC(-x);
	}
public:
	QuantisedGaussianDistribution(T stddev, int fb)
		: m_one(stddev/stddev)
		, m_zero(m_one-m_one)
		, m_stddev(stddev)
		, m_fb(fb)
		, m_dist(0, m_stddev)
	{
		m_step=m_one*pow(2.0,-fb);
	}
	
	//! k==1 -> mean, k==2 -> stddev, k==3 -> skewness, k==4 -> kurtosis, etc.
	virtual T StandardMoment(unsigned k) const
	{
		if(k==0) return 1.0;
		if((k%2)==1) return 0.0;
		
		if(k<m_moments.size()){
			if(m_moments[k]!=0)
				return m_moments[k];
		}
		
		// We want to find x such that x^k*Pmf[x] + Pmf[0] = Pmf[0]
		// Note that in high-precisions this may be a _very_ long way away, so be careful.
		T origin_p=PmfAt(0);
		int start=8;
		while(1){
			start=(start*4)/3;
			T start_p=pow(start*m_step, k) * PmfAt(start);
			T tmp=start_p+origin_p;
			if(tmp==origin_p)
				break;
		}
		
		T acc=m_zero;
		T curr=start*m_step;
		for(int i=start;i>0;i--){
			acc += pow(curr, k) * (m_cdf[i]-m_cdf[i+1]);
			curr -= m_step;
		}
		acc *=2;
		
		T res;
		
		if(k==2){
			res=sqrt(acc);
		}else{
			res=acc / pow(StandardMoment(2), k);
		}
		
		if(m_moments.size()<=k)
			m_moments.resize(k+1, m_zero);
		m_moments[k]=res;
		
		return res;
	}
	
	virtual uint64_t ElementCount() const
	{ return 0; }
	
	virtual int64_t IndexFromRange(const T &x) const
	{
		T q=ldexp(x, m_fb);
		T r=round(q);
		if(q!=r)
			throw std::string("IndexFromRange - Not aligned to grid.");
		return boost::math::tools::real_cast<long>(r);
	}
	
	virtual int64_t ClosestIndexFromRange(const T &x) const
	{
		T q=round(ldexp(x, m_fb));
		return boost::math::tools::real_cast<long>(q);
	}

	virtual T RangeFromIndex(int64_t x) const
	{
		return ldexp(m_one*T((long)x), -m_fb);
	}
		
	virtual bool IsSymmetric() const
	{ return true; }
		
	virtual std::pair<T,T> Support() const
	{
		T inf=std::numeric_limits<T>::infinity();
		return std::pair<T,T>(-inf,inf); 
	}
	
	virtual T Pmf(const T &x) const
	{
		T q=ldexp(x, m_fb);
		int i=boost::math::tools::real_cast<long>(round(q));
		if(q!=i)
			return m_zero;
		return PmfAt(i);
	}
	
	virtual T Cdf(const T &x) const
	{
		T q=ldexp(x, m_fb);
		int64_t i=boost::math::tools::real_cast<long>(floor(q));
		
		return CdfAt(i);
	}
	
	virtual T PmfByIndex(int64_t index) const
	{ return PmfAt(index); }
	
	virtual T CdfByIndex(int64_t index) const
	{ return CdfAt(index); }
	
		
	typedef boost::shared_ptr<QuantisedGaussianDistribution> QuantisedGaussianDistributionPtr;
};

}; // random
}; // flopoco

#endif
