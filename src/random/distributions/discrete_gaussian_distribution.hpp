#ifndef flopoco_random_table_approx_discrete_gaussian_distribution_hpp
#define flopoco_random_table_approx_discrete_gaussian_distribution_hpp

#include <boost/math/distributions/normal.hpp>

#include "distribution.hpp"

namespace flopoco
{
namespace random
{
	
/* The asymptotic approximation for scale is excellent for sigma*2^fb > 2, and pretty much exact for sigma*2^fb > 4

Maxima:	
pmf_scale_exact(sig) := 1.0/(exp(0)+2*sum(exp(-(bfloat(x)^2/(2*sig^2))),x,-round(32*sig),-1));
pmf_scale_asymptotic(sig):=1/(sig*bfloat(sqrt(2*%pi)));
fpprec:64;
for sig:0.5b0 thru 4b0 step 0.5 do print(sig," : ", (pmf_scale_exact(sig)-pmf_scale_asymptotic(sig))/pmf_scale_exact(sig));
Output:
5.0b-1" : "-1.438377206222871346702909896491768299278882392381925898141489348b-2
1.0b0" : "-5.350575982148479362482248080537060646957443172632755077611545768b-9
1.5b0" : "-1.029470036899271965450035408031801859591228429726213862053557989b-19
2.0b0" : "-1.024500455847086035104711735299867411042392431167430232970186196b-34
2.5b0" : "-5.273462610385365704279771360669160093311451328047577604104382393b-54
3.0b0" : "-3.570278789279349475689537591100337425205829189347718849660396188b-65
3.5b0" : "0.0b0
4.0b0" : "0.0b0
	
In terms of moments, the errors in the first few moments are good at similar limits:
pmf_mom_exact(sig,k) := pmf_scale_exact(sig)*2*sum(exp(-(bfloat(x)^2/(2*sig^2)))*x^k,x,-round(32*sig),-1);
pmf_mom_asymptotic(sig,k) :=
    if k=0 then 1
    elseif k=1 then 0
    elseif k=2 then sig*sig
    elseif mod(k,2)=1 then 0
    else sig*sig*(k-1)*pmf_mom_asymptotic(sig,k-2);
log2(x):=log(x)/log(bfloat(2));
errbits(a,e) := if a=e then -inf else ceiling(log2(abs((a-e)/e)));
pmf_errbits(sig,k) := errbits(pmf_mom_asymptotic(sig,k),pmf_mom_exact(sig,k));
for sig:0.5b0 thru 4.0b0 step 0.5 do(
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
"sig="5.0b-1" m2=2^"-2" m4=2^"-2" m8=2^"-2" m16=2^"-2" m32=2^"-4" m64=2^"-2" m128=2^"-2
"sig="1.0b0" m2=2^"-22" m4=2^"-18" m8=2^"-14" m16=2^"-13" m32=2^"-13" m64=2^"-14" m128=2^"-14
"sig="1.5b0" m2=2^"-56" m4=2^"-51" m8=2^"-44" m16=2^"-34" m32=2^"-30" m64=2^"-32" m128=2^"-30
"sig="2.0b0" m2=2^"-105" m4=2^"-99" m8=2^"-90" m16=2^"-76" m32=2^"-59" m64=2^"-58" m128=2^"-55
"sig="2.5b0" m2=2^"-169" m4=2^"-162" m8=2^"-152" m16=2^"-135" m32=2^"-110" m64=2^"-87" m128=2^"-88
"sig="3.0b0" m2=2^"-inf" m4=2^"-214" m8=2^"-213" m16=2^"-208" m32=2^"-179" m64=2^"-140" m128=2^"-126
"sig="3.5b0" m2=2^"-214" m4=2^"-inf" m8=2^"-214" m16=2^"-214" m32=2^"-214" m64=2^"-213" m128=2^"-172
"sig="4.0b0" m2=2^"-214" m4=2^"-214" m8=2^"-213" m16=2^"-214" m32=2^"-212" m64=2^"-212" m128=2^"-212
(this calculations use exact rather than asymptotic scaling, but the results are the same either way).

Note that 2^-212 is about the precision when fpprec:64 (i.e. for 64 decimal digits of precision), so for
sigma*2^fb >= 4 all moments up to 128 are essentially exact, or alternatively, with sigma=1 anything more
with more than 2 fractional bits has perfect moments.

*/

template<class T>
class DiscreteGaussianDistribution
	: public DiscreteDistribution<T>
{
private:	
	T m_stddev, m_one, m_zero;
	int m_fb;

	T m_base, m_step;
	std::vector<T> m_cdf;
public:
	DiscreteGaussianDistribution(T stddev, int fb)
		: m_stddev(stddev)
		, m_one(stddev/stddev)
		, m_zero(m_one-m_one)
		, m_fb(fb)
	{
		if(ldexp(stddev, fb) < 2)
			throw std::string("DiscreteGaussianDistribution - Must have stddev*2^fb >= 2 for the approximation to be valid.");
		
		T one=m_one;
		T zero=m_zero;
		
		// At the moment we arbitrarily go out to +-36sigma. This means that in the +-32sigma range the cdf should be pretty accurate
		int points=1+boost::math::tools::real_cast<long>(round(36*m_stddev*pow(2.0,fb)));

		m_step=one*pow(2.0,-fb);
		m_base=m_step*(points-1);
		m_cdf.resize(points);
		
		T x_scale=1 / (stddev*stddev*2);
		
		if(ldexp(stddev,fb) <4){
			// Calculate exact scale
			
			T acc=zero;
			T curr=m_base;
			for(int i=points-1;i>=0;i--){
				T p=exp(-curr*curr * x_scale);
				acc += p;
				m_cdf.at(i)=acc;
				curr -= m_step;
			}
			assert(curr==-m_step);	// stepped one delta past zero
			
			T scale=1 / (2*acc - 1);	// 1==exp(-0^2/2)
			
			for(int i=points-1;i>=0;i--){
				m_cdf[i] *= scale;
			}
		}else{
			// Use asymptotic scale, as it works very well in this range
			
			T pi=2*acos(zero);
			T pdf_scale=1 / (pow(2.0,fb)*stddev* sqrt(2*pi));
			
			T acc=zero;
			T curr=m_base;
			for(int i=points-1;i>0;i--){
				T p=pdf_scale*exp(-curr*curr * x_scale);
				acc += p;
				m_cdf[i]=acc;
				curr -= m_step;
			}
			assert(2*acc<1);
			assert(curr==0);
			// Doing it this way means that things add up closer to one, compared with m_cdf[0]=acc+pdf_scale*exp(0)
			// which may have accumulated errors. Probably not a big deal.
			m_cdf[0]=1-acc;//=  acc+2*(0.5-acc);
		}
		
		//for(int i=0;i<16;i++){
		//	std::cerr<<i<<", "<<i*m_step<<", "<<m_cdf[i]<<", "<<m_cdf[i]-m_cdf[i+1]<<"\n";
		//}
	}
	
	//! k==1 -> mean, k==2 -> stddev, k==3 -> skewness, k==4 -> kurtosis, etc.
	virtual T StandardMoment(unsigned k) const
	{
		// See comments at start about accuracy of moments. Essentially we have
		// limited ourselves to situations where the discrete gaussian moments are
		// very close to the true gaussian moments (to better than double-precision),
		// so there is no need to calculate the actual moments.  This may need to be
		// revisited, but seems reasonable for now. If necessary, it is cheap to calculate
		// exactly, as the asymptotics only differ when there are very few points covered
		// by the distribution (i.e. stddev*2^fb is small).
		if(k==0)	return 1.0;
		if(k==1)	return 0.0;
		if(k==2) return m_stddev;
		if((k%2)==1) return 0.0;
		T base=1;
		for(unsigned i=4;i<=k;i+=2){
			base=base*(i-1);
		}
		return base;
	}
		
	virtual bool IsSymmetric() const
	{ return true; }
	
	virtual uint64_t ElementCount() const
	{ return 0; } // In principle it is infinite
		
	virtual std::pair<T,T> Support() const
	{ return std::pair<T,T>(-m_base+m_step, m_base-m_step); }
	
	virtual T RangeGranularity() const
	{ return m_step; }
	
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
		T q=ldexp(x, m_fb);
		return boost::math::tools::real_cast<long>(round(q));
	}

	virtual T RangeFromIndex(int64_t x) const
	{
		return ldexp(m_one*(T)long(x), -m_fb);
	}
	
	virtual T Pmf(const T &x) const
	{
		T tmp1=ldexp(x, m_fb);
		T tmp2=round(tmp1);
		if(tmp1!=tmp2)
			return 0;
		long i=boost::math::tools::real_cast<long>(tmp2);
		i=abs(i);
		if((i+1)>=(int)m_cdf.size()){
			throw std::string("DiscreteGaussianDistribution::Pmf - x out of range.");
		}
		return m_cdf[i]-m_cdf[i+1];
	}
	
	virtual T Cdf(const T &x) const
	{
		int i=boost::math::tools::real_cast<long>(floor(ldexp(x, m_fb)));
		i=abs(i);
		if(i>=(int)m_cdf.size()){
			throw std::string("DiscreteGaussianDistribution::Pmf - x out of range.");
		}
		return (x<=0) ? m_cdf[i] : 1-m_cdf[i];
	}
	
	virtual T PmfByIndex(int64_t i) const
	{
		i=std::abs(i);
		if((i+1)>=(int)m_cdf.size()){
			throw std::string("DiscreteGaussianDistribution::Pmf - x out of range.");
		}
		return m_cdf[i]-m_cdf[i+1];
	}
	
	virtual T CdfByIndex(int64_t i) const
	{
		int64_t ai=abs(i);
		if(ai>=(int)m_cdf.size()){
			throw std::string("DiscreteGaussianDistribution::Pmf - x out of range.");
		}
		return (i<=0) ? m_cdf[ai] : 1-m_cdf[ai];
	}
		
	typedef boost::shared_ptr<DiscreteGaussianDistribution> DiscreteGaussianDistributionPtr;
};

}; // random
}; // flopoco

#endif
