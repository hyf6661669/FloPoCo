#ifndef random_base_distribution_hp
#define random_base_distribution_hp

#include <utility>
#include <stdexcept>
#include <vector>

#include <boost/smart_ptr.hpp>

#include "moment_conversions.hpp"

namespace flopoco
{
namespace random
{

template<class T> class Distribution;
template<class T> class ContinuousDistribution;
template<class T> class DiscreteDistribution;
template<class T> class EnumerableDistribution;

template<class T>
class Distribution
	: public boost::enable_shared_from_this<Distribution<T> >
{
public:
	virtual ~Distribution()
	{}
		
	//! k==1 -> mean, k==2 -> stddev, k==3 -> skewness, k==4 -> kurtosis, etc.
	virtual T StandardMoment(unsigned k) const=0;
		
	//! k==1 -> 0, k==2 -> variance, k==3 -> skewness * stddev^3, k==4 -> kurtosis * stddev^4
	virtual T CentralMoment(unsigned k) const
	{
		if(k==0) return 1;
		if(k==1) return 0;
		if(k==2) return pow(StandardMoment(2),2);
		if((k%2) && IsSymmetric()) return 0;
		return StandardMoment(k) * pow(StandardMoment(2),k);
	}
		
	virtual T RawMoment(unsigned k) const
	{
		// The default implementation isn't pretty, but it shouldn't be on the critical path
		std::vector<T> central(k+1);
		for(unsigned i=0;i<k+1;i++)
			central[i]=CentralMoment(i);
		return CentralMomentsToRawMoment(k, StandardMoment(1), &central[0]);
	}
		
	//! Return true if the distribution is definitely symmetric about the mean, false if it is not known.
	virtual bool IsSymmetric() const=0;
		
	virtual std::pair<T,T> Support() const=0;
		
	typedef boost::shared_ptr<Distribution> TypePtr;
};

template<class T>
class ContinuousDistribution
	: public Distribution<T>
{
public:
	virtual T Pdf(const T &x) const=0;
	virtual T Cdf(const T &x) const=0;
	virtual T InvCdf(const T &x) const=0;	

	typedef boost::shared_ptr<ContinuousDistribution> TypePtr;
};


template<class T>
class DiscreteDistribution
	: public Distribution<T>
{
public:
	//! If the set is finite then return a positive number, else return 0
	virtual uint64_t ElementCount() const=0;

	virtual int64_t IndexFromRange(const T &x) const=0;

	virtual T RangeFromIndex(int64_t x) const=0;

	virtual int64_t ClosestIndexFromRange(const T &x) const=0;

	//! If the range is regular (e.g. fixed point), then this is the step. Otherwise 0
	virtual T RangeGranularity() const=0;

	virtual T Pmf(const T &x) const=0;

	virtual T Cdf(const T &x) const=0;
	
	virtual std::pair<int64_t,int64_t> IndexSupport() const
	{
		if(ElementCount()==0)
			throw std::logic_error("IndexSupport - Cannot get support of infinite range.");
		std::pair<T,T> support=this->Support();
		return std::make_pair(IndexFromRange(support.first), IndexFromRange(support.second));
	}
	
	virtual T PmfByIndex(int64_t index) const
	{ return Pmf(RangeFromIndex(index)); }
	
	virtual T CdfByIndex(int64_t index) const
	{ return Cdf(RangeFromIndex(index)); }
	
	virtual std::pair<T,T> ElementByIndex(int64_t index) const
	{ return std::make_pair(RangeFromIndex(index), PmfByIndex(index)); }
	
	virtual void PmfByIndex(int64_t begin, int64_t end, T *pmf) const
	{
		while(begin!=end){
			*pmf++ = PmfByIndex(begin);
			++begin;
		}
	}
	
	virtual void CdfByIndex(int64_t begin, int64_t end, T *cdf) const
	{
		while(begin!=end){
			*cdf++ = CdfByIndex(begin);
			++begin;
		}
	}
	
	virtual void RangeByIndex(int64_t begin, int64_t end, T *range) const
	{
		while(begin!=end){
			*range++ = RangeFromIndex(begin);
			++begin;
		}
	}
	
	virtual void ElementsByIndex(int64_t begin, int64_t end, std::pair<T,T> *range) const
	{
		while(begin!=end){
			*range++ = ElementByIndex(begin);
			++begin;
		}
	}

	typedef boost::shared_ptr<DiscreteDistribution> TypePtr;
};

}; // random
}; // flopoco

#endif
