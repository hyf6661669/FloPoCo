#ifndef random_base_histogram_distribution_hpp
#define random_base_histogram_distribution_hpp

#include "distribution.hpp"
#include "moment_conversions.hpp"

#include <vector>
#include <algorithm>
#include <numeric>

#include "random/utils/mpreal/boost_math_mpreal.hpp"
#include "random/utils/sum.hpp"
#include "boost/math/tools/real_cast.hpp"

namespace flopoco
{
namespace random
{

/* Histogram distribution is disguished from table by having a regularly
	spaced set of discrete points */
template<class T>
class HistogramDistribution
	: public DiscreteDistribution<T>
{
public:	
	typedef boost::shared_ptr<HistogramDistribution> TypePtr;
private:
	bool m_isSymmetric;

	typedef std::vector<T> storage_t;
	storage_t m_elements, m_cdf;
	T m_base, m_step;

	T m_zero, m_one;
	
	mutable std::vector<T> m_standardMoments;
public:

	template<class TC>
	HistogramDistribution(T base, T step, const TC &src)
		: m_isSymmetric(true)
		, m_elements(src)
		, m_base(base)
		, m_step(step)
		, m_zero(src[0]-src[0])
		, m_one(m_zero+1)
	{
		if(src.size()==0)
			throw std::invalid_argument("HistogramDistribution - Table must contain at least one element.");
		
		m_cdf.resize(src.size());
		
		T acc=m_zero;
		for(int i=0;i<(int)src.size();i++){
			if(m_elements[i] < 0){
				throw std::logic_error("HistogramDistribution - Negative element probability.");
			}
			acc+=m_elements[i];
			m_cdf[i]=acc;
		}
		
		if(fabs(acc-1)>1e-12)
			throw std::logic_error("HistogramDistribution - Element probabilities are more than 1e-12 from one.");
		
		T scale=m_one/acc;
		for(unsigned i=0;i<src.size();i++){	
			m_elements[i] = src[i] * scale;
			m_cdf[i] = m_cdf[i] * scale;
		}
		
		for(unsigned i=0;i<src.size()/2;i++){
			if(m_elements[i] != m_elements[m_elements.size()-i-1]){
				m_isSymmetric=false;
				break;
			}
		}
	}
	
	T RangeBase() const
	{ return m_base; }
	
	T RangeStep() const
	{ return m_step; }
	
	virtual T RangeGranularity() const
	{ return m_step;}
	
	virtual T StandardMoment(unsigned k) const
	{
		while(m_standardMoments.size()<=k){
			typename SelectAccumulator<T>::type acc;
			unsigned kc=m_standardMoments.size();
			if((kc%2) && m_isSymmetric){
				m_standardMoments.push_back(m_zero);
			}else{
				acc=m_zero;
				T curr=m_base;
				if(k>1)
					curr=curr-m_standardMoments[1];
				if(m_isSymmetric){
					for(size_t i=0;i<m_elements.size()/2;i++){
						acc += 2 * pow(curr, kc) * m_elements[i];
						curr+=m_step;
					}
					if(m_elements.size()%2){
						acc += pow(curr, kc) * m_elements[m_elements.size()/2];
					}
				}else{
					for(size_t i=0;i<m_elements.size();i++){
						acc += pow(curr, kc) * m_elements[i];
						curr+=m_step;
					}
				}
				if(k==2)
					acc=sqrt(acc);
				if(k>2)
					acc=acc / pow(m_standardMoments[2],k);
				m_standardMoments.push_back(acc);
			}
		}
		return m_standardMoments[k];
	}
	
	virtual bool IsSymmetric() const
	{ return m_isSymmetric; } // this is more strict than symmetric about the mean, but still is valid
		
	virtual std::pair<T,T> Support() const
	{ return std::make_pair(m_base, m_base+m_step*(m_elements.size()-1)); }
		
	virtual T Pmf(const T &x) const
	{
		if(x<m_base)
			return m_zero;
		if(x > (m_base+m_step*(m_elements.size()-1)))
			return m_zero;
		T vi=(x-m_base)/m_step;
		T ii=round(vi);
		if(vi!=ii)
			return m_zero;
		size_t i=boost::math::tools::real_cast<long>(ii);
		assert(i<m_elements.size());
		return m_elements[i];
	}
	
	virtual T Cdf(const T &x) const
	{
		T vi=(x-m_base)/m_step;
		T ii=floor(vi);
		if(ii<0)
			return m_zero;
		if(ii>=m_cdf.size())
			return m_one;
		size_t i=boost::math::tools::real_cast<long>(ii);
		return m_cdf[i];
	}

	virtual uint64_t ElementCount() const
	{ return m_elements.size(); }
	
	virtual int64_t IndexFromRange(const T &x) const
	{
		T vi=(x-m_base)/m_step;
		T ii=round(vi);
		if(vi!=ii)
			throw std::invalid_argument("RangeToIndex - Value is not aligned to range.");
		if((ii<0) || (ii>=m_elements.size()))
			throw std::invalid_argument("RangeToIndex - Value is not aligned to range.");
		return boost::math::tools::real_cast<long>(ii);
	}
	
	virtual int64_t ClosestIndexFromRange(const T &x) const
	{
		T vi=(x-m_base)/m_step;
		T ii=round(vi);
		ii=std::max(T(0), std::min(T(m_elements.size()), ii));
		return boost::math::tools::real_cast<long>(ii);
	}

	virtual T RangeFromIndex(int64_t index) const
	{
		if((index<0) || (index>=m_elements.size()))
			throw std::invalid_argument("IndexToRange - Index is not valid.");
		T xx=index;
		return m_base+xx*m_step;
	}
	
	virtual T PmfByIndex(int64_t index) const
	{
		assert((index>=0) && (index<m_elements.size()));
		return m_elements[index];
	}
	
	virtual T CdfByIndex(int64_t index) const
	{ 
		assert((index>=0) && (index<m_cdf.size()));
		return m_cdf[index];
	}
	
	virtual void PmfByIndex(int64_t begin, int64_t end, T *pmf) const
	{
		if((begin<0) || (end>m_elements.size()))
			throw std::range_error("PmfByIndex - indices out of range.");
		if(begin>end)
			throw std::range_error("PmfByIndex - indices are not ordered.");
		
		std::copy(m_elements.begin()+begin, m_elements.begin()+end, pmf);
	}
	
	virtual void CdfByIndex(int64_t begin, int64_t end, T *cdf) const
	{
		if((begin<0) || (end>m_elements.size()))
			throw std::range_error("PdfByIndex - indices out of range.");
		if(begin>end)
			throw std::range_error("PdfByIndex - indices are not ordered.");
		
		std::copy(m_cdf.begin()+begin, m_cdf.begin()+end, cdf);
	}
};

}; // random
}; // flopoco

#endif
