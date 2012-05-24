#ifndef random_base_histogram_distribution_hpp
#define random_base_histogram_distribution_hpp

#include "distribution.hpp"
#include "moment_conversions.hpp"

#include <vector>
#include <algorithm>
#include <numeric>

namespace flopoco
{
namespace random
{

/* Histogram distribution is disguished from table by having a regularly
	spaced set of discrete points */
template<class T>
class HistogramDistribution
	: public EnumerableDistribution<T>
{
public:	
	typedef boost::shared_ptr<HistogramDistribution> TypePtr;
private:
	bool m_isSymmetric;

	typedef std::vector<T> storage_t;
	storage_t m_elements;
	T m_base, m_step;
	
	mutable std::vector<T> m_standardMoments;
public:

	template<class TC>
	HistogramDistribution(T base, T step, const TC &src)
		: m_isSymmetric(true)
		, m_elements(src)
		, m_base(base)
		, m_step(step)
	{
		if(src.size()==0)
			throw std::invalid_argument("HistogramDistribution - Table must contain at least one element.");
		
		T acc=sum(src.begin(), src.end());
		
		if(fabs(acc-1)>1e-12)
			throw std::logic_error("HistogramDistribution - Element probabilities are more than 1e-12 from one.");
		T scale=1.0/acc;
		for(unsigned i=0;i<src.size();i++){
			if(m_elements[i] < 0){
				throw std::logic_error("HistogramDistribution - Negative element probability.");
			}
			m_elements[i] = src[i] * scale;
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
	
	virtual T StandardMoment(unsigned k) const
	{
		while(m_standardMoments.size()<=k){
			typename SelectAccumulator<T>::type acc;
			unsigned kc=m_standardMoments.size();
			if((kc%2) && m_isSymmetric){
				m_standardMoments.push_back(0);
			}else{
				acc=0.0;
				T curr=m_base;
				for(size_t i=0;i<m_elements.size();i++){
					acc += pow(curr, kc) * m_elements[i];
					curr+=m_step;
				}
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
			return 0;
		if(x > (m_base+m_step*(m_elements.size()-1)))
			return 0;
		T vi=(x-m_base)/m_step;
		T ii=round(vi);
		if(vi!=ii)
			return 0;
		size_t i=(size_t)ii;
		assert(i<m_elements.size());
		return m_elements[i];
	}
	
	virtual T Cdf(const T &x) const
	{
		if(x<m_base)
			return 0;
		if(x > (m_base+m_step*(m_elements.size()-1)))
			return 1;
		
		T vi=(x-m_base)/m_step;
		size_t i=floor(vi);
		
		return sum(m_elements.begin(), m_elements.begin()+i+1);
	}

	virtual uint64_t ElementCount() const
	{ return m_elements.size(); }

	// All legal distributions contain at least one element.
	virtual std::pair<T,T> GetElement(uint64_t index) const
	{ return std::make_pair(m_base+index*m_step, m_elements.at(index)); }

	void GetElements(uint64_t begin, uint64_t end, std::pair<T,T> *dest) const
	{
		if((end<begin) || (end>ElementCount()))
			throw std::range_error("Requested elements are out of range.");
		T curr=m_base+m_step*begin;
		const T *src=&m_elements[begin];
		for(int i=0;i<(end-begin);i++){
			dest=std::make_pair(curr, *src);
			curr+=m_step;
			src++;
		}
	}
	
	void GetProbabilities(uint64_t begin, uint64_t end, T *dest) const
	{
		if((end<begin) || (end>ElementCount()))
			throw std::range_error("Requested elements are out of range.");
		
		std::copy(m_elements.begin()+ begin, m_elements.begin()+end, dest);
	}
};

}; // random
}; // flopoco

#endif
