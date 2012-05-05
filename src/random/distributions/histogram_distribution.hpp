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

template<class T>
class HistogramDistribution
	: public EnumerableDistribution<T>
{
public:	
	typedef boost::shared_ptr<HistoramDistribution> TypePtr;
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
		
		std::sort(m_elements.begin(), m_elements.end(), ProbLessThan);
		T acc=std::sum(m_elements.begin(), m_elements.end(), (T)0.0);
		
		if(fabs(acc-1)>1e-12)
			throw std::logic_error("TableDistribution - Element probabilities are more than 1e-12 from one.");
		T scale=1.0/acc;
		for(unsigned i=0;i<src.size();i++){
			if(m_elements[i] < 0){
				throw std::logic_error("TableDistribution - Negative element probability.");
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
	
	virtual T StandardMoment(unsigned k) const
	{
		while(m_standardMoments.size()<=k){
			std::vector<T> tmp(m_elements.size());
			
			unsigned kc=m_standardMoments.size();
			if((kc%2) && m_isSymmetric){
				m_standardMoments.push_back(0);
			}else{
				// Lazy...
				T curr=base;
				for(unsigned i=0;i<m_elements.size();i++){
					tmp[i] = pow(curr, kc) * m_elements[i];
					curr += step;
				}				
				std::sort(tmp.begin(), tmp.end());
				T ss=std::accumulate(tmp.begin(), tmp.end(), (T)0.0);
			}
		}
	}
	
	virtual bool IsSymmetric() const
	{ return m_isSymmetric; } // this is more strict than symmetric about the mean, but still is valid
		
	virtual std::pair<T,T> Support() const
	{ return std::make_pair(m_elements.front().first, m_elements.back().second); }
		
	virtual T Pmf(const T &x) const
	{
		// Elements must have non-negative probability, so we always have (x-eps,p) < (x,-1) < (x,p)
		typename storage_t::const_iterator it=std::lower_bound(m_elements.begin(), m_elements.end(), std::pair<T,T>(x,-1.0));
		T acc=0.0;
		while( it!=m_elements.end() ? it->first==x : false ){
			acc+=it->second;
			++it;
		}
		return acc;
	}
	
	virtual T Cdf(const T &x) const
	{
		if(x>=m_elements.back().first)
			return 1.0;
		T acc=0.0;
		for(unsigned i=0;i<m_elements.size();i++){
			if(x<m_elements[i].first)
				return acc;
			acc+=m_elements[i].second;
		}
		assert(0);
	}

	virtual uint64_t ElementCount() const
	{ return m_elements.size(); }

	// All legal distributions contain at least one element.
	virtual std::pair<T,T> GetElement(uint64_t index) const
	{ return m_elements.at(index); }

	void GetElements(uint64_t begin, uint64_t end, std::pair<T,T> *dest) const
	{
		if((end>begin) || (end>=ElementCount()))
			throw std::range_error("Requested elements are out of range.");
		std::copy(m_elements.begin()+begin, m_elements.end()+end, dest);
	}
	
	TypePtr ApplyPolynomial(const std::vector<T> &poly) const
	{
		storage_t elts(m_elements.begin(), m_elements.end());
		for(unsigned i=0;i<elts.size();i++){
			T xx=elts[i].first;
			T acc=poly.back();
			for(int j=(int)poly.size()-2;j>=0;j--){
				acc=poly[j] + acc * xx;
			}
			elts[i].first=acc;
		}
		return boost::make_shared<TableDistribution>(elts);
	}
};

}; // random
}; // flopoco

#endif
