#ifndef random_base_table_distribution_hpp
#define random_base_table_distribution_hpp

#include "distribution.hpp"
#include "moment_conversions.hpp"

#include "random/utils/sum.hpp"

#include <boost/make_shared.hpp>

#include <vector>
#include <algorithm>
#include <numeric>

namespace flopoco
{
namespace random
{

/* This implements a distribution based on a table of (x,p) pairs. It is a bit odd,
	as it internally remembers all the distinct pairs, but when viewed as a distribution
	we have to collapse all identical x values down to a single point.
*/
template<class T>
class TableDistribution
	: public DiscreteDistribution<T>
{
public:	
	typedef boost::shared_ptr<TableDistribution> TypePtr;
private:
	bool m_isSymmetric;
	int m_fracBits;
	T m_zero, m_one;

	// Elements are sorted as (x,p). Repeated x is allowed. Must have 0<=p<=1, so zero probability elements are allowed.
	typedef std::vector<std::pair<T,T> > storage_t;
	storage_t m_elements;
	std::vector<T> m_range, m_pdf, m_cdf;
	
	mutable std::vector<T> m_rawMoments;

	// First sort on x, then sort on decreasing p if x<0 or increasing p if x>0
	static bool EltLessThan(const std::pair<T,T> &a, const std::pair<T,T> &b)
	{
		if(a.first<b.first)
			return true;
		if(a.first>b.first)
			return false;
		if(a.first < 0)
			return a.second<b.second;
		else
			return a.second>b.second;
	};
	
	static bool ProbLessThan(const std::pair<T,T> &a, const std::pair<T,T> &b)
	{
		return a.second<b.second;
	};
	
	void CompleteInit(bool sorted, const T &acc)
	{
		unsigned n=m_elements.size();
		
		if(acc!=m_one){
			T scale=m_one/acc;
			for(size_t i=0;i<m_elements.size();i++){
				m_elements[i].second *=scale;
			}
		}
		
		if(!sorted)
			std::sort(m_elements.begin(), m_elements.end(), EltLessThan);

		if((n%2) && (m_elements[n/2].first!=0)){
			//fprintf(stderr, "elts[n/2]=%lg\n", m_elements[n/2].first);
			m_isSymmetric=false;
		}
		
		if(m_isSymmetric){
			for(size_t i=0;i<n/2 && m_isSymmetric;i++){
				if(m_elements[i].first!=-m_elements[n-i-1].first){
					//fprintf(stderr, "first\n");
					m_isSymmetric=false;
				}
				if(m_elements[i].second!=m_elements[n-i-1].second){
					//fprintf(stderr, "second : e[%d] = %lg, e[%d]=%lg,  diff=%lg\n", i, m_elements[i].second, n-i-1, m_elements[n-i-1].second,  m_elements[i].second-m_elements[n-i-1].second);
					m_isSymmetric=false;
				}
			}
		}
		
		m_range.resize(m_elements.size(), m_zero);
		m_pdf.resize(m_elements.size(), m_zero);
		m_cdf.resize(m_elements.size(), m_zero);
		
		int dest_i=-1;
		T dest_x=m_elements[0].first-1;
		T running_acc=m_zero;
		for(int i=0;i<(int)m_elements.size();i++){
			if(m_elements[i].first!=dest_x){
				dest_i++;
				dest_x=m_elements[i].first;
				m_range[dest_i]=dest_x;
				
				if(m_fracBits!=INT_MAX){
					T tmp=ldexp(round(ldexp(dest_x,m_fracBits)),-m_fracBits);
					if(tmp!=dest_x)
						throw std::string("Table is being created with fixed-point specification, but points are not aligned.");
				}
			}
			T p=m_elements[i].second;
			running_acc += p;
			m_pdf[dest_i] += p;
			m_cdf[dest_i] = running_acc;
		}
		
		m_range.resize(dest_i+1);
		m_cdf.resize(dest_i+1);
		m_pdf.resize(dest_i+1);
	}
public:
	template<class TC>
	TableDistribution(const TC &src, int fracBits=INT_MAX)
		: m_isSymmetric(true)
		, m_fracBits(fracBits)
		, m_zero(src.begin()->second-src.begin()->second)
		, m_one(m_zero+1)
	{
		if(src.size()==0)
			throw std::invalid_argument("TableDistribution - Table must contain at least one element.");
		
		m_elements.assign(src.begin(), src.end());
		T acc=sum_weighted_powers(m_elements.begin(), m_elements.end(), 0);
		if(fabs(acc-1)>1e-12)
			throw std::logic_error("TableDistribution - Element probabilities are more than 1e-12 from one.");
		
		CompleteInit(false, acc);
	}

	TableDistribution(const T *begin, const T *end, int fracBits=INT_MAX)
		: m_isSymmetric(true)
		, m_fracBits(fracBits)
		, m_zero(*begin-*begin)
		, m_one(m_zero+1)
	{
		if(end<=begin)
			throw std::invalid_argument("TableDistribution - Table must contain at least one element.");
		
		size_t n=end-begin;
		
		m_elements.reserve(n);
		T p=m_one;
		p=p/(int)n;
		bool sorted=true;
		for(int i=0;i<(end-begin);i++){
			m_elements.push_back(std::make_pair(begin[i],p));
			if(i!=0)
				sorted=sorted && EltLessThan(*(m_elements.end()-1), m_elements.back());
		}

		CompleteInit(sorted, m_one);
	}
	
	//! If the table has a specific resolution this will be returned, otherwise INT_MAX is returned
	int FixedPointResolution() const
	{
		return m_fracBits;
	}
	
	virtual T RangeGranularity() const
	{
		if(m_fracBits==INT_MAX)
			return m_zero;
		else
			return ldexp(m_one, -m_fracBits);
	}
	
	T RawMoment(unsigned k) const
	{
		if(k>=m_rawMoments.size()){
			while(k>=m_rawMoments.size()){
				int kk=m_rawMoments.size();
				if((kk%2) && m_isSymmetric){
					m_rawMoments.push_back(m_zero);
				}else{
					m_rawMoments.push_back(sum_weighted_powers(m_elements.begin(), m_elements.end(), (int)kk));
				}
			}
		}
		return m_rawMoments[k];
	}

	virtual T StandardMoment(unsigned k) const
	{
		if(k==0)
			return m_one;
		if(k==1)
			return RawMoment(k);
		if(k==2)
			return sqrt(CentralMoment(k));
		return CentralMoment(k) / pow(CentralMoment(2), k/2.0);
	}
	
	T CentralMoment(unsigned k) const
	{
		RawMoment(k);	// force calculate if necessary
		return RawMomentsToCentralMoment(k, &m_rawMoments[0]);
	}
	
	virtual bool IsSymmetric() const
	{ return m_isSymmetric; } // this is more strict than symmetric about the mean, but still is valid
		
	virtual std::pair<T,T> Support() const
	{ return std::make_pair(m_elements.front().first, m_elements.back().first); }
		
	virtual T Pmf(const T &x) const
	{
		// Elements must have non-negative probability, so we always have (x-eps,p) < (x,-1) < (x,p)
		typename std::vector<T>::const_iterator it=std::lower_bound(m_range.begin(), m_range.end(), x);
		if(it==m_range.end())
			return m_zero;
		if(*it!=x)
			return m_zero;
		return m_pdf[it-m_range.begin()];
	}
	
	virtual T Cdf(const T &x) const
	{
		typename std::vector<T>::const_iterator it=std::lower_bound(m_range.begin(), m_range.end(), x);
		if(it==m_range.end())
			return m_one;
		if(*it != x){
			if(it==m_range.begin())
				return m_zero;
			else
				return m_cdf[it-m_range.begin()-1];
		}
		return m_cdf[it-m_range.begin()];
	}

	virtual uint64_t ElementCount() const
	{ return m_elements.size(); }
	
	virtual int64_t IndexFromRange(const T &x) const
	{
		typename std::vector<T>::const_iterator it=std::lower_bound(m_range.begin(), m_range.end(), x);
		if(it==m_range.end())
			throw std::string("IndexFromRange - Value ")+boost::lexical_cast<std::string>(x)+" is out of range.";
		if(*it!=x)
			throw std::string("IndexFromRange - Value ")+boost::lexical_cast<std::string>(x)+" does not occur.";
		return it-m_range.begin();
	}
	
	virtual int64_t ClosestIndexFromRange(const T &x) const
	{
		typename std::vector<T>::const_iterator it=std::lower_bound(m_range.begin(), m_range.end(), x);
		if(it==m_range.end())
			return m_range.size()-1;
		return it-m_range.begin();
	}

	virtual T RangeFromIndex(int64_t x) const
	{
		return m_range.at(x);
	}

	std::vector<std::pair<T,T> > GetElements() const
	{
		return m_elements;
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
