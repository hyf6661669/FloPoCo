#ifndef flopoco_random_distributions_build_chi2_partition_hpp
#define flopoco_random_distributions_build_chi2_partition_hpp

#include <algorithm>
#include <numeric>

#include "random/utils/binomial_sampler.hpp"

#include <boost/math/distributions/chi_squared.hpp>

#include <random/utils/mpreal/boost_math_mpreal.hpp>

namespace flopoco{
namespace random{
	
struct Chi2Res
{
	int dof;
	double n;
	double statistic;
	double pvalue;
};

template<class T>
class Chi2Partition
	: public boost::enable_shared_from_this<Chi2Partition<T> >
{
private:
	T m_zero, m_one;
	std::vector<T> m_boundaries;	// contains n+1 boundaries where m_boundary[i] <= bucket[i] < m_boundary[i]
	std::vector<T> m_probabilities;	// probabilities for the n buckets

	void CheckBuckets()
	{
		T sum=std::accumulate(m_probabilities.begin(), m_probabilities.end(), T());
		if(abs(sum-1) > 1e-12)
			throw std::string("Chi2Partition - Sum of bucket probs does not add up to 1.");
		
		m_zero=sum-sum;
		m_one=1+m_zero;
		
		if(m_boundaries.size()-1 != m_probabilities.size())
			throw std::string("Chi2Partition - Boundary and probability sizes do not match.");
		for(int i=1;i<(int)m_boundaries.size();i++){
			if(m_boundaries[i-1] >= m_boundaries[i])
				throw std::string("Chi2Partition - boundaries do not increase monotonically.");
		}
	}
public:
	Chi2Partition(const std::vector<T> &boundaries, const std::vector<T> &probabilities)
		: m_boundaries(boundaries)
		, m_probabilities(probabilities)
	{
		CheckBuckets();
	}
	
	Chi2Partition(const std::vector<T> &boundaries, typename DiscreteDistribution<T>::TypePtr source)
		: m_boundaries(boundaries)
	{
		if(m_boundaries.size()<3)
			throw std::string("Chi2Parition - Need at least three boundaries to get two buckets.");
		
		if(! (source->IsSymmetric() && source->StandardMoment(1)==0))
			throw std::string("Only symmetric zero-mean distributions supported at the moment.");
		
		m_probabilities.resize(m_boundaries.size()-1);
		T prev=source->Cdf(m_boundaries[0]);
		for(int i=0;i<(int)m_probabilities.size()/2;i++){
			T curr=source->Cdf(source->RangePrev(m_boundaries[i+1]));
			m_probabilities[i]=curr-prev;
			prev=curr;
		}
		prev=0;
		for(int i=(int)m_probabilities.size()-1;i>=(int)m_probabilities.size()/2;i--){
			T curr=source->Cdf(-m_boundaries[i]);
			m_probabilities[i]=curr-prev;
			prev=curr;
		}
		
		/*T sum=0;
		for(int i=0;i<m_probabilities.size();i++){
			sum+=m_probabilities[i];
			std::cerr<<m_boundaries[i]<<", "<<m_boundaries[i+1]<<", "<<m_probabilities[i]<<", "<<sum<<"\n";
		}*/
		
		CheckBuckets();
	}
	
	static boost::shared_ptr<const Chi2Partition> CreateEqualProbability(typename DiscreteDistribution<T>::TypePtr source, int n)
	{
		if(! (source->IsSymmetric() && source->StandardMoment(1)==0))
			throw std::string("Only symmetric zero-mean distributions supported at the moment.");
		
		T one=source->Pmf(0);
		one=(one+1)/(one+1);
		T zero=one-one;
		zero=zero;
		
		std::vector<T> boundaries;
		boundaries.push_back(source->Support().first);
		//std::cerr<<0<<", "<<boundaries.back()<<"\n";
		T step=one/n;
		for(int i=1;i<=n/2;i++){
			T p=i*step;
			T x=source->InvCdf(p);
			if(x!=boundaries.back()){
				boundaries.push_back(x);
				//std::cerr<<i<<", "<<boundaries.back()<<"\n";
			}
		}
		//if((n%2)==0)
		//	boundaries.push_back(zero);
		for(int i=boundaries.size()-1;i>0;i--){
			if(boundaries[i]!=0){
				T x=source->RangeFromIndex(source->IndexFromRange(-boundaries[i])+1) ;
				boundaries.push_back(x);
			}
		}
		boundaries.push_back(-boundaries.front());
		
		return boost::make_shared<Chi2Partition>(boundaries, source);
	}
	
	/*! Create a partition with ~n buckets, where the first bucket is [-inf,firstUpperX) and the last bucket is [lastLowerX,+inf) */
	static boost::shared_ptr<const Chi2Partition> CreateEqualRange(typename DiscreteDistribution<T>::TypePtr source, T firstUpperX, T lastLowerX, int n)
	{
		if(! (source->IsSymmetric() && source->StandardMoment(1)==0))
			throw std::string("Only symmetric zero-mean distributions supported at the moment.");
		
		if(n<3)
			throw std::string("Need to create at least three buckets.");
		
		T zero=(firstUpperX-firstUpperX)+(lastLowerX-lastLowerX)+source->Pmf(0)-source->Pmf(0);
		T one=zero+1;
		
		assert(firstUpperX<lastLowerX);
		T step=(one*lastLowerX - one*firstUpperX)/(n-2);
		
		std::vector<T> boundaries;
		boundaries.push_back(source->Support().first);
		
		//std::cerr<<"[-inf,"<<firstUpperX<<"), ... , ["<<lastLowerX<<",+inf), n="<<n<<", step="<<step<<"\n";
		
		T curr=firstUpperX;
		for(int i=1;i<n-1;i++){
			//std::cerr<<i<<" : "<<curr<<"\n";
			T x=source->RangeFromIndex(source->ClosestIndexFromRange(curr));
			if(x!=boundaries.back() && x < lastLowerX)
				boundaries.push_back(x);
			curr+=step;
		}
		boundaries.push_back(lastLowerX);
		boundaries.push_back(source->Support().second);
		
		return boost::make_shared<Chi2Partition>(boundaries, source);
	}
	
	unsigned NumBuckets() const
	{ return m_probabilities.size(); }
	
	void Dump(std::ostream &dst) const
	{
		for(unsigned i=0;i<m_probabilities.size();i++){
			dst<<i<< " : ["<<m_boundaries[i]<<".."<<m_boundaries[i+1]<<") = "<<m_probabilities[i]<<"\n";
		}
	}
	
	T MinChi2SampleSize() const
	{
		T min_p=*std::min_element(m_probabilities.begin(), m_probabilities.end());
		assert(min_p!=0);
		return 10.0 / min_p;
	}
	
	boost::shared_ptr<Chi2Partition<T> > Resample(typename DiscreteDistribution<T>::TypePtr source) const
	{
		return boost::make_shared<Chi2Partition<T> >(m_boundaries, source);
	}
	
	template<class TRng>
	std::vector<T> GenerateSample(T n, TRng &rng) const
	{
		std::vector<T> counts(m_probabilities.size(), m_zero);
		
		//std::cerr<<"Sorting\n";
		std::vector<std::pair<T,int> > reorder;
		for(int i=0;i<(int)m_probabilities.size();i++){
			reorder.push_back(std::make_pair(m_probabilities[i], i));
		}
		std::sort(reorder.begin(), reorder.end());
		
		//std::cerr<<"Generating\n";
		T gone=m_zero;
		for(int i=0;i<(int)reorder.size();i++){
			assert(gone<=1);
			T p=reorder[i].first;
			T bin_p=p / (1-gone);	// Hrmm, not too happy about that
			
			T got;
			if(i==(int)reorder.size()-1){
				if( abs(bin_p-1) > 1e-10)
					throw std::string("GenerateSample - Loss of precision while updating probs. Move to higher precision?");
				got=n;
			}else{
				//std::cerr<<i<<", p="<<p<<", bin_p="<<bin_p<<", n="<<n<<"\n";
				
				got=GenerateBinomialSample<T,TRng>(bin_p, n, rng, true);
				//std::cerr<<"  got="<<got<<"\n";
			}
			
			assert(got <= n);
			assert(round(got)==got);
			counts[reorder[i].second]=got;
			n -= got;
			
			if(n==0)
				break;
			
			gone += p;
		}
		
		assert(n==0);
		
		return counts;
	}
	
	template<class TC>
	Chi2Res Chi2Test(const std::vector<TC> &counts) const
	{
		if(counts.size()!=m_probabilities.size())
			throw std::string("Chi2Test - Wrong number of counts.");
		
		T sum=std::accumulate(counts.begin(), counts.end(), m_zero);
		//std::cerr<<"Sum="<<sum<<"\n";
		
		T acc=T();
		for(int i=0;i<(int)m_probabilities.size();i++){
			T o=counts[i], e=m_probabilities[i]*sum;
			//std::cerr<<"o="<<o<<", e="<<e<<"\n";
			if(e < 10)
				throw std::string("Chi2Test - One bucket has expected count < 10");
			T diff=o-e;
			acc += diff*diff / e;
		}
		
		Chi2Res res;
		res.dof=m_probabilities.size()-1;
		res.n=boost::math::tools::real_cast<double>(sum);
		res.statistic=boost::math::tools::real_cast<double>(acc);
		res.pvalue=cdf(complement(boost::math::chi_squared_distribution<double>(m_probabilities.size()-1), res.statistic)); 
		return res;
	}
	
	boost::shared_ptr<const Chi2Partition> AdaptForSampleSize(T n) const
	{
		T lowest=*std::min_element(m_probabilities.begin(), m_probabilities.end());
		if(lowest * n > 10)
			return boost::enable_shared_from_this<Chi2Partition<T> > ::shared_from_this();
		
		
		std::vector<T> boundaries;
		std::vector<T> probabilities;
		boundaries.push_back(m_boundaries[0]);
		
		T curr=m_zero;
		
		for(int i=0;i<(int)m_probabilities.size();i++){
			curr += m_probabilities[i];
			if( curr * n > 10){
				probabilities.push_back(curr);
				boundaries.push_back(m_boundaries[i+1]);
				curr=m_zero;
			}
		}
		if(curr>0){
			probabilities.back()+=curr;
			boundaries.back()=m_boundaries.back();
		}
		
		assert(boundaries.back()==m_boundaries.back());
		
		
		return boost::make_shared<Chi2Partition>(boundaries, probabilities);
	}
};
	
}; // random
}; // flopoco

#endif
