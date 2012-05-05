#ifndef flopoco_random_distributions_sum_distribution_hpp
#define flopoco_random_distributions_sum_distribution_hpp

namespace flopoco
{
namespace random
{
	template<class T>
	EnumerableDistribution<T>::TypePtr ConvolveDistributions(
		EnumerableDistribution<T>::TypePtr a,
		EnumerableDistribution<T>::TypePtr b,
		const TF &f
	){
		if(a->ElementCount() * b->ElementCount() > pow(2.0,40))
			throw std::runtime_error("Attempt to perform additive convolution with more than 2^40 combinations.");
		
		std::vector<std::pair<T,T> > va, vb;
		a->GetElements(0, a->ElementCount(), &va[0]);
		b->GetElements(0, b->ElementCount(), &vb[0]);
				
		boost::unordered_map<T,T> acc;
		
		for(size_t ia=0;ia<va.size();ia++){
			std::pair<T,T> ca=va[ia];
			if(ca.second==0)
				continue;
			for(size_t ib=0;ib<vb.size();ib++){
				std::pair<T,T> cb;
				acc[ f(ca.first,c.first) ] += ca.second * cb.second;
			}
		}
	}
	
	namespace 
	
	template<class T>
	EnumerableDistribution<T>::TypePtr AddDistributions(
		EnumerableDistribution<T>::TypePtr a,
		EnumerableDistribution<T>::TypePtr b
	){
		if(a->ElementCount() < b->ElementCount)
			return AddDistributions(b,a);
		return ConvolveDistributions(a,b,std::plus<T>());
	}
	
};
};

#endif
