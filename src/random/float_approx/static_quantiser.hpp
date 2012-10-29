#ifndef flopoco_random_float_approx_static_quantiser_hpp
#define flopoco_random_float_approx_static_quantiser_hpp

#include "Operator.hpp"
#include "FixedPointFunctions/Function.hpp"

#include "random/utils/mpfr_vec.hpp"

namespace flopoco
{
namespace random
{

class StaticQuantiser	: public Operator
{
private:	
	int m_wValue;
	std::vector<mpz_class> boundaries;

	void Build(int wValue);	// Given fixed-point boundaries, build the operator
	void AddLevel(int level);
	mpz_class GetBoundary(int level, int index);
public:
	/*! Given an input value x, and a set of n+1 bucket boundaries
		find 0<i<n where boundaries[i] <=x < boundaries[i+1].
		Output has log2ceil(n) buckets
	*/
	StaticQuantiser(Target *target, int wValue, const std::vector<mpz_class> &_boundaries, map<string, double> inputDelays = emptyDelayMap);

	~StaticQuantiser();

	void emulate(TestCase * tc);

	static Operator *BuildFloatQuantiser(Target *target, int wE, int WF, int n, const Function *f, map<string, double> inputDelays = emptyDelayMap);
};

void StaticQuantiser_registerFactory();

}; 
};

#endif
