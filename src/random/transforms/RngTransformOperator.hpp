#ifndef random_transform_rng_transform_operator_hpp
#define random_transform_rng_transform_operator_hpp

#include "Operator.hpp"

namespace flopoco
{
namespace random
{

class RngTransformOperator
	: public flopoco::Operator
{
protected:
	RngTransformOperator(Target* target, map<string, double> inputDelays = emptyDelayMap)
		: Operator(target, inputDelays)
	{}
public:
	virtual unsigned uniformInputBits() const=0;
	virtual std::string uniformInputName() const =0;
};

};
};

#endif
