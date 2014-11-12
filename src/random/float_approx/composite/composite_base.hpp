#ifndef flopoco_random_float_approx_composite_base_hpp
#define flopoco_random_float_approx_composite_base_hpp

#include "Operator.hpp"

namespace flopoco
{
namespace random
{
namespace float_approx
{

class CompositeBase
	: public Operator
{
protected:
	// Must inherit to use

	std::string EncodeFPConstant(std::string constant, int wE, int wF);

    std::string EncodeFPConstant(double constant, int wE, int wF);

	void MakeFPNegate(std::string dst, std::string src);
	
	// Output will be "name", instance will be "name"+"_inst"
	Operator *MakeFPExp(std::string name, std::string src);

    Operator *MakeFPSquarer(std::string name, std::string src);

	// Output will be "name", instance will be "name"+"_inst"
    Operator *MakeFPMultiplier(std::string name, std::string srcA, std::string srcB);
	
	Operator *MakeFPConstMult(std::string name, std::string src, std::string constant);

	Operator *MakeFPInverse(std::string name, std::string src);

	Operator *MakeFPConstAdd(std::string name, std::string src, std::string constant);

	CompositeBase(
			Target *target
		);
	
	// No public members
};

}; // float_approx
}; // random
}; // flopoco

#endif
