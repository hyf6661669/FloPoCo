#ifndef random_float_approx_hpp
#define random_float_approx_hpp

#include "Operator.hpp"

#include "utils.hpp"

namespace flopoco
{
namespace random
{
  
Operator *CreateFloatApproxOperator(
    Target* target,
    int wDomainE, int wDomainF, mpfr_t domainMin, mpfr_t domainMax,
    int wRangeE, int wRangeF,
    const Function &f,
    unsigned degree
    );

void FloatApprox_registerFactory();

};  // random
};  // flopoco

#endif
