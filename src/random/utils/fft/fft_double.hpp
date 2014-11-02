#ifndef flopoco_random_utils_fft_double_hpp
#define flopoco_random_utils_fft_double_hpp

#include <cassert>

#include "mpreal.h"

#include "random/utils/fft/fft.hpp"

#include "gsl/gsl_fft_complex.h"

namespace flopoco
{
namespace random
{
	inline void fft_radix2_complex(double *data, unsigned nn)
	{ gsl_fft_complex_radix2_forward(data, 1, nn); }
	
	inline void ifft_radix2_complex(double *data, unsigned nn)
	{ gsl_fft_complex_radix2_inverse(data, 1, nn); }
};
};

#endif
