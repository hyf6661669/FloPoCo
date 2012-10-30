#ifndef random_utils_conversions_hpp
#define random_utils_conversions_hpp

#ifdef HAVE_NTL
#include "NTL/quad_float.h"
#endif

#include "mpreal.h"

namespace flopoco
{
namespace random
{
	template<class T>
	T convert(const double x);
	
	template<>
	double convert<double>(const double x)
	{ return x; }
	
#ifdef HAVE_NTL
	
	template<class T>
	T convert(const NTL::quad_float &x);
	
	template<>
	NTL::quad_float convert<NTL::quad_float>(const NTL::quad_float &x)
	{ return x; }
	
	template<>
	NTL::quad_float convert<NTL::quad_float>(const double x)
	{ return NTL::to_quad_float(x); }
	
	template<>
	double convert<double>(const NTL::quad_float &x)
	{ return NTL::to_double(x); }
#endif
	
	template<class T>
	T convert(const mpfr::mpreal &x);
	
	template<>
	double convert<double>(const mpfr::mpreal &x)
	{ return x.toDouble(); }
	
	template<>
	mpfr::mpreal convert<mpfr::mpreal>(const mpfr::mpreal &x)
	{ return x; }
	
	template<>
	mpfr::mpreal convert<mpfr::mpreal>(const double x)
	{ return mpfr::mpreal(x); }
	
}; // random
}; // flopoco

#endif
