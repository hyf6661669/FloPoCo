#ifndef flopoco_random_mpreal_h
#define flopoco_random_mpreal_h

#include "boost/multiprecision/mpfr.hpp"

namespace mpfr{
	
	typedef boost::multiprecision::mpfr_float mpreal;
	
	inline mpfr_ptr get_mpfr_ptr(mpfr::mpreal &x)
	{ return x.backend().data(); }

	inline mpfr_srcptr get_mpfr_srcptr(const mpfr::mpreal &x)
	{ return x.backend().data(); }

	inline mpreal create_zero(unsigned digits)
	{
		assert(digits>0);

		mpreal res;
		res=0;
		mpfr_set_prec(get_mpfr_ptr(res), digits);
		return res;
	}

	inline mpreal calc_pi(unsigned digits)
	{
		mpreal res=create_zero(digits);
		mpfr_const_pi(get_mpfr_ptr(res), MPFR_RNDN);
		return res;
	}

	inline mpfr_prec_t get_mpfr_prec(const mpreal &x)
	{
		return mpfr_get_prec(get_mpfr_srcptr(x));
	}

	inline mpreal sqr(const mpreal &a)
	{
		mpreal tmp(a);
		mpfr_sqr(get_mpfr_ptr(tmp),get_mpfr_srcptr(a), MPFR_RNDN);
		return tmp;
	}

	inline mpreal log2(const mpreal &a)
	{
		mpreal tmp(a);
		mpfr_log2(get_mpfr_ptr(tmp), get_mpfr_srcptr(a), MPFR_RNDN);
		return tmp;
	}
};

using mpfr::get_mpfr_ptr;
using mpfr::get_mpfr_srcptr;
using mpfr::get_mpfr_prec;
using mpfr::sqr;
using mpfr::log2;


#endif
