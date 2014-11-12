#ifndef flopoco_random_poly_fixed_format_t_hpp
#define flopoco_random_poly_fixed_format_t_hpp

#include "mpreal.h"
#include "gmpxx.h"

namespace flopoco
{
namespace random
{
	
/*! Represents fixed-point types as (isSigned,msb,lsb)
	The type has bits 2^msb..2^lsb inclusive, so width is always (msb-lsb+1).
	isSigned indicates if msb is negative or not.
	The lsb is dominant, so zero-width (if necessary) has msb=lsb-1.
*/
struct fixed_format_t{
	bool isSigned;
	int msb;
	int lsb;
	int width() const { return msb-lsb+1; }
};

bool operator==(const fixed_format_t &a, const fixed_format_t &b);

fixed_format_t ParseFixedFormat(const std::string &spec);

mpfr::mpreal DecodeRaw(const fixed_format_t &fmt, mpz_class raw, int prec);
mpz_class EncodeRaw(const fixed_format_t &fmt, mpfr::mpreal raw);

// Return the smallest format which can encode the given value
fixed_format_t FixedFormatFromValue(const mpfr::mpreal &x);

//! Extend the current type to include the specified value
fixed_format_t Union(const fixed_format_t &curr, const mpfr::mpreal &x);

std::ostream &operator<<(std::ostream &dst, const fixed_format_t &fmt);

}; // random
}; // flopoco

#endif
