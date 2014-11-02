#include "fixed_format_t.hpp"

#include <assert.h>
#include <sollya.h>

#include <boost/lexical_cast.hpp>

namespace flopoco{
namespace random{
	
bool operator==(const fixed_format_t &a, const fixed_format_t &b)
{
	return (a.isSigned==b.isSigned) && (a.msb==b.msb) && (a.lsb==b.lsb);
}

mpfr::mpreal DecodeRaw(const fixed_format_t &fmt, mpz_class raw, int prec)
{
	if((raw<0) || ((int)mpz_sizeinbase(raw.get_mpz_t(),2) > fmt.width()))
		throw std::string("DecodeRaw - raw value out of range.");
	
	if(fmt.isSigned && mpz_tstbit(raw.get_mpz_t(), fmt.width()-1)){
		raw=raw - (mpz_class(1)<<fmt.width());
	}
	
	mpfr::mpreal res(0.0, prec);
	if(0!=mpfr_set_z(get_mpfr_ptr(res), raw.get_mpz_t(), MPFR_RNDN))
		throw std::string("DecodeRaw - Loss of precision.");
	
	mpfr_mul_2si(get_mpfr_ptr(res), get_mpfr_ptr(res), fmt.lsb, MPFR_RNDN);
	
	return res;
}

mpz_class EncodeRaw(const fixed_format_t &fmt, mpfr::mpreal x)
{
	mpfr_mul_2si(get_mpfr_ptr(x), get_mpfr_ptr(x), -fmt.lsb, MPFR_RNDN);
	
	mpz_class res;
	mpfr_get_z(res.get_mpz_t(), get_mpfr_ptr(x), MPFR_RNDN);
	
	mpz_class maxVal=mpz_class(1)<<fmt.width();
	
	if(res < 0)
		res += maxVal;
	
	if((res<0) || (res>=maxVal))
		throw std::string("EncodeRaw - Raw value out of range.");
	
	return res;
}

fixed_format_t ParseFixedFormat(const std::string &x)
{
	std::string left=x;
	if(left.size()<5)
		throw std::string("ParseFixedFormat('"+x+"') - Need at least five characters.");
	
	fixed_format_t res;
	
	if(left[0]=='U'){
		res.isSigned=false;
	}else if(left[0]=='S'){
		res.isSigned=true;
	}else{
		throw std::string("ParseFixedFormat('"+x+"') - Format should start with 'U' or 'S'.");
	}
	
	if(left[1]!=';')
		throw std::string("ParseFixedFormat('"+x+"') - 'U' or 'S' must be followed by ';'");
	
	left=left.substr(2,-1);
	
	int split=left.find(';');
	if(split==-1)
		throw std::string("ParseFixedFormat('"+x+"') - Couldn't find second ';'");
	
	std::string msbStr=left.substr(0,split);
	std::string lsbStr=left.substr(split+1, -1);
	
	if(msbStr.size()==0)
		throw std::string("ParseFixedFormat('"+x+"') - Couldn't convert msb to a number");
	res.msb=boost::lexical_cast<int>(msbStr);	
	
	if(lsbStr.size()==0)
		throw std::string("ParseFixedFormat('"+x+"') - Couldn't convert lsb to a number");
	res.lsb=boost::lexical_cast<int>(lsbStr);	
	
	return res;
}

}; // random
}; // flopoco
