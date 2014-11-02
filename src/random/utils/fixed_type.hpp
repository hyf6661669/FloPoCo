#ifndef flopoco_random_utils_fixed_type_hpp
#define flopoco_random_utils_fixed_type_hpp

#include "mpfr.h"
#include "gmpxx.h"

#include <boost/smart_ptr.hpp>
#include <sstream>

namespace flopoco
{
namespace random
{

template<class TInt>
class iops
{
public:
	static unsigned NumBits(const TInt &x, bool isSigned);
	static bool ToMpfr(mpfr_t res, const TInt &x, int exp=0);	// res=float(x*2^exp)
	static bool FromMpfr(TInt &res, const mpfr_t x, int exp=0);	// res=round(x*2^exp)
	static int GetTwosComplementBit(const TInt &x, int bit);
};

template<>
class iops<mpz_class>
{
public:
	static unsigned NumBits(const mpz_class &x, bool isSigned)
	{
		if(!isSigned){
			assert(x>=0);
			return mpz_sizeinbase(x.get_mpz_t(),2);
		}else{
			if(x==0)
				return 0;
			if(x<0){
				mpz_class tmp(1-x);
				return 1+mpz_sizeinbase(tmp.get_mpz_t(),2);	// For -2^{k} <= x < -2^{k-1} need k+1 bits  -> 2^{k-1} < -x <= 2^k -> 2^{k-1} <= -x+1 = 1-x < 2^k
			}
			return 1+mpz_sizeinbase(x.get_mpz_t(),2);
		}
	}
	
	static bool FromMpfr(mpz_class &res, const mpfr_t src, int exp=0)
	{
		if(exp==0){
			return 0!=mpfr_get_z(res.get_mpz_t(), src, MPFR_RNDN);
		}else{
			mpfr_t tmp;
			mpfr_init2(tmp, mpfr_get_prec(src));
			mpfr_mul_2si(tmp, src, exp, MPFR_RNDN);
			bool rounded=0!=mpfr_get_z(res.get_mpz_t(), src, MPFR_RNDN);
			mpfr_clear(tmp);
			return rounded;
		}
	}
	
	static bool ToMpfr(mpfr_t res, const mpz_class &src, int exp=0)
	{ return 0!=mpfr_set_z_2exp(res, src.get_mpz_t(), exp, MPFR_RNDN); }

	static int GetTwosComplementBit(const mpz_class &x, int bit)
	{ return mpz_tstbit(x.get_mpz_t(), bit); }
};


/* The data-type consists of a type specification
	(number of bits, shift, and whether it is signed or not),
	and a range specification (the range of actual numbers).
	The type specification will always encompass the range
	specification, but the range may be much smaller than the
	type. e.g. the type could by uint8, but the range might only
	be [3,7].

	The different types of value that are understood are:
	- Real : A floating-point representation of a number
	- Raw : An integer value with the same LSB position as the type
	- Fixed : An integer value with a potentially different LSB position (specified seperately).
	
	Conversion functions return true if it was not possible to exactly convert
	values.	
*/
class FixedType
{
public:
	typedef mpz_class integer_t;
private:
	typedef iops<integer_t> io;

	// This is in a holding class as I originally supported segmented ranges
	// but decided it wasn't needed. For now it is just premature optimisation...
	struct info_t
	{
		friend class FixedType;
		
		info_t(int lsb, const integer_t &a, const integer_t &b)
			: m_isSigned(a<0)
			, m_width(std::max(io::NumBits(a,m_isSigned), io::NumBits(b, m_isSigned)))
			, m_lsb(lsb)
		{
			UpdateMask();
		}
		
		info_t(bool isSigned, int width, int lsb)
			: m_isSigned(isSigned)
			, m_width(width)
			, m_lsb(lsb)
		{
			if(width>0){
				if(isSigned){
					m_minRaw=-1;
					m_minRaw<<=(width-1);
					m_maxRaw=(-m_minRaw)-1;
				}else{
					m_maxRaw=1;
					m_maxRaw=(m_maxRaw<<width)-1;
				}
				UpdateMask();
			}
		}
		
		void UpdateMask()
		{
			m_mask=1;
			m_mask<<=m_width;
			m_mask-=1;
		}
		
		// This is always a widening conversion. If the sign/width/lsb specify something that
		// is larger than the min/max specify, then they won't be changed.
		void ExtendRangeRaw(const mpz_class &x)
		{
			if((x<m_minRaw) || (m_maxRaw<x)){
				m_minRaw=std::min(m_minRaw, x);
				m_maxRaw=std::max(m_maxRaw, x);
				
				assert(m_lsb!=INT_MAX);	// Somebody must have specified an lsb already
				if((m_minRaw < 0) && !m_isSigned){
					m_isSigned=true;
					m_width++;	// Make sure we don't narrow the unsigned part
				}
				m_width = std::max((unsigned)m_width, std::max(io::NumBits(m_minRaw, m_isSigned), io::NumBits(m_maxRaw, m_isSigned)) );
				UpdateMask();
			}
		}
		
		void ExtendRangeFixed(int lsb, const integer_t &x)
		{
			if(m_lsb==INT_MAX){
				m_lsb=lsb;
				ExtendRangeRaw(x);
			}else if(lsb < m_lsb){
				int diff=m_lsb-lsb;
				m_minRaw <<= diff;
				m_maxRaw <<= diff;
				m_width+=diff;	// Don't narrow the type
				m_lsb=lsb;
				ExtendRangeRaw(x);
			}else{
				ExtendRangeRaw( x<<(lsb-m_lsb) );
			}
		}		
		
		bool MaskToType(integer_t &x) const
		{
			integer_t orig(x);
			x &= m_mask;
			return (x!=orig) || (x < m_minRaw) || (x>m_maxRaw);
		}
		
		bool m_isSigned;
		int m_width;
		int m_lsb;
		mpz_class m_minRaw, m_maxRaw;
		mpz_class m_mask;
	};
	typedef boost::shared_ptr<info_t> info_ptr_t;
	
	void MakeUnique()
	{
		if(!m_p.unique())
			m_p=boost::make_shared<info_t>(*m_p);
	}
	
	info_ptr_t m_p;
public:
	FixedType()
		: m_p(boost::make_shared<info_t>(false, 0, INT_MAX))
	{}
	
	FixedType(bool isSigned, int width, int lsb)
		: m_p(boost::make_shared<info_t>(isSigned, width, lsb))
	{
		if(width<0)
			throw std::invalid_argument("FixedType - width cannot be zero.");
	}
	
	FixedType(int lsb, integer_t minRaw, integer_t maxRaw)
		: m_p(boost::make_shared<info_t>(lsb, minRaw, maxRaw))
	{}
	
	bool IsSigned() const
	{ return m_p->m_isSigned; }
	
	int Width() const
	{ return m_p->m_width; }
	
	int LSB() const
	{ return m_p->m_lsb; }
	
	int MSB() const
	{ return LSB()+Width() - (IsSigned()?1:0); }
	
	
	const integer_t &MinRaw() const
	{ return m_p->m_minRaw; }
	
	const integer_t &MaxRaw() const
	{ return m_p->m_maxRaw; }	
	
	bool MinReal(mpfr_t res) const
	{ return RealFromRaw(res, MinRaw()); }
	
	bool MaxReal(mpfr_t res) const
	{ return RealFromRaw(res, MaxRaw()); }

	bool RawFromReal(integer_t &dst, const mpfr_t src) const
	{ return io::FromMpfr(dst, src, -LSB());	}
	
	bool RealFromRaw(mpfr_t dst, const integer_t &src) const
	{ return io::ToMpfr(dst, src, LSB());	}
	
	bool RawFromFixed(integer_t &dst, int lsb, integer_t src) const
	{
		int diff=lsb-LSB();
		if(diff > 0){
			dst=src<<diff;
			return false;
		}else{
			dst=src>>(-diff);
			return (dst<<(-diff)) != src;
		}
	}
	
	bool FixedFromRaw(integer_t &dst, int lsb, integer_t src) const
	{
		int diff=LSB()-lsb;
		if(diff > 0){
			dst=src<<diff;
			return false;
		}else{
			dst=src>>(-diff);
			return (dst<<(-diff)) != src;
		}	
	}
	
	bool IsInRangeRaw(const integer_t &x) const
	{ return (MinRaw() >= x) && (MaxRaw() <= x); }
	
	// Note that a number which is not representable is considered to be out of range
	bool IsInRangeFixed(int lsb, const integer_t &x) const
	{
		if(LSB()==INT_MAX){
			assert(MinRaw()==0 && MaxRaw()==0);
			return x==0;	// If there is no type, then only zero is in range
		}else{
			integer_t tmp;
			if(RawFromFixed(tmp, lsb, x))
				return false;
			return IsInRangeRaw(tmp);
		}
	}
	
	void ExtendRangeRaw(const integer_t &x)
	{
		if(IsInRangeRaw(x))
			return;
		MakeUnique();
		m_p->ExtendRangeRaw(x);
	}
	
	void ExtendRangeFixed(int lsb, const integer_t &x)
	{
		if(IsInRangeRaw(x))
			return;
		MakeUnique();
		m_p->ExtendRangeFixed(lsb, x);
	}
	
	bool ExtendRangeReal(int lsb, mpfr_t x)
	{
		integer_t tmp;
		bool rounded=io::FromMpfr(tmp, x, -lsb);
		if(IsInRangeFixed(lsb, tmp))
			return rounded;
		MakeUnique();
		m_p->ExtendRangeFixed(lsb, tmp);
		return rounded;
	}
	
	// Uses the current lsb as the lsb of the input
	bool ExtendRangeReal(mpfr_t x)
	{
		assert(LSB()!=INT_MAX);
		return ExtendRangeReal(LSB(), x);
	}
	
	FixedType GetEnclosingType() const
	{ return FixedType(LSB(), MinRaw(), MaxRaw()); }	// Could be more efficient
	
	std::string VHDLType() const
	{
		std::stringstream acc;
		acc<<(IsSigned()?"signed":"unsigned")<<"("<<Width()-1<<" downto 0)";
		return acc.str();
	}
	
	template<class TInt>
	std::string VHDLConstant(const TInt &value) const
	{
		std::stringstream acc;
		acc<<"\"";
		for(int i=Width()-1;i>=0;i--){
			acc<< io::GetTwosComplementBit(value, i);
		}
		acc<<"\"";
		return acc.str();
	}
	
	enum rounding_mode_t{
		ROUND_TRUNCATE,
		ROUND_HALF_UP
	};
	
	
	/* Truncate the incoming fixed-point number to this format
		If the output is out of the types range, then the result will behave as if
		pure bit truncation has happened, and the function will return true. If
		the output out of the specified range (but within the type), then the result
		will be correct, but it will still return true
	*/
	bool RoundTruncateFixed(integer_t &dst, int lsb, integer_t src) const
	{
		if(lsb<LSB()){
			dst=src>>(LSB()-lsb);
		}else{
			dst=src<<(lsb-LSB());
		}
		return m_p->MaskToType(dst);
	}
	
	/* As with truncate fixed, but use the method floor(x+0.5) for the conversion, i.e.
		truncate to (m_width+1), add 1, then truncate to m_width,
		or equivalently  ( src[LSB()-1] ? src+1 : src )
	*/
	bool RoundHalfUpFixed(integer_t &dst, int lsb, integer_t src) const
	{
		if(lsb<LSB()){	// The input has more precision
			int places=LSB()-lsb;	// the number of digits to lose
			dst=src>>(places-1);	// drop all but one of the digits (the shift may be zero)
			dst++;	// Now add 0.5
			dst=src>>1;	// chop off the final digit
		}else{
			dst=(src<<(lsb-LSB()));	// the precision is equal or lower
		}
		return m_p->MaskToType(dst);
	}
	
	
	bool RoundFixed(integer_t &dst, int lsb, integer_t src, rounding_mode_t mode) const
	{
		if(mode==ROUND_TRUNCATE){
			return RoundTruncateFixed(dst, lsb, src);
		}else if(mode==ROUND_HALF_UP){
			return RoundHalfUpFixed(dst, lsb, src);
		}else{
			throw std::invalid_argument("RoundFixed - Unknown rounding mode.");
		}
	}
};

}; // random
}; // flopoco

#endif
