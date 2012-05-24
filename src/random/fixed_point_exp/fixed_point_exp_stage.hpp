#ifndef flopoco_random_fixed_point_exp_stage_hpp
#define flopoco_random_fixed_point_exp_stage_hpp

#include "Operator.hpp"
#include "utils.hpp"

namespace flopoco
{
namespace random
{

class FixedPointExpStage
	: public flopoco::Operator
{
public:
	typedef double value_t;	

	struct iterate_range_t
	{
		value_t current, delta;
		
		iterate_range_t(value_t _current, value_t _delta)
			: current(_current)
			, delta(_delta)
		{}
		
		bool operator==(const iterate_range_t &o) const
		{ return current==o.current; }
		
		iterate_range_t &operator++()
		{ current+=delta; return *this; }
		
		const value_t &operator*() const
		{ return current; }
	};
	
	struct residual_type
	{
		value_t expMu;
		value_t expSigma;
		bool isSigned;
		int msb;	// power of the highest non-sign bit
		int lsb;		// power of the smallest bit
		value_t rangeMin;
		value_t rangeMax;
		
		// msb==lsb -> single bit wide, msb==lsb-1 -> zero bits
		int Width() const
		{ return (msb-lsb+1) + (isSigned?1:0);
			
		bool IsInRange(value_t x) const
		{ return (rangeMin<=x) && (x<=rangeMax); }
		
		uint64_t ToBits(const value_t &x) const;
		value_t FromBits(uint64_t x) const;
		
		typedef iterate_range_t iterator;
		
		iterate_range_t begin_msbs(int msbs) const
		{
			assert(msbs <= Width());
			value_t delta=power2(msb-msbs);
			value_t curr=isSigned ? -power2(msb) : 0;
			return iterate_range_t(curr, delta);
		}
			
		iterate_range_t end_msbs(int msbs) const
		{
			assert(msbs <= Width());
			value_t delta=power2(msb-msbs);
			return iterate_range(power2(msb), delta);
		}
	};
	
	
	
	value_t RoundFloat(const value_t &x, int bits) const
	{
		assert(bits>=1);
		int e;
		value_t f=frexp(x,&e); // in range [0.5,1)
		return ldexp(f, round(ldexp(f,bits)), e-bits);
	}

protected:
	residual_type m_inputResidualType;
	result_type m_inputResultType;

	FixedPointExpStage(Target* target,
		residual_type inputResidualType,
		result_tye inputResultType
	)
		: m_inputResidualType(inputResidualType)
		, m_inputResultType(inputResultType)
	{}
	
public:
    virtual ~FixedPointExpStage()
	{};
	
	value_t ReferenceExp(const value_t &x) const
	{ return exp(m_expMu+m_expSigma*x); }
	
	const residual_type &InputResidualType() const
	{ return m_inputResidualType; }
	
	result_type InputResultType() const
	{ return m_inputResultType; }
	
	virtual residual_type OutputResidualType() const =0;
	virtual result_type OutputResultType() const=0;
	
	virtual std:pair<value_t,value_t> Execute(const value_t &residual, const value_t &result) const=0;
};

}; // random
}; // flopoco

#endif
