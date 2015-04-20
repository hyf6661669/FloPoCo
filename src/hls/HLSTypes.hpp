#ifndef flopoco_hls_types_hpp
#define flopoco_hls_types_hpp

#include <memory>
#include <cassert>

#include <Operator.hpp>

namespace flopoco
{

class HLSType;
typedef std::shared_ptr<HLSType> HLSTypePtr;

/*
 * So the type heirarchy is:
 *
 * Type
 * +-Void  : A null type that cannot be instantiated. Used for declarations with no
 * |         associated run-time value, like calls to operators
 * |
 * +-Fixed : Any unsigned or twos complement type.
 * |         isFixed() == true
 * |         isSigned() - Is it unsigned or twos complement
 * |         getMSB() - The weight of the MSB (weight of the sign bit in signed)
 * |         getLSB() - The weight of the LSB
 * |         |
 * |         +- Int - Backward compatibility, getLSB()==0
 * |         |
 * |         +- Bool - Backwards compatibility, getLSB()==0, getWidth()==0
 * |
 * +-Float : Flopoco style float with three fields
 *           isFloat() == true
 *           getExponentWidth()
 *           getFractionWidth()
 *
 * The fixed type has certain sub-types determined by type properties
 * isInt : This is an unsigned or twos complement integer, isFixed() && getLSB()==0
 * isUnsigned : This is an unsigned integer,  isInt() && !isSigned()
 * isSigned : This is a twos complement integer, isInt() && isSigned()
 * isBool : This is a single-bit integer, isInt() && getWidth()==1
 */


class HLSTypeFixed;

class HLSType
{
private:
	int m_width;
protected:
	HLSType(int width)
		: m_width(width)
	{}
public:
	virtual ~HLSType()
	{}
	
	int getWidth() const
	{ return m_width; }

	virtual bool equals(const HLSTypePtr &o) const=0;

	virtual std::string getName() const =0;

	//! Is this an unsigned or twos complement fixed-point number?
	virtual bool isFixed() const
	{ return false; }

	virtual std::shared_ptr<HLSTypeFixed> asFixed() const
	{
		assert(!isFixed());
		throw std::runtime_error("asFixed - this is not a fixed-point type.");
	}

	//! Is this an unsigned or twos complement integer number?
	virtual bool isInt() const
	{ return false; }

	//! Can this represent signed values?
	virtual bool isSigned() const
	{ return false; }

	//! Is this a single-bit logical value?
	virtual bool isBool() const
	{ return false; }

	//! Is this a flopoco float?
	virtual bool isFloat() const
	{ return false; }

	virtual std::shared_ptr<HLSTypeFixed> asFloat() const
	{
		assert(!isFloat());
		throw std::runtime_error("asFixed - this is not a float type.");
	}
};

class HLSTypeVoid
  : public HLSType
{
public:
  HLSTypeVoid()
    : HLSType(0)
  {}
	
  virtual bool equals(const HLSTypePtr &o) const
  { return !!std::dynamic_pointer_cast<HLSTypeVoid>(o); }

  virtual std::string getName() const
  { return "Void"; }

  static std::shared_ptr<HLSTypeVoid> create()
  { return std::make_shared<HLSTypeVoid>(); }
};


class HLSTypeFixed
	: public HLSType
{
private:
	bool m_isSigned;
	int m_msb, m_lsb;
public:
	HLSTypeFixed(bool _signed, int msb, int lsb)
		: HLSType(msb-lsb+1)
		, m_isSigned(_signed)
		, m_msb(msb)
		, m_lsb(lsb)
	{}

	bool isSigned() const
	{ return m_isSigned; }

	int getMSB() const
	{ return m_msb; }

	int getLSB() const
	{ return m_lsb; }

	virtual bool isInt() const
	{ return m_lsb==0; }

	virtual bool isBool() const
	{ return isInt() && getWidth()==1; }

	virtual bool equals(const HLSTypePtr &o) const
	{
		auto p=std::dynamic_pointer_cast<HLSTypeFixed>(o);
		if(!p)
			return false;
		return p->isSigned()==isSigned() && p->getLSB()==getLSB() && p->getMSB()==getMSB();
	}

	virtual std::string getName() const
	{
		std::stringstream tmp;
		tmp<<"Int<signed="<<(isSigned()?"1":"0")<<",width="<<getWidth()<<">";
		return tmp.str();
	}

	static std::shared_ptr<HLSTypeFixed> create(bool _signed, int msb, int lsb)
	{ return std::make_shared<HLSTypeFixed>(_signed, msb, lsb); }
};

class HLSTypeInt
	: public HLSTypeFixed
{
public:
	HLSTypeInt(bool _signed, int w)
		: HLSTypeFixed(_signed, w-1, 0)
	{}

	mpz_class minValue() const
	{
		if(isSigned()){
			mpz_class tmp(1);
			tmp=tmp<<(getWidth()-1);
			return -tmp;
		}else{
			return mpz_class(0);
		}
	}

	mpz_class maxValue() const
	{
		mpz_class tmp(1);
		tmp=tmp<<(isSigned()?(getWidth()-1):getWidth());;
		return tmp-1;
	}

	virtual std::string getName() const
	{
		std::stringstream tmp;
		tmp<<"Int<signed="<<(isSigned()?"1":"0")<<",width="<<getWidth()<<">";
		return tmp.str();
	}

	static std::shared_ptr<HLSTypeInt> create(bool _signed, int w)
	{ return std::make_shared<HLSTypeInt>(_signed, w); }
};

class HLSTypeBool
	: public HLSTypeFixed
{
public:
	HLSTypeBool()
		: HLSTypeFixed(false, 0, 0)
	{}

	virtual std::string getName() const
	{
		return "bool";
	}

	static std::shared_ptr<HLSTypeBool> create()
	{ return std::make_shared<HLSTypeBool>(); }
};

class HLSTypeFloat
	: public HLSType
{
private:
	int m_wE, m_wF;

public:
	HLSTypeFloat(int wE, int wF)
		: HLSType(3+wE+wF)
		, m_wE(wE)
		, m_wF(wF)
	{}

	virtual bool isFloat() const
	{ return true; }

	int getExponentWidth() const
	{ return m_wE; }

	int getFractionWidth() const
	{ return m_wF; }

	virtual bool equals(const HLSTypePtr &o) const
	{
		auto p=std::dynamic_pointer_cast<HLSTypeFloat>(o);
		if(!p)
			return false;
		return p->getFractionWidth()==getFractionWidth() && p->getExponentWidth()==getExponentWidth();
	}

	virtual std::string getName() const
	{
		std::stringstream tmp;
		tmp<<"Float<wE="<<getExponentWidth()<<",wF="<<getFractionWidth()<<">";
		return tmp.str();
	}

	static std::shared_ptr<HLSTypeFloat> create(int wE, int wF)
	{ return std::make_shared<HLSTypeFloat>(wE, wF); }
};

  HLSTypePtr makeHLSType(const Signal &sig);


}; // flopoco

#endif
