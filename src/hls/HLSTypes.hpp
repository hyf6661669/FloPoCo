#ifndef flopoco_hls_types_hpp
#define flopoco_hls_types_hpp

#include <memory>

namespace flopoco
{

class HLSType;
typedef std::shared_ptr<HLSType> HLSTypePtr;


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
};


class HLSTypeInt
	: public HLSType
{
private:
	bool m_isSigned;

public:
	HLSTypeInt(bool _signed, int w)
		: HLSType(w)
		, m_isSigned(_signed)
	{}

	bool isSigned() const
	{ return m_isSigned; }

	virtual bool equals(const HLSTypePtr &o) const
	{
		auto p=std::dynamic_pointer_cast<HLSTypeInt>(o);
		if(!p)
			return false;
		return p->isSigned()==isSigned() && p->getWidth()==getWidth();
	}

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

	static std::shared_ptr<HLSTypeInt> create(bool _signed, int w)
	{ return std::make_shared<HLSTypeInt>(_signed, w); }
};


class HLSTypeBool
	: public HLSType
{
public:
	HLSTypeBool()
		: HLSType(1)
	{}

	virtual bool equals(const HLSTypePtr &o) const
	{
		auto p=std::dynamic_pointer_cast<HLSTypeBool>(o);
		if(!p)
			return false;
		return true;
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

	static std::shared_ptr<HLSTypeFloat> create(int wE, int wF)
	{ return std::make_shared<HLSTypeFloat>(wE, wF); }
};


}; // flopoco

#endif
