
class fixed_type_t
{
private:
	struct info_t
	{
		friend fixed_type_t;
		
		info_t(int _lsb, const mpz_class &_min, const mpz_class &_max)
			: m_lsb(_lsb)
			, m_min(_min)
			, m_max(_max)
		{
			assert(_min <= _max);
			m_isSigned=_min < 0;
			if(m_isSigned){
					assert(0);
			}else{
				assert(0);
			}
		}
		
		info_t(bool isSigned, int width, int lsb)
			: _isSigned(isSigned)
			, _width(width)
			, _lsb(lsb)
		{
			if(isSigned){
				mpz_ui_pow_ui(_min.get(), 2, width-1);
				_max=(-_min)-1;
			}else{
				mpz_ui_pow_ui(_max.get(), 2, width);
				_max=_max-1;
			}
		}
		
		bool m_isSigned;
		int m_width;
		int m_lsb;
		mpz_class m_min, mpz_class m_max;
	};
	typedef boost::shared_ptr_t<info_t> info_ptr_t;
	
	info_ptr_t m_p;
public:
	fixed_type_t(int _lsb, const mpz_class &_min, const mpz_class &_max)
		: m_p(boost::make_shared<info_t>(_lsb, _min, _max))
	{}

	fixed_type_t(bool isSigned, int width, int lsb)
		: m_p(boost::make_shared<info_t>(isSigned, width, lsb))
	{}
	
	bool IsSigned() const
	{ return m_p->m_isSigned; }
	
	int Width() const
	{ return m_p->m_width; }
	
	int LSB() const
	{ return m_p->m_lsb; }
	
	int MSB() const
	{ return LSB()+Width() - (IsSigned()?1:0); }
	
	const mpz_class &MinRaw() const
	{ return m_p->m_min; }
	
	const mpz_class &MaxRaw() const
	{ return m_p->m_max; }
	
	
	void MinReal(mpfr_t res) const
	{ mpfr_set_z_2exp(res, MinRaw().get(), LSB(), MPFR_RNDN); }
	
	void MaxReal(mpfr_t res) const
	{ mpfr_set_z_2exp(res, MaxRaw().get(), LSB(), MPFR_RNDN); }
	
	std::string VHDLType() const
	{
		std::stringstream acc;
		acc<<(IsSigned()?"signed":"unsigned")<<"("<<Width()-1<<" downto 0)";
		return acc.str();
	}
	
	std::string VHDLConstant(const mpz_class &value) const
	{
		std::stringstream acc;
		acc<<"\"";
		for(int i=Width()-1;i>=0;i--){
			acc<<mpz_tstbit(value.get(), i);
		}
		acc<<"\"";
	}
	
	std::string VHDLConstant(mpfr_t x) const
	{
		mpz_class tmp;
		RealToRaw(tmp, x);
		return VHDLConstant(tmp);
	}
};
