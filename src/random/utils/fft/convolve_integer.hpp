

class FieldProvider
{
	typedef ... element_t;
	
	template<class T>
	void import_element(element_t &dst, const T & src);
	
	template<class T>
	void export_element(T &dst, const element_t & src);
	
	void mul(element_t &dst, const element_t &a, const element_t &b);
	
	void add(element_t &dst, const element_t &a, const element_t &b);
};

template<class T>
class GenericPrimeField
{
private:	
	T m_p;
public:
	GenericPrimeField(T p)
		: m_p(p)
	{}

	typedef T element_t;
	typedef struct{
		T real;
		T imag;
	}complex_t;
	
	template<class TSrc>
	void import_element(element_t &dst, const TSrc & src)
	{ dst=src; }
	
	template<class TDst>
	void export_element(TDst &dst, const element_t & src)
	{ dst=src; }
	
	void mul(element_t &dst, const element_t &a, const element_t &b)
	{
		assert((0<=a) && (a<m_p));
		assert((0<=b) && (b<m_p));
		
		dst=(a*b)%m_p;
	}
	
	void add(element_t &dst, const element_t &a, const element_t &b)
	{
		assert((0<=a) && (a<m_p));
		assert((0<=b) && (b<m_p));
		
		dst=a+b;
		if(dst >= m_p)
			dst -= m_p;
	}
	
	void add(complex_t &dst, const complex_t &a, const complex_t &b)
	{
		dst=a+b;
	}
	
	typedef complex_t root_t;
	
	void root_of_unity(
};

