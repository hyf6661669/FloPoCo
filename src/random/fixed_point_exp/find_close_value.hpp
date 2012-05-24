

class no_close_value_error : public std::exception
{};	
	
template<class T>
class fixed_point_t
{
private:
	bool isSigned;
	int msb;
	int lsb;
	
public:
	class range_it_t
	{
	private:
		T m_curr, m_delta;
	public:
		range_it_t(const T &curr, const T &delta)
			: m_curr(curr)
			, m_delta(delta)
		{}
		
		bool operator==(const range_it_t &o) const
		{ return m_curr==o.m_curr; }
		
		range_it_t operator ++() const
		{ m_curr+=m_delta; return *this; }
		
		const T &operator*() const
		{ return m_curr; }
	};		
	
	bool IsSigned() const
	{ return 

	unsigned Width() const
	{ return msb-lsb+1; } 
	
	T MinValue() const
	{ return isSigned ? power2<T>(msb) : (T)0.0; }
	
	T MaxValue() const
	{ return power2<T>(msb)-power2<T>(lsb); }
	
	T DeltaValue() const
	{ return power2<T>(lsb); }
	
	typedef range_it_t const_iterator;
	
	const_iterator begin() const
	{ return const_iterator(MinValue(), DeltaValue()); }
	
	const_iterator end() const
	{ return const_iterator(MaxValue()+DeltaValue(), DeltaValue()); }
};

	
template<class T>
struct close_value_params_t
{
	T sigma;
	T mu;
	int fractionBits;
	bool residualSigned;
	int residualMsb; 
	int residualLsb;
};

template<class T>
struct close_value_entry_t
{
private:
	boost::shared_ptr<T> m_params;
	int m_exponent, m_fraction;
	T m_offset, m_error;
public:
	close_value_entry_t(boost::shared_ptr<T> params, T offset)
		: m_params(params)
		, m_offset(offset)
	{}
	
	const T &GetValue(const T &x) const
	{ return exp(m_params->sigma*(m_x-m_offset)+m_params->mu); }
};

//! Find a close appproximation to exp(x*SIGMA+MU)
/*! The return value  is (value,residual), where:
	value*exp(residual*SIGMA+MU) = exp(x*SIGMA+MU),
	residual is a fixed-point number with the specified type,
	and "value" can be represented with a floating-point
	value with "bits" worth of mantissa.
*/
template<class T>
close_value_entry_t<T> FindCloseValue(
	boost::shared_ptr<close_value_params_t> params,
	const T &x
){
	T target=exp(params->sigma*x+params->mu);
	
	bool bestValid=false
	T bestOffset=0, T bestError=0;
	T offset=0, delta=pow(2.0, params->residualLsb), max=pow(2.0, params->residualMsb);
	while(base<max){
		T v=exp(params->sigma*(x-offset)+params->mu);
		T error=(target-v)/target;
		if(params->maxRelError.first <= error && error <=params->maxRelError.second){
			if(!bestValid || (fabs(error) < fabs(bestError))){
				bestOffset=offset;
				bestError=error;
			}
			bestValid=true;
			if(params->stoppingRelError.first <= error && error <= params->stoppingRelError.second){
				break;
			}
		}
		if(params->residualSigned){
			v=exp(params->sigma*(x+offset)+params->mu);
			error=(target-v)/target;
			if(params->maxRelError.first <= error && error <=params->maxRelError.second){
				if(!bestValid || (fabs(error) < fabs(bestError))){
					bestOffset=-offset;
					bestError=error;
				}
				bestValid=true;
				if(params->stoppingRelError.first <= error && error <= params->stoppingRelError.second){
					break;
				}
			}
		}
		offset+=delta;
	}
	if(!bestValid)
		throw no_close_value_error();
	params->worstRelError.first=std::min(
	return close_value_entry_t(params, bestOffset);
}

template<class T>
struct close_value_table_t
{
private:
	boost::shared_ptr<T> m_entryParams;
	
	bool m_inputSigned;
	int m_inputMsb, m_inputLsb;
public:
	close_value_table_t(
		boost::shared_ptr<T> inputParams
	){
		: m_params(params)
		, m_offset(offset)
	{
		
	}
	
	
};
