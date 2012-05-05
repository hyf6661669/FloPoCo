

class LadderAcc
{
	std::vector<int64_t> m_sums;
	
	void Add(unsigned e, int64_t v)
	{
		if(e >= m_sums.size())
			m_sums.resize(e+1, 0);
		
		int64_t acc=m_sums[e]+v;
		if( acc >= (1ll<<62)){
			acc -= 1ll<<60;
			v=1;
		}else if(acc <= -1ll<<60){
			acc += 1ll<<60;
			v=-1;
		}else{
			v=0;
		}
		m_sums[e]=acc;
		if(v)
			Add(e+60, v);
	}
	
	void Add(double x)
	{
		int e;
		double m=frexp(x, &e);
		
		e += BIAS;
		assert(e>=0);
		
		int64_t im=(int64_t)ldexp(x, 53);
		assert(m==ldexp((double)im, -53));
		
		Add(e, im);
	}
	
	void Normalise()
	{
		int i=0;
		while(i+1<m_sums.size()){
			while( std::abs(2*m_sums[i]) > std::abs(m_sums[i])){
				
			}
		}
	}
};

template<
	class TIt,
	class TF,
	class TRes=TF::result_type
>
TRes Sum(TIt begin, TIt end, const TF &transform)
{
	
}
