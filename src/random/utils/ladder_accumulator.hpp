#ifndef flopoco_random_utils_ladder_accumulator_hpp
#define flopoco_random_utils_ladder_accumulator_hpp

#include <cfloat>
#include <cmath>
#include <cassert>
#include <cstdio>
#include <stdint.h>
#include <vector>
#include <iostream>
#include <stdexcept>

namespace flopoco
{
namespace random
{

class LadderAccumulator
{
	enum{ BIAS= 1023 };
	
private:
	std::vector<int64_t> m_sums;
	
	void Add(unsigned e, int64_t v)
	{
		//fprintf(stderr, "  Add(%d,%lld)\n", e, v);
		
		if(e >= m_sums.size())
			m_sums.resize(e+1, 0);
		
		int64_t acc= (m_sums[e]+=v);
		if( acc >= (1ll<<60)){
			m_sums[e] -= 1ll<<60;
			Add(e+60, +1);
		}else if(acc <= -1ll<<60){
			m_sums[e] += 1ll<<60;
			Add(e+60, -1);
		}
	}
	
	/*! Reduce each level so that the MSB is in {-1,0,+1}, and all other bits
		have the same sign and are zero or one.
		\retval True if the sum is non-zero
	*/
	bool Normalise()
	{
		/*if(m_sums.size()>=1){
			if((m_sums[0]!=0) || (m_sums[1]!=0)){
				throw std::invalid_argument("A denormal has been added to the LAdderAccumulator, which isn't supported at the moment.");
			}
		}*/
		
		int msb_i=-1;
		
		// Propagate everything up to determine the sign of the result
		int i=0;
		while(i<(int)m_sums.size()){
			int64_t v=m_sums[i];
			//fprintf(stderr, "i=%d, v=%lld\n", i, v);
			if(v>1){
				Add(i+1, v/2);
				v=v%2;
			}else if(v<-1){
				Add(i+1, -((-v)/2));
				v=-((-v)%2);
			}
			m_sums[i]=v;
			if(v!=0)
				msb_i=i;
			i++;
		}
		
		if(msb_i==-1)
			return false;
		
		//fprintf(stderr, "  PropogatedSum\n");
		
		// We definitely know whether it is positive, negative, or zero
		// by looking at the MSB, however, there may still be opposite
		// sign bits
		
		// Make sure everything has the same sign
		int64_t sign=m_sums[msb_i];
		//fprintf(stderr, "  MSB=%lld\n", sign);
		
		// P:N -> 0:P
		// PZ:N -> 0P:P
		// PZZ:N -> 0PP:P
		i=m_sums.size()-1;
		while(i>=0){
			if(std::abs(m_sums[i])!=0 && m_sums[i]!=sign){
				assert(m_sums[i]==-sign);
				int o=i+1; // Gauranteed to be at least one bit above
				while(m_sums[o]==0){
					o++;
				}
				assert(m_sums[o]==sign);
				m_sums[o]=0;
				o--;
				while(o>=i){
					m_sums[o]=sign;
					o--;
				}
			}
			i--;
		}
		
		//fprintf(stderr, "End Normalise\n");
		return true;
	}
public:
	LadderAccumulator()
		: m_sums(2050, 0)
	{}

	void Clear()
	{
		std::fill(m_sums.begin(), m_sums.end(), 0);
	}

	/*void Add(double x)
	{
		int e;
		double m=frexp(x, &e);
		
		e += BIAS;
		assert(e>=0);
		
		int64_t im=(int64_t)ldexp(m, 53);
		assert(m==ldexp((double)im, -53));
		
		AddFast(e, im);
	}*/
	
	void Add(double x)
	{
		assert((-DBL_MAX <= x) && (x==x) && (x<=DBL_MAX));
		
		uint64_t vx=*(uint64_t*)&x;
		
		unsigned e=((vx>>52)&0x7FF)+1;
		assert((x==0) || (e>1)); // we don't deal with denormals
		
		int64_t m=(vx&0xFFFFFFFFFFFFFul);
		
		if(e==1){
			if(m==0)
				return;
			e++; // denormal number
		}else{
			m |= (1ull<<52);
		}
		
		if(vx&0x8000000000000000ull){
			int64_t tmp=(m_sums[e] -= m);
			if(tmp<= (-1ll<<60)){
				m_sums[e] += (1ll<<60);
				Add(e+60, -1);
			}
		}else{
			int64_t tmp=(m_sums[e] += m);
			if(tmp>= (1ll<<60)){
				m_sums[e] -= (1ll<<60);
				Add(e+60, +1);
			}
		}
	}
	
	void Add(const LadderAccumulator &a)
	{
		for(unsigned i=0;i<a.m_sums.size();i++){
			Add(i, a.m_sums[i]);
		}
	}
	
	void Sub(const LadderAccumulator &a)
	{
		for(unsigned i=0;i<a.m_sums.size();i++){
			Add(i, -a.m_sums[i]);
		}
	}
	
	double SumDouble()
	{
		//fprintf(stderr, "Begin Sum\n");
		
		if(!Normalise())
			return 0.0;
		
		// We now have the whole thing exactly summed, and normalised
		// into a ladder of values in {-1,0,+1}, so we can fill in the 53 bits
		// from top to bottom
		
		int e=(m_sums.size()-1);
		while((e>=0) ? m_sums[e]==0 : false)
			e=e-1;
		
		if(e<0)
			return 0.0;
			
		double one=ldexp(1, e-BIAS-53); // Represents LSB at current exponent
		
		double acc=0.0;
		int64_t lsb=0;
		for(int i=0;i<53;i++){
			acc += m_sums[e] * one;
			//std::cerr<<"  acc="<<acc<<", s="<<m_sums[e]<<", one="<<one<<", e="<<e<<"\n";
			lsb=m_sums[e];
			e=e-1;
			if(e<0)
				break;
			one=one/2;	// Will this will go denormal?
		}
		
		if(e<0)	// No bits to round, drop out now
			return acc;
		
		// This rounding is almost guaranteed to be implemented wrong :)
		
		if(m_sums[e]==0)
			return acc;	// digit after binary point is zero, rounding towards zero
		e=e-1;
		
		while(e>=0){
			if(m_sums[e]!=0){
				// Ok, so we round away from zero
				return nextafter(acc, acc < 0 ? -DBL_MAX : +DBL_MAX);
			}
			e=e-1;
		}
		
		// We have exactly 0.5. Round to even, although it's not clear what that means for stats
		if(lsb){
			// round away from zero 
			return nextafter(acc, acc < 0 ? -DBL_MAX : +DBL_MAX);
		}else{
			// round towards zero
			return acc;
		}
	}
	
	void operator=(const double &x)
	{
		Clear();
		Add(x);
	}
	
	void operator+=(const double &x)
	{ Add(x); }
	
	void operator-=(const double &x)
	{ Add(-x); }
	
	LadderAccumulator operator+(const LadderAccumulator &o) const
	{
		LadderAccumulator res(*this);
		res.Add(o);
		return res;
	}
	
	LadderAccumulator operator-(const LadderAccumulator &o) const
	{
		LadderAccumulator res(*this);
		res.Sub(o);
		return res;
	}
	
	operator double() // note that this is non-const
	{ return SumDouble(); }
	
	
	
};

}; // random
}; // flopoco

#endif
