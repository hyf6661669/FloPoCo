#ifndef flopoco_random_utils_fft_big_hpp
#define flopoco_random_utils_fft_big_hpp

namespace flopoco
{
namespace random
{
	namespace detail{
		
		inline bool IsBinaryPower(unsigned n)
		{
			if(n==0)
				return false;
			while(n>1)
				n=n>>1;
			return n==1;
		}
		
		inline unsigned NextBinaryPower(unsigned n)
		{
			unsigned res=1;
			while(res<n){
				unsigned next=res<<1;
				if(next<res)
					throw std::string("NextBinaryPower - Overflow.");	// Happens for all ones...
				res=next;
			}
			return res;
		}
	}
	
	// General purpose in-place complex->complex fft
	template<class T>
	void fft(T *data, unsigned n);

	// General purpose in-place complex->complex ifft
	template<class T>
	void ifft(T *data, unsigned n);
};
};

#endif
