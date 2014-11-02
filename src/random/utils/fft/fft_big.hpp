#ifndef flopoco_random_utils_fft_big_hpp
#define flopoco_random_utils_fft_big_hpp

namespace flopoco
{
namespace random
{
	namespace detail{
		
		bool IsBinaryPower(unsigned n)
		{
			if(n==0)
				return false;
			while(0==(n%2))
				n=n/2;
			return n==1;
		}
		
		unsigned NextBinaryPower(unsigned n)
		{
			unsigned curr=1;
			while(n>curr){
				unsigned prev=curr;
				curr<<=1;
				if(curr<prev)
					throw std::string("NextBinaryPower - wrapped!");
			}
			return curr;
		}
		
		template<class T>
		T Pi();
		
		template<class T>
		T twice(const T &x);
		
		template<class T>
		T square(const T &x);
		
		// FFT from numerical recipes. Input is as a zero-based array (i.e. normal C, not NumRec style)
		// TODO : What is the license on numerical recipes...
		template<class T>
		void four1(T *data, unsigned nn, bool inverse)
		{
			if(!IsBinaryPower(nn))
				throw std::string("four1 - FFT Input data size must be power of 2.");
			
			const T pi(Pi<T>());
			
			data=data-1;	// Move to stupid one-based index
			
			unsigned n, mmax, m, j, istep, i;
			T wtemp, wr, wpr, wpi, wi, theta;
			T tempr, tempi;
			
			n=nn<<1;
			j=1;
			for(i=1;i<n;i+=2){
				if(j>i){
					std::swap(data[j], data[i]);
					std::swap(data[j+1], data[i+1]);
				}
				m=n>>1;
				while(m>=2 && j>m){
					j-=m;
					m>>=1;
				}
				j+=m;
			}
			
			mmax=2;
			while(n>mmax){
				istep=mmax<<1;
				theta=twice(pi)/mmax;
				if(inverse)
					theta=-theta;
				wtemp=sin(half(theta));
				wpr= -twice(square(wtemp));
				wpi=sin(theta);
				wr=1.0;
				wi=0.0;
				for(m=1;m<mmax;m+=2){
					for(i=m;i<=n;i+=istep){
						j=i+mmax;
						tempr=wr*data[j]-wi*data[j+1];
						tempi=wr*data[j+1]+wi*data[j];
						data[j]=data[i]-tempr;
						data[j+1]=data[i+1]-tempi;
						data[i]+=tempr;
						data[i+1]+=tempi;
					}
					wtemp=wr;
					wr=wr*wpr-wi*wpi+wr;
					wi=wi*wpr+wtemp*wpi+wi;
				}
				mmax=istep;
			}
			
			if(inverse){
				T scale;
				scale=1.0;
				scale=scale/nn;
				for(unsigned i=1;i<=2*nn;i++){
					data[i]*=scale;
				}
			}
		}
}; // detail
};
};

#endif
