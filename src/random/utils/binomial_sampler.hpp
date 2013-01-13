#ifndef flopoco_random_utils_binomial_sampler_hpp
#define flopoco_random_utils_binomial_sampler_hpp

namespace flopoco
{
namespace random
{

template<class TCalc,class TRng>
class BinomialSampler
{
private:
	TCalc m_zero, m_one;
	TCalc p, n;
	TCalc r, q;
	TCalc M;
	TCalc p1, p2, p3, p4;
	TCalc x_M, x_L, x_R;
	TCalc c, lambda_L, lambda_R;
	TCalc nrq;

	TRng &m_rng;

	void Setup(TCalc _p, TCalc _n, bool checkInt)
	{
		//std::cerr<<"Setup(_p="<<_p<<", _n="<<_n<<", checkInt="<<checkInt<<")\n";
		
		m_zero=(_p-_p) + (_n-_n);
		m_one=m_zero+1;
		
		_p=_p*m_one;
		_n=_n*m_one;
		
		if(checkInt){
			TCalc tmp=2*_n+1;
			
			if( (tmp==2*_n) || (tmp-1!=2*_n) )
				throw std::string("BinomialSampler - TCalc does not behave as integer for given sample size.");
		}
		
		p=_p;
		n=_n;
		
		if( (p==0) || (p==1)){
			//std::cerr<<"  p==0 or 1\n";
			return;
		}
		
		r=std::min(p,1-p);
		q=1-r;
		
		if(n*r<30){
			//std::cerr<<"  small_mean\n";
			return;
		}
		
		TCalc f_M=n*r+r;
		M=floor(f_M);
		p1=floor(2.195*sqrt(n*r*q)-4.6*q)+0.5;
		x_M=M+0.5;
		x_L=x_M-p1;
		x_R=x_M+p1;
		c=0.134+20.5/(15.3+M);
		TCalc a=(f_M-x_L)/(f_M-x_L*r);
		lambda_L=a*(1+a/2);
		if(x_R*q ==0){
			std::cerr<<"p="<<p<<", n="<<n<<"\n";
			std::cerr<<"r="<<r<<", q="<<q<<"\n";
			std::cerr<<"f_M="<<f_M<<", x_M="<<x_M<<"\n";
			std::cerr<<"p1="<<p1<<"\n";
			std::cerr<<"x_L="<<x_L<<", x_R="<<x_R<<"\n";
		}
		a=(x_R-f_M)/(x_R*q);
		lambda_R=a*(1+a/2);
		p2=p1*(1+2*c);
		p3=p2+c/lambda_L;
		p4=p3+c/lambda_R;
		nrq=n*r*q;
		
		//std::cerr<<"  done setup\n";
	}
	
	TCalc log_f(TCalc y)
	{
		TCalc x1=y+1, f1=M+1, z=n+1-M, w=n-y+1;
		TCalc x2=x1*x1, f2=f1*f1, z2=z*z, w2=w*w;
		TCalc acc = x_M *log(f1/x1);
		acc += (n-M+0.5)*log(z/w);
		acc+=(y-M)*log(w*r/(x1*q));
		acc+=(13860.-(462.-(132.-(99.-140./f2)/f2)/f2)/f2)/f1/166320;
		acc+=(13860.-(462.-(132.-(99.-140./z2)/z2)/z2)/z2)/z/166320;
		acc+=(13860.-(462.-(132.-(99.-140./x2)/x2)/x2)/x2)/x1/166320;
		acc+=(13860.-(462.-(132.-(99.-140./w2)/w2)/w2)/w2)/w/166320;
		return acc;
	}
	
	/* The stirling approximation should get more accurate for large n,
		so for our purposes (huge samples) it should be fine. */
	TCalc log_f_opt(TCalc y)
	{
		TCalc x1=y+1, f1=M+1, z=n+1-M, w=n-y+1;
		assert(x1!=0);
		assert(f1!=0);
		assert(z!=0);
		assert(w!=0);
		TCalc ix1=1.0/x1, if1=1.0/f1, iz=1.0/z, iw=1.0/w;
		
		TCalc ix2=ix1*ix1, if2=if1*if1, iz2=iz*iz, iw2=iw*iw;
		
		assert(x1*q!=0);
		
		TCalc acc = x_M *log(f1*ix1);
		acc += (n-M+0.5)*log(z*iw);
		acc+=(y-M)*log(w*r/(x1*q));
		acc+=(13860.-(462.-(132.-(99.-140.*if2)*if2)*if2)*if2)*if1*(1.0/166320);
		acc+=(13860.-(462.-(132.-(99.-140.*iz2)*iz2)*iz2)*iz2)*iz*(1.0/166320);
		acc+=(13860.-(462.-(132.-(99.-140.*ix2)*ix2)*ix2)*ix2)*ix1*(1.0/166320);
		acc+=(13860.-(462.-(132.-(99.-140.*iw2)*iw2)*iw2)*iw2)*iw*(1.0/166320);
		return acc;
	}
	
	TCalc log_f_a(TCalc y)
	{
		TCalc acc=0.0;
		acc += gamma_log(M+1);
		acc += gamma_log(n-M+1);
		acc += -gamma_log(y+1);
		acc += -gamma_log(n-y+1);
		acc += log(p/(1-p)) * (y-M);
		return acc.Sum();
	}
	
	TCalc small_mean()
	{
		//std::cerr<<"small_mean\n";
		
		TCalc qn=pow(q,n);
		assert(q!=0);
		TCalc rr=r/q;
		TCalc g=rr*(n+1);
		
		const int MAX_ITER=10000;
		int iterations=0;
		
		while(1){	// 14
			TCalc ix;
			ix=0;
			TCalc f=qn;
			TCalc u=m_rng();
			
			while(1){ // 15
				if(u < f){
					//std::cerr<<"  done\n";
					return (p>0.5) ? n-ix : ix;
				}
				
				iterations++;
				if(iterations>MAX_ITER)
					throw std::string("BinomialSample - small_mean seems to have locked up.");
					
				if(ix > 110){
					break;	// goto 14
				}
				u=u-f;
				ix=ix+1;
				assert(ix!=0);
				f=f*(g/ix-rr);
				// goto 15
			}
		}
	}
public:
	BinomialSampler(TCalc p, TCalc n, TRng &r, bool checkInt=true)
		: m_rng(r)
	{
		Setup(p,n, checkInt);
	}
	
	void Reset(TCalc p, TCalc n, bool checkInt=true)
	{
		Setup(p,n,checkInt);
	}

	TCalc operator()()
	{
		//std::cerr<<"operator(), x_L="<<x_L<<"\n";
		
		TCalc y;
		
		if(p==1){
			return n;
		}
		if(p==0){
			return m_zero;
		}
		
		//volatile double pp=to<double>(p);
		//pp=pp;
		//volatile double nn=to<double>(n);
		//nn=nn;
		
		if(n*r<30)
			return small_mean();
		
		//std::cerr<<"  x_M="<<x_M<<"\n";
		
		const int MAX_ITER=10000;
		int iterations=0;
		
		while(1){
			if(iterations++ > MAX_ITER)
				throw std::string("BinomialSampler - Loop seems to have locked up. Increase precision?");
			
			TCalc u=m_rng()*p4;
			TCalc v=m_rng();
			//TCalc u=0.5*p4, v=0.5;
			//std::cerr<<" u="<<u<<", v="<<v<<",  p4="<<p4<<"\n";
			
			if(u<=p1){
				y=floor(x_M-p1*v+u);
				// goto 6
				//std::cerr<<"  Tri : "<<y<<"\n";
				break;
			}else{
				if(u<=p2){
					// Parallelograms
					//std::cerr<<"Para, x_L="<<x_L<<", u="<<u<<", p1="<<p1<<", c="<<c<<"\n";
					assert(c!=0);
					TCalc x=x_L+(u-p1)/c;
					//std::cerr<<"  x="<<x<<"\n";
					assert(p1!=0);
					v=v*c+1-abs(M-x+0.5)/p1;
					//std::cerr<<"  cond\n";
					if((v>1) || (v==0))
						continue;
					y=floor(x);
					//std::cerr<<"  Parllelogram : "<<y<<"\n";
				}else if(u<=p3){
					//std::cerr<<"LeftTail\n";
					// left tail
					assert(lambda_L!=0);
					y=floor(x_L+log(v)/lambda_L);
					if(y<0)
						continue;
					v=v*(u-p2)*lambda_L;
					//std::cerr<<"  Left : "<<y<<"\n";
				}else{
					//std::cerr<<"RightTail\n";
					// right tail
					assert(lambda_R!=0);
					y=floor(x_R-log(v)/lambda_R);
					if(y>n)
						continue;
					v=v*(u-p3)*lambda_R;
					//std::cerr<<"  Right : "<<y<<"\n";
				}
				
				TCalc k=abs(y-M);
				if((k<20) || (k>(n*r*q)/2-1)){
					// Recursive f(y)
					TCalc s=r/q;
					TCalc a=s*(n+1);
					TCalc F;
					F=1.0;
					if(M<y){
						for(TCalc i=M+1;i<=y;i++){
							assert(i!=0);
							F=F*(a/i-s);
						}
					}else{
						for(TCalc i=y+1;i<=M;i++){
							assert(i!=0);
							assert(a/i-s!=0);
							F=F/(a/i-s);
						}
					}
					if(v<=F)
						break;	// SUccess
					
					// Else try again
				}else{
					// Squeeze
					assert(nrq!=0);
					TCalc rho=(k/nrq)*((k*(k/3+0.625)+1.0/6)/nrq+0.5);
					TCalc t=-k*k/(2*nrq);
					TCalc A=log(v);
					if(A<t-rho)
						break;	// Success;
					if(A>t+rho)
						continue;	// Start again
					
					// Accurate
					if(A <= log_f_opt(y))
						break;	// success
					
					// try again
				}
			}
		}
		// Step 6
		if(p>0.5){
			y=n-y;
		}
		return y;
	}
	
	static TCalc Sample(TCalc p, TCalc n, TRng &r, bool checkInt=false)
	{
		BinomialSampler<TCalc,TRng> tmp(p,n,r,checkInt);
		return tmp();
	}
};

template<class TCalc,class TRng>
TCalc GenerateBinomialSample(TCalc p, TCalc n, TRng &r, bool checkInt=false)
{
	return BinomialSampler<TCalc,TRng>::Sample(p,n,r,checkInt);
}

}; // random
}; // flopoco

#endif
