#include "ladder_accumulator.hpp"

#define BOOST_TEST_MODULE AccumulationTests
#include <boost/test/unit_test.hpp>

#include <time.h>

#include "NTL/quad_float.h"
#include "NTL/RR.h"

double getTime()
{
	struct timespec ts;
	clock_gettime(CLOCK_THREAD_CPUTIME_ID,  &ts);
	return ts.tv_sec+1e-9*ts.tv_nsec;
}

double dbl(double x)
{ return x; }

double dbl(flopoco::random::LadderAccumulator x)
{ return x; }

double dbl(const NTL::quad_float &x)
{ return NTL::to_double(x); }

double dbl(const NTL::RR &x)
{ return NTL::to_double(x); }

template<class T>
std::pair<double,double> TimeAcc(T acc)
{
	double TT=0, vv=0.0;
	double n=65536;
	while(TT<2){
		n=n*2;
		acc=0.0;
		volatile double S=getTime();
		for(double i=0;i<n;i++){
			acc += i;
		}
		vv=dbl(acc);
		volatile double E=getTime();
		
		TT=E-S;
		std::cerr<<"n="<<n<<", TT="<<TT<<"\n";
	}
	return std::make_pair(TT/n,vv);
}

BOOST_AUTO_TEST_CASE(TestLadderAcc2)
{
	srand48(1);
	
	NTL::RR::SetPrecision(2200);	// Can do any double-precision calcs exactly
	
	for(int i=0;i<10;i++){
		NTL::RR alt;
		flopoco::random::LadderAccumulator acc;
		
		for(int j=0;j<10;j++){
			double x=ldexp(drand48(), (lrand48()%1000)-500);
			acc += x;
			alt += x;
		}
		
		// Check they are rounding the right way
		BOOST_CHECK_EQUAL(acc.SumDouble(), NTL::to_double(alt));
	}
}

BOOST_AUTO_TEST_CASE(TestLadderAcc)
{
	srand48(1);
	
	NTL::RR::SetPrecision(128);
	
	flopoco::random::LadderAccumulator acc;
	
	acc.Clear();
	acc.Add(nextafter(0.0,1.0));
	BOOST_CHECK_EQUAL(acc.SumDouble(), nextafter(0.0,1.0));
	acc.Add(nextafter(0.0,1.0));
	BOOST_CHECK_EQUAL(acc.SumDouble(), 2*nextafter(0.0,1.0));
	acc.Add(-1);
	BOOST_CHECK_EQUAL(acc.SumDouble(), -1);
	acc.Add(1);
	acc.Add(nextafter(0.0,1.0));
	BOOST_CHECK_EQUAL(acc.SumDouble(), 3*nextafter(0.0,1.0));
	
	acc.Clear();
	acc.Add(DBL_MIN);
	BOOST_CHECK_EQUAL(acc.SumDouble(), DBL_MIN);
	acc.Add(-DBL_MAX);
	BOOST_CHECK_EQUAL(acc.SumDouble(), -DBL_MAX);
	acc.Add(DBL_MIN);
	acc.Add(DBL_MAX);
	BOOST_CHECK_EQUAL(acc.SumDouble(), 2*DBL_MIN);
	
	for(int i=0;i<10;i++){
		acc.Clear();
		
		double tmp=drand48()-drand48();
		acc.Add(tmp);
		BOOST_CHECK_EQUAL(acc.SumDouble() , tmp);
	}
	
	for(int i=0;i<10;i++){
		acc.Clear();
		
		double tmp=drand48()-drand48();
		acc.Add(tmp);
		
		std::vector<double> o;
		for(int j=0;j<1000;j++){
			double tmp2=log(drand48());
			acc.Add(tmp2);
			o.push_back(tmp2);
		}
		while(o.size()>0){
			int ii=lrand48()%o.size();
			acc.Add(-o[ii]);
			o.erase(o.begin()+ii);
		}
		BOOST_CHECK_EQUAL(acc.SumDouble() , tmp);
	}
	
	for(int i=0;i<10;i++){
		acc.Clear();
		for(int j=0;j<1000000;j++){
			acc.Add(sin(j));
		}
		for(int j=0;j<1000000;j++){
			acc.Add(-sin(j));
		}
		BOOST_CHECK(acc.SumDouble()==0);
	}
	
	for(int i=1;i<5;i++){
		double tmp=0;
		acc.Clear();
		for(int j=-100000;j<100000;j++){
			double v=pow(sin(j),i)*pow(0.99999,j);
			acc.Add(v);
			tmp+=v;
		}
		std::cerr<<" double acc = "<<tmp<<", ladder = "<<acc.SumDouble()<<", rel="<<(tmp-acc.SumDouble())/acc.SumDouble()<<"\n";
		for(int j=100000-1;j>=-100000;j--){
			double v=-pow(sin(j),i)*pow(0.99999,j);
			acc.Add(v);
			tmp+=v;
		}
		std::cerr<<" double acc = "<<tmp<<", ladder = "<<acc.SumDouble()<<"\n";
		BOOST_CHECK(acc.SumDouble()==0);
	}
}


BOOST_AUTO_TEST_CASE(TimeLadderAcc)
{
	
	std::pair<double,double> xx_double=TimeAcc((double)0);
	std::cerr<<"double, time="<<xx_double.first<<", acc="<<xx_double.second<<"\n";
	std::pair<double,double> xx_ladder=TimeAcc(flopoco::random::LadderAccumulator());
	std::cerr<<"ladder, time="<<xx_ladder.first<<", acc="<<xx_ladder.second<<"\n";
	std::pair<double,double> xx_quad=TimeAcc(NTL::to_quad_float(0.0));
	std::cerr<<"quad_float, time="<<xx_quad.first<<", acc="<<xx_quad.second<<"\n";
	std::pair<double,double> xx_rr=TimeAcc(NTL::to_RR(0.0));
	std::cerr<<"RR128, time="<<xx_quad.first<<", acc="<<xx_quad.second<<"\n";
	std::cerr<<"  ladder penalty = "<<xx_ladder.first/xx_double.first<<"\n";
	std::cerr<<"  quad penalty = "<<xx_quad.first/xx_double.first<<"\n";
	std::cerr<<"  RR128 penalty = "<<xx_rr.first/xx_double.first<<"\n";
	
}
