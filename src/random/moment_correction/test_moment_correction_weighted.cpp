#include "correct_distribution.hpp"

#include "random/distributions/gaussian_distribution.hpp"
#include "random/distributions/table_distribution.hpp"

#define BOOST_TEST_MODULE MomentCorrectionTestsWeighted
#include <boost/test/unit_test.hpp>

using namespace flopoco::random;

std::vector<std::pair<double,double> > MakeTriApprox(ContinuousDistribution<double>::TypePtr target, int n)
{
	std::vector<double> tri(2*n, 0.0);
	for(unsigned i=0;i<n;i++){
		double p=(i+1)/(n*(double)n);
		tri[2*i]=p;
		tri[2*i+1]=p;
	}
	tri.resize(2*n-1);
	
	std::vector<std::pair<double,double> > x(2*tri.size());
	
	double cdf=0.0;
	for(int i=0;i<tri.size();i++){
		double p=tri[i];
		cdf+=p;
		//fprintf(stderr, "%lg, %lg\n", p, cdf);
		x[i].first=target->InvCdf(cdf/2);
		x[i].second=p/2;
		x[x.size()-i-1].first=x[i].first==0 ? 0 : -x[i].first;
		x[x.size()-i-1].second=p/2;
	}
	
	///for(unsigned i=0;i<x.size();i++){
	//	fprintf(stderr, "%d, %lg, %lg\n", i, x[i].first, x[i].second);
	//}
	
	return x;
}

BOOST_AUTO_TEST_CASE(SymmetricCubicTest)
{
	ContinuousDistribution<double>::TypePtr target=boost::make_shared<GaussianDistribution<double> >();
	
	int n=8;
	std::vector<std::pair<double,double> > tab=MakeTriApprox(target, n);
	
	TableDistribution<double>::TypePtr current=boost::make_shared<TableDistribution<double> >(tab);
	
	std::vector<double> poly=FindSymmetricCubicPolynomialCorrection<double>(
		current,
		target
	);
	if(poly.size()==0){
		std::cerr<<"No solution\n";
	}else{
		std::cerr<<"Poly = "<<poly[0];
		for(unsigned i=1;i<poly.size();i++){
			std::cerr<<" + "<<poly[i]<<"*x^"<<i;
		}
		std::cerr<<"\n";
	}
	
	TableDistribution<double>::TypePtr corrected=current->ApplyPolynomial(poly);
	
	for(int i=0;i<=6;i++){
		fprintf(stderr, "%12.10lg, %12.10lg, %12.10lg\n", target->RawMoment(i), current->RawMoment(i), corrected->RawMoment(i));
	}
}

BOOST_AUTO_TEST_CASE(SymmetricQuinticTest)
{
	ContinuousDistribution<double>::TypePtr target=boost::make_shared<GaussianDistribution<double> >();
	
	int n=32;
	std::vector<std::pair<double,double> > tab=MakeTriApprox(target, n);
	
	TableDistribution<double>::TypePtr current=boost::make_shared<TableDistribution<double> >(tab);
	
	std::vector<double> poly=FindSymmetricQuinticPolynomialCorrection<double>(
		current,
		target
	);
	if(poly.size()==0){
		std::cerr<<"No solution\n";
	}else{
		std::cerr<<"Poly = "<<poly[0];
		for(unsigned i=1;i<poly.size();i++){
			std::cerr<<" + "<<poly[i]<<"*x^"<<i;
		}
		std::cerr<<"\n";
	}
	
	TableDistribution<double>::TypePtr corrected=current->ApplyPolynomial(poly);
	
	for(int i=0;i<=10;i++){
		fprintf(stderr, "%12.10lg, %12.10lg, %12.10lg\n", target->RawMoment(i), current->RawMoment(i), corrected->RawMoment(i));
	}
}

BOOST_AUTO_TEST_CASE(SymmetricHepticTest)
{
	ContinuousDistribution<double>::TypePtr target=boost::make_shared<GaussianDistribution<double> >();
	
	int n=128;
	std::vector<std::pair<double,double> > tab=MakeTriApprox(target, n);
	
	TableDistribution<double>::TypePtr current=boost::make_shared<TableDistribution<double> >(tab);
	
	std::vector<double> poly=FindSymmetricHepticPolynomialCorrection<double>(
		current,
		target
	);
	if(poly.size()==0){
		std::cerr<<"No solution\n";
	}else{
		std::cerr<<"Poly = "<<poly[0];
		for(unsigned i=1;i<poly.size();i++){
			std::cerr<<" + "<<poly[i]<<"*x^"<<i;
		}
		std::cerr<<"\n";
	}
	
	TableDistribution<double>::TypePtr corrected=current->ApplyPolynomial(poly);
	
	for(int i=0;i<=14;i++){
		fprintf(stderr, "%12.10lg, %12.10lg, %12.10lg\n", target->RawMoment(i), current->RawMoment(i), corrected->RawMoment(i));
	}
}
