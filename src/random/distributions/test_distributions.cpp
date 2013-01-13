#include "distribution.hpp"
#include "clt_distribution.hpp"
#include "table_distribution.hpp"
#include "histogram_distribution.hpp"
#include "gaussian_distribution.hpp"
#include "discrete_gaussian_distribution.hpp"
#include "quantised_gaussian_distribution.hpp"

#include "sum_distribution.hpp"

#include "random/distributions/chi2_partition.hpp"

#include "moment_conversions.hpp"

#define BOOST_TEST_MODULE TestDistributions
#include <boost/test/unit_test.hpp>

using namespace flopoco::random;


BOOST_AUTO_TEST_CASE(HistogramDistributionTests)
{
	typedef HistogramDistribution<double>::TypePtr HistogramDistributionPtr;
	
	std::vector<double> data;
	data.push_back(0.25);
	data.push_back(0.50);
	data.push_back(0.25);
	
	HistogramDistributionPtr ptr=boost::make_shared<HistogramDistribution<double> >(-1, 1, data);
	
	BOOST_CHECK(ptr->ElementCount()==3);
	
	BOOST_CHECK_EQUAL(ptr->Pmf(-1.5), 0);
	BOOST_CHECK_CLOSE(ptr->Pmf(-1), 0.25, 1e-10);
	BOOST_CHECK_EQUAL(ptr->Pmf(-.5), 0);
	BOOST_CHECK_CLOSE(ptr->Pmf(0), 0.5, 1e-10);
	BOOST_CHECK_EQUAL(ptr->Pmf(.5), 0);
	BOOST_CHECK_CLOSE(ptr->Pmf(+1), 0.25, 1e-10);
	BOOST_CHECK_EQUAL(ptr->Pmf(1.5), 0);
	
	BOOST_CHECK_EQUAL(ptr->Cdf(-1.5), 0);
	BOOST_CHECK_EQUAL(ptr->Cdf(-1), 0.25);
	BOOST_CHECK_EQUAL(ptr->Cdf(-.5), 0.25);
	BOOST_CHECK_EQUAL(ptr->Cdf(0), 0.75);
	BOOST_CHECK_EQUAL(ptr->Cdf(.5), 0.75);
	BOOST_CHECK_EQUAL(ptr->Cdf(+1), 1);
	BOOST_CHECK_EQUAL(ptr->Cdf(1.5), 1.0);
	
	typename DiscreteDistribution<double>::TypePtr tri=SelfAddHistogramDistributions<double>(ptr,2);
	
	BOOST_CHECK_CLOSE(tri->Pmf(-2.5), 0, 1e-10);
	BOOST_CHECK_CLOSE(tri->Pmf(-2), 0.0625, 1e-10);
	BOOST_CHECK_CLOSE(tri->Pmf(-1.5), 0, 1e-10);
	BOOST_CHECK_CLOSE(tri->Pmf(-1), 0.25, 1e-10);
	BOOST_CHECK_CLOSE(tri->Pmf(-0.5), 0, 1e-10);
	BOOST_CHECK_CLOSE(tri->Pmf(0), 0.375, 1e-10);
	BOOST_CHECK_CLOSE(tri->Pmf(.5), 0, 1e-10);
	BOOST_CHECK_CLOSE(tri->Pmf(1), 0.25, 1e-10);
	BOOST_CHECK_CLOSE(tri->Pmf(1.5), 0, 1e-10);
	BOOST_CHECK_CLOSE(tri->Pmf(2), 0.0625, 1e-10);
	BOOST_CHECK_CLOSE(tri->Pmf(2.5), 0, 1e-10);
}


BOOST_AUTO_TEST_CASE(CLTDistributionTests)
{
	typedef typename TableDistribution<double>::TypePtr TableDistributionPtr;
	
	TableDistributionPtr clt=flopoco::random::MakeCLTDistribution<double>(3, 2);
	
	std::vector<std::pair<double,double> > elements=clt->GetElements();
	
	for(unsigned i=0;i<elements.size();i++){
		std::pair<double,double> v=elements[i];
		//fprintf(stderr, "%lg, %lg\n", v.first, v.second);
	}
	
	//fprintf(stderr, "IsSym=%d\n", clt->IsSymmetric()?1:0);
}


double ReferenceRawMoment(unsigned n, const double *v, unsigned k)
{
	double acc=0;
	for(unsigned i=0;i<n;i++){
		acc += pow(v[i],(int)k);
	}
	return acc / n;
}

double ReferenceCentralMoment(unsigned n, const double *v, unsigned k)
{
	double mean=std::accumulate(v, v+n, 0.0)/n;
	double acc=0;
	for(unsigned i=0;i<n;i++){
		acc += pow(v[i]-mean,(int)k);
	}
	return acc / n;
}

BOOST_AUTO_TEST_CASE(MomentConversionTests)
{
	// Standardised normal moments, come through unscathed
	double normal[]={1,0,1,0,3,0,15,0,105};
	for(unsigned i=0;i<=8;i++){
		BOOST_CHECK(RawMomentsToCentralMoment(i, normal)==normal[i]);
	}
	
	int K=8;
	for(int i=0;i<10;i++){
		int n=1+(lrand48()%10000);
		std::vector<double> values(n);
		for(int j=0;j<n;j++){
			values[j]=drand48()*10;
		}
		double mean=std::accumulate(values.begin(), values.end(), 0.0)/n;
		
		std::vector<double> raw(K+1, 0), ref(K+1,0);
		for(int j=0;j<=K;j++){
			raw[j]=ReferenceRawMoment(n, &values[0], j);
			ref[j]=ReferenceCentralMoment(n, &values[0], j);
		}
		
		BOOST_CHECK_CLOSE(raw[0], 1, 1e-10);
		BOOST_CHECK_CLOSE(raw[1], mean, 1e-10);
		for(int j=2;j<=K;j++){
			BOOST_CHECK(abs(RawMomentsToCentralMoment(j, &raw[0])-ref[j]) < 1e-10);
			BOOST_CHECK(abs(CentralMomentsToRawMoment(j, raw[1], &ref[0])-raw[j]) < 1e-10);
		}
		
		TableDistribution<double> tab(&values[0], &values[n]);
		
		for(int j=1;j<=K;j++){
			BOOST_CHECK(abs(raw[j]-tab.RawMoment(j)) < 1e-10);
			BOOST_CHECK(abs(ref[j]-  tab.CentralMoment(j))< 1e-10);
		}
	}
}

BOOST_AUTO_TEST_CASE(TableDistributionTests)
{
	typedef TableDistribution<double>::TypePtr TableDistributionPtr;
	
	double data[4]={-1 , +1, 0};
	
	TableDistributionPtr ptr=boost::make_shared<TableDistribution<double> >(data, data+3);
	
	BOOST_CHECK(ptr->ElementCount()==3);
	
	BOOST_CHECK_CLOSE(ptr->Pmf(-1), 1.0/3, 1e-10);
	BOOST_CHECK(ptr->Pmf(-0.5)==0.0);
	BOOST_CHECK_CLOSE(ptr->Pmf(0), 1.0/3, 1e-10);
	BOOST_CHECK(ptr->Pmf(0.5)==0.0);
	BOOST_CHECK_CLOSE(ptr->Pmf(+1), 1.0/3, 1e-10);
	
	BOOST_CHECK(ptr->Cdf(-1.1)==0.0);
	BOOST_CHECK_CLOSE(ptr->Cdf(-1), 1.0/3, 1e-10);
	BOOST_CHECK_CLOSE(ptr->Cdf(-0.5), 1.0/3, 1e-10);
	BOOST_CHECK_CLOSE(ptr->Cdf(0), 2.0/3, 1e-10);
	BOOST_CHECK_CLOSE(ptr->Cdf(0.5), 2.0/3.0, 1e-10);
	BOOST_CHECK(ptr->Cdf(+1)==1.0);
	BOOST_CHECK(ptr->Cdf(1.1)==1.0);
	
	BOOST_CHECK(ptr->StandardMoment(0)==1.0);
	BOOST_CHECK(ptr->StandardMoment(1)==0.0);
	BOOST_CHECK_CLOSE(ptr->StandardMoment(2), sqrt(2.0/3.0), 1e-10);
	BOOST_CHECK(ptr->StandardMoment(3)==0.0);
	BOOST_CHECK_CLOSE(ptr->StandardMoment(4), 1.5, 1e-10);
	BOOST_CHECK(ptr->StandardMoment(5)==0.0);
	BOOST_CHECK_CLOSE(ptr->StandardMoment(6), 2.25, 1e-10);
}


BOOST_AUTO_TEST_CASE(GaussianDistributionTests)
{
	typedef ContinuousDistribution<double>::TypePtr ContinuousDistributionPtr;
	
	ContinuousDistributionPtr ptr=boost::make_shared<GaussianDistribution<double> >();
	
	BOOST_CHECK(ptr->StandardMoment(0)==1.0);
	BOOST_CHECK(ptr->StandardMoment(1)==0.0);
	BOOST_CHECK(ptr->StandardMoment(2)==1.0);
	BOOST_CHECK(ptr->StandardMoment(3)==0.0);
	BOOST_CHECK(ptr->StandardMoment(4)==3.0);
	BOOST_CHECK(ptr->StandardMoment(5)==0.0);
	BOOST_CHECK(ptr->StandardMoment(6)==15.0);
	BOOST_CHECK(ptr->StandardMoment(7)==0.0);
	BOOST_CHECK(ptr->StandardMoment(8)==105.0);
	
	BOOST_CHECK(ptr->RawMoment(0)==1.0);
	BOOST_CHECK(ptr->RawMoment(1)==0.0);
	BOOST_CHECK(ptr->RawMoment(2)==1.0);
	BOOST_CHECK(ptr->RawMoment(3)==0.0);
	BOOST_CHECK(ptr->RawMoment(4)==3.0);
	BOOST_CHECK(ptr->RawMoment(5)==0.0);
	BOOST_CHECK(ptr->RawMoment(6)==15.0);
	BOOST_CHECK(ptr->RawMoment(7)==0.0);
	BOOST_CHECK(ptr->RawMoment(8)==105.0);
	
	BOOST_CHECK(ptr->CentralMoment(0)==1.0);
	BOOST_CHECK(ptr->CentralMoment(1)==0.0);
	BOOST_CHECK(ptr->CentralMoment(2)==1.0);
	BOOST_CHECK(ptr->CentralMoment(3)==0.0);
	BOOST_CHECK(ptr->CentralMoment(4)==3.0);
	BOOST_CHECK(ptr->CentralMoment(5)==0.0);
	BOOST_CHECK(ptr->CentralMoment(6)==15.0);
	BOOST_CHECK(ptr->CentralMoment(7)==0.0);
	BOOST_CHECK(ptr->CentralMoment(8)==105.0);
	
	BOOST_CHECK_CLOSE(ptr->Pdf(0.0), 0.398942280401433, 1e-10);
	BOOST_CHECK_CLOSE(ptr->Pdf(1.0), 0.241970724519143, 1e-10);
	BOOST_CHECK_CLOSE(ptr->Pdf(-1.0), 0.241970724519143, 1e-10);
	
	BOOST_CHECK_CLOSE(ptr->Cdf(0.0), 0.5, 1e-10);
	BOOST_CHECK_CLOSE(ptr->Cdf(1.0), 0.841344746068543, 1e-10);
	BOOST_CHECK_CLOSE(ptr->Cdf(-1.0),  0.1586552539314573, 1e-10);
	
	BOOST_CHECK_CLOSE(ptr->InvCdf(0.5), 0.0, 1e-10);
	BOOST_CHECK_CLOSE(ptr->InvCdf(0.841344746068543), 1.0, 1e-10);
	BOOST_CHECK_CLOSE(ptr->InvCdf(0.1586552539314573), -1.0 , 1e-10);	
	
	
	ptr=boost::make_shared<GaussianDistribution<double> >(1.0,2.0);
	
	BOOST_CHECK(ptr->StandardMoment(0)==1.0);
	BOOST_CHECK(ptr->StandardMoment(1)==1.0);
	BOOST_CHECK(ptr->StandardMoment(2)==2.0);
	BOOST_CHECK(ptr->StandardMoment(3)==0.0);
	BOOST_CHECK(ptr->StandardMoment(4)==3.0);
	BOOST_CHECK(ptr->StandardMoment(5)==0.0);
	BOOST_CHECK(ptr->StandardMoment(6)==15.0);
	BOOST_CHECK(ptr->StandardMoment(7)==0.0);
	BOOST_CHECK(ptr->StandardMoment(8)==105.0);
	
	BOOST_CHECK(ptr->CentralMoment(0)==1.0);
	BOOST_CHECK(ptr->CentralMoment(1)==0.0);
	BOOST_CHECK(ptr->CentralMoment(2)==1.0*4);
	BOOST_CHECK(ptr->CentralMoment(3)==0.0);
	BOOST_CHECK(ptr->CentralMoment(4)==3.0*16);
	BOOST_CHECK(ptr->CentralMoment(5)==0.0);
	BOOST_CHECK(ptr->CentralMoment(6)==15.0*64);
	BOOST_CHECK(ptr->CentralMoment(7)==0.0);
	BOOST_CHECK(ptr->CentralMoment(8)==105.0*256);
	
	// These came from integration in sage
	BOOST_CHECK(ptr->RawMoment(0)== 1);
	BOOST_CHECK(ptr->RawMoment(1)== 1);
	BOOST_CHECK(ptr->RawMoment(2)== 5);
	BOOST_CHECK(ptr->RawMoment(3)== 13);
	BOOST_CHECK(ptr->RawMoment(4)== 73);
	BOOST_CHECK(ptr->RawMoment(5)==281);
	BOOST_CHECK(ptr->RawMoment(6)== 1741);
	BOOST_CHECK(ptr->RawMoment(7)== 8485);
	BOOST_CHECK(ptr->RawMoment(8)== 57233);

	
	BOOST_CHECK_CLOSE(ptr->Pdf(0.0*2.0+1.0), 0.398942280401433/2.0, 1e-10);
	BOOST_CHECK_CLOSE(ptr->Pdf(1.0*2.0+1.0), 0.241970724519143/2.0, 1e-10);
	BOOST_CHECK_CLOSE(ptr->Pdf(-1.0*2.0+1.0), 0.241970724519143/2.0, 1e-10);
	
	BOOST_CHECK_CLOSE(ptr->Cdf(0.0*2.0+1.0), 0.5, 1e-10);
	BOOST_CHECK_CLOSE(ptr->Cdf(1.0*2.0+1.0), 0.841344746068543, 1e-10);
	BOOST_CHECK_CLOSE(ptr->Cdf(-1.0*2.0+1.0),  0.1586552539314573, 1e-10);
	
	BOOST_CHECK_CLOSE(ptr->InvCdf(0.5), 0.0*2.0+1.0, 1e-10);
	BOOST_CHECK_CLOSE(ptr->InvCdf(0.841344746068543), 1.0*2.0+1.0, 1e-10);
	BOOST_CHECK_CLOSE(ptr->InvCdf(0.1586552539314573), -1.0*2.0+1.0 , 1e-10);	
}

BOOST_AUTO_TEST_CASE(DiscreteGaussianDistributionTests)
{
	typedef DiscreteGaussianDistribution<double>::TypePtr DiscreteGaussianDistributionPtr;
	
	for(double sigma=0.2;sigma<=5.0;sigma*=1.5){
		for(int fb=-1;fb<=3;fb++){
			if(ldexp(sigma,fb)<2)
				continue;
			
			DiscreteGaussianDistributionPtr ptr=boost::make_shared<DiscreteGaussianDistribution<double> >(sigma, fb);
	
			BOOST_CHECK_CLOSE(ptr->StandardMoment(0),1.0, 1e-8);
			BOOST_CHECK_CLOSE(ptr->StandardMoment(2),sigma, 1e-8);
			BOOST_CHECK_CLOSE(ptr->StandardMoment(4),3.0, 1e-8);
			BOOST_CHECK_CLOSE(ptr->StandardMoment(6),15.0, 1e-8);
			BOOST_CHECK_CLOSE(ptr->StandardMoment(8),105.0,1e-8);
			BOOST_CHECK_EQUAL(ptr->StandardMoment(1),0.0);
			BOOST_CHECK_EQUAL(ptr->StandardMoment(3),0.0);
			BOOST_CHECK_EQUAL(ptr->StandardMoment(5),0.0);
			
			double x_scale=1/(sigma*sigma*2);
			double p_scale=1/(pow(2.0,fb)*sigma*sqrt(2*3.1415926535897932384626433832795));
			BOOST_CHECK_CLOSE(ptr->Pmf(0), p_scale, 1e-10);
			BOOST_CHECK_CLOSE(ptr->Pmf(-1), exp(-x_scale)*p_scale, 1e-10);
			BOOST_CHECK_CLOSE(ptr->Pmf(1), exp(-x_scale)*p_scale, 1e-10);
			BOOST_CHECK_CLOSE(ptr->Pmf(-2), exp(-4*x_scale)*p_scale, 1e-10);
			BOOST_CHECK_CLOSE(ptr->Pmf(2), exp(-4*x_scale)*p_scale, 1e-10);
		}
	}
}

BOOST_AUTO_TEST_CASE(QuantisedGaussianDistributionFixedTests)
{
	typedef QuantisedGaussianDistribution<double>::TypePtr QuantisedGaussianDistributionPtr;
	
	QuantisedGaussianDistributionPtr ptr=boost::make_shared<QuantisedGaussianDistribution<double> >(2, 1);
	
	const double values[]={
		5.3312349751089005e-5,
1.4448072588124639e-4,
3.6907845427508468e-4,
8.8902529910844796e-4,
0.002020137489946,
0.0043324483630126,
0.0087744750957384,
0.016793306448449,
0.030396361765261,
0.05208127941522,
0.084565722351336,
0.13029451713681,
0.19078695285251,
0.2659855290487,
0.35383023332728,
0.45026177516989,
0.54973822483011,
0.64616976667272,
0.7340144709513,
0.80921304714749,
0.86970548286319,
0.91543427764866,
0.94791872058478,
0.96960363823474,
0.98320669355155,
0.99122552490426,
0.99566755163699,
0.99797986251005,
0.99911097470089,
0.99963092154572,
0.99985551927412,
0.99994668765025,
0.99998146326215
	};
	
	int i=0;
	for(double x=-8;x<8;x+=0.5){
		BOOST_CHECK_CLOSE(ptr->Cdf(x), values[i], 1e-10);
		++i;
	}
}


BOOST_AUTO_TEST_CASE(QuantisedGaussianDistributionTests)
{
	srand48(0);
	
	typedef QuantisedGaussianDistribution<double>::TypePtr QuantisedGaussianDistributionPtr;
	
	for(double sigma=0.2;sigma<=5.0;sigma*=1.5){
		for(int fb=-1;fb<=3;fb++){
			QuantisedGaussianDistributionPtr ptr=boost::make_shared<QuantisedGaussianDistribution<double> >(sigma, fb);
			
			boost::math::normal_distribution<double> dist(0, sigma);
			double delta=pow(2.0, -fb);
			
			BOOST_CHECK_EQUAL(ptr->IndexFromRange(delta)-ptr->IndexFromRange(0), 1);
			
			for(int j=0;j<100;j++){
				double x=boost::math::quantile(dist, drand48());
				x=ldexp(round(ldexp(x,fb)),-fb);
				
				double pup=boost::math::cdf(dist, x+delta/2) ;
				double pdown=boost::math::cdf(dist, x-delta/2);
				
				BOOST_CHECK_CLOSE(ptr->Pmf(x), pup-pdown, 1e-10);
				BOOST_CHECK_EQUAL(ptr->Pmf(x+delta/2), 0.0);
				
				BOOST_CHECK_EQUAL(ptr->Pmf(x), ptr->Pmf(-x));
				
				BOOST_CHECK_CLOSE(ptr->Cdf(x-delta/4), pdown, 1e-10);
				BOOST_CHECK_CLOSE(ptr->Cdf(x), pup, 1e-10);
				BOOST_CHECK_CLOSE(ptr->Cdf(x+delta/4), pup, 1e-10);
				
				BOOST_CHECK_CLOSE(ptr->Cdf(x) + ptr->Cdf(-(x+delta)), 1, 1e-10);
			}			
		}
	}
	
	QuantisedGaussianDistributionPtr sig4=boost::make_shared<QuantisedGaussianDistribution<double> >(4.0, 0);
	BOOST_CHECK_CLOSE(sig4->RawMoment(2), 16.08333333333333, 1e-10);
	BOOST_CHECK_CLOSE(sig4->RawMoment(4), 776.0125000000001, 1e-10);
	BOOST_CHECK_CLOSE(sig4->RawMoment(6), 62403.00223214286, 1e-10);
	BOOST_CHECK_CLOSE(sig4->RawMoment(8), 7025313.000434028, 1e-10);
	
	QuantisedGaussianDistributionPtr sig16=boost::make_shared<QuantisedGaussianDistribution<double> >(16.0, 0);
	BOOST_CHECK_CLOSE(sig16->RawMoment(2), 256.0833333333333, 1e-10);
	BOOST_CHECK_CLOSE(sig16->RawMoment(4), 196736.0125, 1e-10);
	BOOST_CHECK_CLOSE(sig16->RawMoment(6), 2.5190404800223213e8, 1e-10);
	BOOST_CHECK_CLOSE(sig16->RawMoment(8), 4.5155894068800043e11, 1e-10);
	
	QuantisedGaussianDistributionPtr sig4p2=boost::make_shared<QuantisedGaussianDistribution<double> >(4.0, 2);
	BOOST_CHECK_CLOSE(sig4p2->RawMoment(2), 256.0833333333333/16.0, 1e-10);
	BOOST_CHECK_CLOSE(sig4p2->RawMoment(4), 196736.0125/(16.0*16.0), 1e-10);
	BOOST_CHECK_CLOSE(sig4p2->RawMoment(6), 2.5190404800223213e8 / pow(4.0,6), 1e-10);
	BOOST_CHECK_CLOSE(sig4p2->RawMoment(8), 4.5155894068800043e11 / pow(4.0,8), 1e-10);
}


BOOST_AUTO_TEST_CASE(QuantisedGaussianDistributionMprealTests)
{
	typedef QuantisedGaussianDistribution<mpfr::mpreal>::TypePtr QuantisedGaussianDistributionPtr;
	
	int prec=128;
	
	srand48(0);
	
	for(double sigma=0.2;sigma<=5.0;sigma*=1.5){
		for(int fb=-1;fb<=3;fb++){
			QuantisedGaussianDistributionPtr ptr=boost::make_shared<QuantisedGaussianDistribution<mpfr::mpreal> >(mpfr::mpreal(sigma,prec), fb);
			
			boost::math::normal_distribution<double> dist(0, sigma);
			double delta=pow(2.0, -fb);
			
			for(int j=0;j<100;j++){
				double x=boost::math::quantile(dist, drand48());
				x=ldexp(round(ldexp(x,fb)),-fb);
				
				double pup=boost::math::cdf(dist, x+delta/2) ;
				double pdown=boost::math::cdf(dist, x-delta/2);
				
				BOOST_CHECK_CLOSE(ptr->Pmf(x), pup-pdown, 1e-10);
				BOOST_CHECK_EQUAL(ptr->Pmf(x+delta/2), 0.0);
				
				BOOST_CHECK_CLOSE(ptr->Cdf(x), pup, 1e-10);
			}			
		}
	}
}

struct rng_drand48
{
	double operator()()
	{ return drand48(); }
};

BOOST_AUTO_TEST_CASE(Chi2PartitionTests)
{
	srand48(0);
	
	rng_drand48 rng;
	
	for(int fb=4;fb<=8;fb++){
		for(int k=4;k<=8192;k*=2){
			DiscreteDistribution<double>::TypePtr dist=boost::make_shared<QuantisedGaussianDistribution<double> >(1.0, fb);
			boost::shared_ptr<const Chi2Partition<double> > part(Chi2Partition<double>::CreateEqualProbability(dist, k));
			
			for(int i=8;i<=50;i++){
				double n=pow(2.0, i);
				if(n>part->MinChi2SampleSize()){
					std::vector<double> sample=part->GenerateSample<rng_drand48>(n, rng);
					Chi2Res res=part->Chi2Test(sample);
					//std::cerr<<"statistic="<<res.statistic<<", pvalue="<<res.pvalue<<"\n";
					BOOST_CHECK(std::max(res.pvalue,1-res.pvalue) > 1e-6);
				}
			}
		}
	}
}

struct rng_mpfr
{
	int m_prec;
	 gmp_randstate_t m_rng;
	
	rng_mpfr(int prec)
		: m_prec(prec)
	{
		gmp_randinit_default(m_rng);
	}
		
	~rng_mpfr()
	{
		gmp_randclear(m_rng);
	}
	
	mpfr::mpreal operator()()
	{ 
		mpfr::mpreal res(0, m_prec);
		mpfr_urandomb(res.mpfr_ptr(), m_rng);
		return res;
	}
};

BOOST_AUTO_TEST_CASE(Chi2PartitionMprealTests)
{
	int prec=128;
	rng_mpfr rng(prec);
	
	mpfr::mpreal one(1, prec);
	
	for(int fb=4;fb<=8;fb+=2){
		for(int k=4;k<=4096;k*=8){
			DiscreteDistribution<mpfr::mpreal>::TypePtr dist=boost::make_shared<QuantisedGaussianDistribution<mpfr::mpreal> >(one, fb);
			boost::shared_ptr<const Chi2Partition<mpfr::mpreal> > part(Chi2Partition<mpfr::mpreal>::CreateEqualProbability(dist, k));
			
			for(int i=16;i<=96;i+=16){
				mpfr::mpreal n=pow(one*2, i);
				if(n>part->MinChi2SampleSize()){
					std::vector<mpfr::mpreal> sample=part->GenerateSample<rng_mpfr>(n, rng);
					Chi2Res res=part->Chi2Test(sample);
					//std::cerr<<"statistic="<<res.statistic<<", pvalue="<<res.pvalue<<"\n";
					BOOST_CHECK(std::max(res.pvalue,1-res.pvalue) > 1e-6);
				}
			}
		}
	}
}

BOOST_AUTO_TEST_CASE(Chi2PartitionRegularTests)
{
	srand48(0);
	
	rng_drand48 rng;
	
	for(int fb=4;fb<=4;fb++){
		for(int k=16;k<=16;k*=2){
			DiscreteDistribution<double>::TypePtr dist=boost::make_shared<QuantisedGaussianDistribution<double> >(1.0, fb);
			boost::shared_ptr<const Chi2Partition<double> > part(Chi2Partition<double>::CreateEqualRange(dist, -4.0, 4.0, k));
			
			for(int i=8;i<=50;i++){
				double n=pow(2.0, i);
				if(n>part->MinChi2SampleSize()){
					std::vector<double> sample=part->GenerateSample<rng_drand48>(n, rng);
					Chi2Res res=part->Chi2Test(sample);
					//std::cerr<<"statistic="<<res.statistic<<", pvalue="<<res.pvalue<<"\n";
					BOOST_CHECK(std::max(res.pvalue,1-res.pvalue) > 1e-6);
				}
			}
		}
	}
}

BOOST_AUTO_TEST_CASE(Chi2RePartitionRegularTests)
{
	srand48(0);
	
	rng_drand48 rng;
	
	for(int fb=4;fb<=4;fb++){
		for(int k=16;k<=1024;k*=2){
			DiscreteDistribution<double>::TypePtr dist=boost::make_shared<QuantisedGaussianDistribution<double> >(1.0, fb);
			boost::shared_ptr<const Chi2Partition<double> > part(Chi2Partition<double>::CreateEqualRange(dist, -4.0, 4.0, k));
			
			for(int i=4;i<=50;i++){
				double n=pow(2.0, i);
				boost::shared_ptr<const Chi2Partition<double> > part_adapt=part->AdaptForSampleSize(n);
				
				if(part_adapt->NumBuckets()>1){
					part_adapt->Dump(std::cerr);
					
					std::vector<double> sample=part_adapt->GenerateSample<rng_drand48>(n, rng);
					Chi2Res res=part_adapt->Chi2Test(sample);
					std::cerr<<"statistic="<<res.statistic<<", pvalue="<<res.pvalue<<"\n";
					BOOST_CHECK(std::max(res.pvalue,1-res.pvalue) > 1e-6);
				}
			}
		}
	}
}
