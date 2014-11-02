#include "correct_distribution.hpp"

#include "random/distributions/gaussian_distribution.hpp"
#include "random/distributions/make_table_approximation.hpp"
#include "random/distributions/clt_distribution.hpp"

#include "random/utils/anneal_gsl.hpp"

#include "gsl/gsl_cdf.h"

using namespace flopoco::random;

class HadamardEnergyFunction
	: public VectorEnergyFunction
{
public:
	unsigned m_n;
	int m_fb;
	unsigned m_k;

	std::vector<double> m_guess;

	int m_expectedBase;
	std::vector<double> m_expectedPdf, m_expectedCdf;

	std::vector<double> m_buckets1024; // Equally probability bucket boundaries

	int ToInt(double v) const
	{ return ldexp(v, m_fb); }
	
	double Fix(double v) const
	{ return ldexp(round(ldexp(v,m_fb)), -m_fb); }

	HadamardEnergyFunction(unsigned n, double fb, int k)
		: m_n(n)
		, m_fb(fb)
		, m_k(k)
	{
		double scale=1.0/sqrt(k);
		
		double std=0.0;
		m_guess.resize(n/2);
		for(int i=0;i<(int)n/2;i++){
			double u=(i+0.5)/n;
			m_guess[i]=gsl_cdf_ugaussian_Pinv(u);
			std+=2*m_guess[i]*m_guess[i];
		}
		std=sqrt(std/n);
		for(int i=0;i<(int)n/2;i++){
			m_guess[i]=std::fabs(Fix(m_guess[i]/std*scale));
		}
		for(int i=n/2-1;i>0;i--){
			m_guess[i]=m_guess[i]-m_guess[i-1];
		}
		
		m_expectedBase=-16*(1<<m_fb);
		m_expectedPdf.resize(-2*m_expectedBase+1);
		m_expectedCdf.resize(-2*m_expectedBase+1);
		
		double pdfScale=1.0/(pow(2.0,m_fb)*sqrt(2*3.1415926535897932384626433832795));
		for(int i=0;i<(int)m_expectedPdf.size();i++){
			double x=ldexp(i+m_expectedBase, -m_fb);
			m_expectedCdf[i] = gsl_cdf_ugaussian_P(x+ldexp(1.0,-m_fb));
			m_expectedPdf[i] = exp(-x*x/2) * pdfScale;
		}
		
		m_buckets1024.resize(1025);
		for(int i=0;i<=1024;i++){
			m_buckets1024[i]=gsl_cdf_ugaussian_Pinv(i/1024.0);
		}
	}

	virtual state_t Step(const gsl_rng *r, const state_t &src, double step_size) const
	{
		std::vector<double> res(src);
		for(int i=0;i<(int)src.size();i++){
			double delta=(gsl_rng_uniform(r)-gsl_rng_uniform(r)) * step_size;
			double tmp=Fix(res[i]+delta);
			//fprintf(stderr, "%lg\n", tmp-res[i]);
			res[i]=tmp;
		}
		return res;
	}
	
	struct metrics_t
	{
		double pdfChi2;
		double cdfChi2;
		double ksStat;
		double adStat;
		double stddev;
		double kurtosis;
		double baseStddev;
		double baseKurt;
		double buckets1024;
	};
	
	metrics_t CalcMetrics(double gotBase, int nGot, const double *gotPdf, const double *truePdf, const double *trueCdf) const
	{
		metrics_t res;
		res.pdfChi2=0;
		res.cdfChi2=0;
		res.ksStat=0;
		res.adStat=0;
		res.stddev=0;
		res.kurtosis=0;
		res.buckets1024=0;
		double oc=0;
		double x=gotBase, xstep=ldexp(1,-m_fb);
		x=x-xstep;
		for(int i=0;i<=nGot/2;i++){
			x+=xstep;
			double op=gotPdf[i], ep=truePdf[i];
			oc += op;
			double ec=trueCdf[i];
			res.pdfChi2 += 2*(op-ep)*(op-ep)/ep;
			res.cdfChi2 += 2*(oc-ec)*(oc-ec)/ec;
			res.ksStat= std::max(res.ksStat, std::fabs(oc-ec));
			res.adStat += 2*(oc-ec)*(oc-ec) / (ec*(1-ec)); 
			res.stddev += 2*x*x * op;
			res.kurtosis += 2*x*x*x*x * op;
			if(x==0){
				res.pdfChi2 -= (op-ep)*(op-ep)/ep;
				res.adStat -= (oc-ec)*(oc-ec) / (ec*(1-ec)); 
			}
		}
		assert(x==0);
		res.kurtosis /= res.stddev*res.stddev;
		res.stddev=sqrt(res.stddev);
		
		const double *pUpper=&m_buckets1024[1];
		double tAcc=0, cAcc=0;
		x=gotBase-xstep;
		for(int i=0;i<nGot;i++){
			x+=xstep;
			double p=gotPdf[i];
			while(x>*pUpper){
				tAcc+= pow(cAcc-1.0/1024,2)*1024;
				cAcc=0;
				pUpper++;
			}
			cAcc+=p;
		}
		tAcc+= pow(cAcc-1.0/1024,2)*1024;
		res.buckets1024=tAcc;
		
		return res;
	}

	metrics_t EnergyFull(const std::vector<double> &s) const
	{
		std::vector<double> tmp(s);
		double baseStddev=0, baseKurt=0;
		double x=0;
		for(int i=0;i<s.size();i++){
			x+=s[i];
			tmp.push_back(-x);
			tmp.push_back(x);
			baseStddev+=2*x*x;
			baseKurt+=2*x*x*x*x;
		}
		std::sort(tmp.begin(), tmp.end());
		baseStddev/=(tmp.size());
		baseKurt/=(tmp.size()*baseStddev);
		baseStddev=sqrt(baseStddev);
		
		int minv=ToInt(tmp.front());
		int maxv=ToInt(tmp.back());
		
		std::vector<double> starter(maxv-minv+1, 0.0);
		for(int i=0;i<(int)tmp.size();i++){
			starter[ToInt(tmp[i])-minv]+=1.0/tmp.size();
		}
		
		std::vector<double> full=self_convolve(starter, m_k);
		int fullBase=minv*m_k;
		
		//double acc=sum(full.begin(), full.end());
		//fprintf(stderr, "  sum(full=1+%lg\n", 1-acc);
		
		const double *pExpectedPdf=&m_expectedPdf[fullBase-m_expectedBase];
		const double *pExpectedCdf=&m_expectedCdf[fullBase-m_expectedBase];
		
		assert(fullBase+full.size()-1 == -fullBase);
		metrics_t mm= CalcMetrics(ldexp(fullBase,-m_fb), full.size(), &full[0], pExpectedPdf, pExpectedCdf);
		mm.baseStddev=baseStddev;
		mm.baseKurt=baseKurt;
		
		/*std::vector<std::pair<double,double> > ww(full.size());
		for(int i=0;i<full.size();i++){
			ww[i]=std::make_pair(ldexp(fullBase+i,-m_fb), full[i]);
		}
		double trueMean=sum_weighted_powers(ww.begin(), ww.end(), 1);
		double trueStd=sqrt(sum_weighted_powers(ww.begin(), ww.end(), 2));
		fprintf(stderr, "  trueMean=%lg, trueStd=1+%lg, gotStd=1+%lg\n", trueMean, 1-trueStd, 1-mm.stddev);*/
		return mm;
	}
	
	virtual double Energy(const std::vector<double> &s) const
	{
		metrics_t m=EnergyFull(s);
		return m.buckets1024;
	}
	
	virtual void Print(const state_t &a) const
	{
		metrics_t mm=EnergyFull(a);
		fprintf(stderr, "pchi2=%lg, cchi2=%lg, ks=%lg, ad=%lg, sttdev=%.8lg kurt=%lg, baseStd=%.8lg, baseKurt=%lg, buckets1024=%lg\n",
			mm.pdfChi2, mm.cdfChi2, mm.ksStat, mm.adStat, mm.stddev, mm.kurtosis, mm.baseStddev, mm.baseKurt, mm.buckets1024
		);
	}
	
	virtual unsigned IterationsPerTemperature() const
	{ return 10; }
	
	virtual double InitialTemperature() const
	{ return 0.25; }
	
	virtual double MinimumTemperature() const
	{ return 1e-10; }
	
	virtual double MaxStepSize() const
	{ return pow(2.0, -(m_fb-3)); }

	virtual double MinStepSize() const
	{ return pow(2.0, -m_fb); }
	
	virtual double StepScale() const
	{ return 0.25; }
	
	virtual double CoolingFactor() const
	{ return 1.01; }
};

int main(int argc, char *argv[])
{
	try{
		int n=32;
		int k=4;
		int s_w=8, t_w=10;
		std::vector<double> curr;
		
		while(s_w<=t_w){
			s_w++;
		
			HadamardEnergyFunction ef(n, s_w, k);
			if(curr.size()!=0){
				double eg=ef.Energy(ef.m_guess), eo=ef.Energy(curr);
				std::cerr<<" eGuess="<<eg<<", eCurr="<<eo<<"\n";
				if(eg > eo ){
					ef.m_guess=curr;
				}
			}
		
			std::cerr<<"init = "<<ef.Energy(ef.m_guess)<<"\n";
		
			std::pair<std::vector<double>,double> res=ef.Execute(ef.m_guess, false);
			double acc=0;
			for(int i=res.first.size()-1;i>=0;i--){
				acc-=res.first[i];
				fprintf(stdout, "%d, %.16lg\n", i, acc);
			}
			acc=0;
			for(int i=0;i<res.first.size();i++){
				acc+=res.first[i];
				fprintf(stdout, "%d, %.16lg\n", i+res.first.size(), acc);
			}
			ef.Print(res.first);
			
			curr=res.first;
		}
			
		return 0;
	}catch(std::exception &e){
		std::cerr<<"Caught Exception : "<<e.what()<<"\n";
		return 1;
	}catch(...){
		std::cerr<<"Caught unexpected exception.";
		return 1;
	}
}
