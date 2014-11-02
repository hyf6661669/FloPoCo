#ifndef random__utils_anneal_gsl_hpp
#define random__utils_anneal_gsl_hpp

#include "gsl_utils.hpp"

#include "gsl/gsl_siman.h"

#include <stdexcept>
#include <iostream>
#include <cmath>

namespace flopoco
{
namespace random
{
namespace detail
{	
/*
	class TF{
		typedef ... state_t; // must be copyable and copy-constructible
	
		double Energy(const state_t &a);
		double Distance(const state_t &a, const state_t &b) const;
		state_t Step(gsl_rng *r, const state_t &src, double step_size) const;
		void Print(const state_t &a);	
	
		double BoltzmanConstant() const;
		double InitialTemperature() const;
		double CoolingFactor() const;
		double MinimumTemperature() const;
		
		double MaxStepSize() const;
		unsigned TriesPerStep() const;
		unsigned IterationsPerTemperature() const;
	};

*/
template<class TF>
class AnnealGSLImpl
{
private:
	TF &m_f;
	AnnealGSLImpl *m_me;

	typedef typename TF::state_t state_t;

	struct state_wrapper_t
	{
		AnnealGSLImpl *parent;
		state_t state;
		double energy;
	};
	
	static state_wrapper_t *get(void *s)
	{
		assert( ((state_wrapper_t*)s)->parent->m_me == ((state_wrapper_t*)s)->parent);
		return (state_wrapper_t*)s;
	}

	static void CopyFunc(void *src, void *dst)
	{ *get(dst) = *get(src); }
	
	static void *CopyConstructFunc(void *src)
	{ return (void*)new state_wrapper_t(*get(src)); }
	
	static void DestroyFunc(void *src)
	{ delete get(src); }
	
	static void StepFunc(const gsl_rng *r, void *xp, double step_size)
	{
		state_t ns = get(xp)->parent->m_f.Step(r, get(xp)->state, step_size); 
		if(get(xp)->state == ns){
			return;
		}else{
			get(xp)->state = ns;
			get(xp)->energy=-DBL_MAX;
		}
	}
	
	static double DistanceFunc(void *a, void *b)
	{ return get(a)->parent->m_f.Distance(get(a)->state, get(b)->state); }
	
	static void PrintFunc(void *src)
	{ get(src)->parent->m_f.Print(get(src)->state); }
	
	static double EnergyFunc(void *src)
	{
		if(get(src)->energy==-DBL_MAX){
			get(src)->energy = get(src)->parent->m_f.Energy(get(src)->state);
		}
		return get(src)->energy;
	}
	
public:
	AnnealGSLImpl(TF &f)
		: m_f(f)
		, m_me(0)
	{}
	
	state_t Go(const state_t &init, bool quiet=true)
	{
		gsl_rng *rng=gsl_rng_alloc (gsl_rng_default);
		
		double currEnergy=m_f.Energy(init);
		
		state_wrapper_t res={
			this,
			init,
			currEnergy
		};
		fprintf(stderr, "init.size=%d\n", init.size());
		gsl_siman_params_t params;
		params.k=m_f.BoltzmanConstant();
		params.t_initial=std::max(m_f.MinimumTemperature(), std::min(m_f.InitialTemperature(), currEnergy*10));
		params.mu_t=m_f.CoolingFactor();
		params.t_min=m_f.MinimumTemperature();
		params.step_size=m_f.MaxStepSize();
		params.n_tries=m_f.TriesPerStep();
		params.iters_fixed_T=m_f.IterationsPerTemperature();
		
		while(params.step_size>=m_f.MinStepSize()){
			m_me=this;
			gsl_siman_solve (
				rng,
				(void *)&res, /*x0_p,*/
				EnergyFunc, //gsl_siman_Efunc_t Ef,
				StepFunc, //gsl_siman_step_t take_step,
				DistanceFunc, //gsl_siman_metric_t distance,
				quiet ? 0 : PrintFunc, //gsl_siman_print_t print_position,
				CopyFunc, //gsl_siman_copy_t copyfunc,
				CopyConstructFunc, //gsl_siman_copy_construct_t copy_constructor,
				DestroyFunc, //gsl_siman_destroy_t destructor,
				0, /*size_t element_size,*/
				params
			);
			if(m_f.Energy(res.state) < m_f.TargetError())
				break;
			params.step_size *= m_f.StepScale();
			params.t_initial=m_f.StepScale();
		}
		
		gsl_rng_free(rng);
		rng=0;
		
		return res.state;
	}
};

}; // detail

class VectorEnergyFunction
{
public:
	virtual ~VectorEnergyFunction()
	{}
	
	typedef std::vector<double> state_t; // must be copyable and copy-constructible
		
	virtual double Energy(const state_t &s) const=0;

	static double EuclidianDistance(const state_t &a, const state_t &b)
	{
		double acc=0;
		for(size_t i=0;i<a.size();i++){
			double d=a[i]-b[i];
			acc+=d*d;
		}
		return sqrt(acc);
	}
	
	static double ManhattanDistance(const state_t &a, const state_t &b)
	{
		double acc=0;
		for(size_t i=0;i<a.size();i++){
			acc+=std::fabs(a[i]-b[i]);
		}
		return acc;
	}
	
	virtual double Distance(const state_t &a, const state_t &b) const
	{ return EuclidianDistance(a,b); }
	
	static state_t AbsoluteTriangleStep(const gsl_rng *r, const state_t &src, double step_size)
	{
		state_t res(src);
		for(size_t i=0;i<res.size();i++){
			double delta=gsl_rng_uniform(r)-gsl_rng_uniform(r);
			res[i] += delta*step_size;
		}
		return res;
	}
	
	static state_t RelativeTriangleStep(const gsl_rng *r, const state_t &src, double step_size)
	{
		state_t res(src);
		for(size_t i=0;i<res.size();i++){
			double delta=exp( (gsl_rng_uniform(r)-gsl_rng_uniform(r)) * step_size);
			res[i] *= delta;
		}
		return res;
	}
	
	virtual state_t Step(const gsl_rng *r, const state_t &src, double step_size) const
	{ return AbsoluteTriangleStep(r, src, step_size); }
	
	virtual void Print(const state_t &a) const
	{
		fprintf(stdout, "<");
		for(size_t i=0;i<a.size();i++){
			fprintf(stdout, i==0?"%lg":",%lg", a[i]);
		}
		fprintf(stdout, ">");
	}

	virtual double BoltzmanConstant() const
	{ return 1; }
	
	virtual double InitialTemperature() const
	{ return 0.1; }
	
	virtual double CoolingFactor() const
	{ return 1.01; }
	
	virtual double MinimumTemperature() const
	{ return 1e-6; }
	
	virtual double MaxStepSize() const
	{ return 1.0; }
	
	virtual double MinStepSize() const
	{ return pow(2.0,-32); }
	
	virtual double StepScale() const
	{ return 0.5; }
	
	virtual unsigned TriesPerStep() const
	{ return 10; }
	
	virtual unsigned IterationsPerTemperature() const
	{ return 50; }
	
	virtual double TargetError() const
	{ return 1e-20; }
	
	std::pair<state_t,double> Execute(const state_t &init, bool quiet=true)
	{
		detail::AnnealGSLImpl<VectorEnergyFunction> impl(*this);
		std::vector<double> res=impl.Go(init, quiet);
		return std::make_pair(res, Energy(res));
	}
};


}; // random
}; // flopoco


#endif
