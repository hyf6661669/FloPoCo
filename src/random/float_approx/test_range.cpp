#include "range.hpp"
#include <stdlib.h>

#include "range_polys.hpp"
#include "static_quantiser.hpp"

#include "Targets/Virtex5.hpp"

using namespace flopoco::random::float_approx;

extern "C" void blockSignals();

int main(int argc, char *argv[])
{
	try{
		initTool();
		invalidateRecoverEnvironment();
		blockSignals();	// Someone manages to turn this back on...
		setDisplayMode(DISPLAY_MODE_DECIMAL);
		setToolPrecision(160);
		
		std::string funcdesc="exp(x)";
		if(argc > 1){
			funcdesc=argv[1];
		}
		
		const char *sDomainStart="0.125", * sDomainEnd="8";
		if(argc> 2){
			sDomainStart=argv[2];
		}
		if(argc>3){
			sDomainEnd=argv[3];
		}
		
		int domainWF=16, rangeWF=16;
		if(argc>4){
			domainWF=atoi(argv[4]);
			rangeWF=domainWF;
		}
		if(argc>5){
			rangeWF=atoi(argv[5]);
		}
		
		int degree=3;
		if(argc>6){
			degree=atoi(argv[6]);
		}
		
		flopoco::Function f(funcdesc);
		unblockSignals();
		
		mpfr_t domainStart, domainFinish, domainEnd;
		mpfr_inits2(domainWF, domainStart, domainEnd, domainFinish, (mpfr_ptr)0);
		
		mpfr_set_str(domainStart, sDomainStart, 0, MPFR_RNDN);
		mpfr_set_str(domainEnd, sDomainEnd, 0, MPFR_RNDN);
		mpfr_set(domainFinish, domainEnd, MPFR_RNDN);
		mpfr_nextbelow(domainFinish);
		
		mpfr_fprintf(stderr, "Function=%s, domainWF=%u, rangeWF=%u\ndomain=[%Rf,%Rf)=[%Rf,%Rf]\ndegree=%u\n\n",
			funcdesc.c_str(), domainWF, rangeWF, domainStart, domainEnd, domainStart, domainFinish, degree
		); 

		Range r(f, domainWF, rangeWF, domainStart, domainFinish);
		r.dump(stdout);
		
		r.make_monotonic_or_range_flat();
		r.dump(stdout);
		
		r.flatten_domain();
		r.dump(stdout);
		
		r.flatten_range(true);
		r.dump(stdout);
		
		Range::segment_it_t curr=r.m_segments.begin();
		while(curr!=r.m_segments.end()){
			mpfr_t x, frac, fracSF, res;
			mpfr_init2(x, domainWF);
			mpfr_init2(frac, rangeWF);
			mpfr_init2(fracSF, getToolPrecision());
			mpfr_init2(res, rangeWF);
			
			int eR=mpfr_get_exp(curr->rangeStart), eD=mpfr_get_exp(curr->domainStart);
			
			sollya_node_t f=r.m_function.getSollyaNode();
			sollya_node_t sf=r.get_scaled_flat_function(eD, eR);
			
			mpfr_set(x, curr->domainStart, MPFR_RNDN);
			mpfr_mul_2si(x, x, -eD, MPFR_RNDN);
			mpfr_sub_d(x, x, 0.5, MPFR_RNDN);
			assert(mpfr_equal_p(x, curr->domainStartFrac));
			
			r.eval_scaled_flat_function(frac, x, eD, eR);
			
			mpfr_fprintf(stdout, " Start: %Rg -> %Rf -> f(.) -> %Rf\n", curr->domainStart, x, frac);
			mpfr_fprintf(stdout, "  got : %Rf\n", frac);
			evaluateFaithful(fracSF, sf, x, getToolPrecision());
			mpfr_fix(fracSF, r.m_rangeWF);
			mpfr_fprintf(stdout, "  sollya : %Rf\n", fracSF);
			// TODO : fix this!
			//assert(mpfr_equal_p(frac, fracSF));
			
			mpfr_set(res, frac, MPFR_RNDN);
			mpfr_add_d(res, res, 0.5, MPFR_RNDN);
			mpfr_mul_2si(res, res, eR, MPFR_RNDN);
			
			//mpfr_fprintf(stdout, "   result, res=%Rf, rangeStart=%Rf\n", res, curr->rangeStart);
			assert(mpfr_equal_p(res, curr->rangeStart));
			
			
			mpfr_set(x, curr->domainFinish, MPFR_RNDN);
			mpfr_mul_2si(x, x, -eD, MPFR_RNDN);
			mpfr_sub_d(x, x, 0.5, MPFR_RNDN);
			assert(mpfr_equal_p(x, curr->domainFinishFrac));
			
			r.eval_scaled_flat_function(frac, x, eD, eR);
			
			mpfr_fprintf(stdout, " Finish: %Rg -> %Rf -> f(.) -> %Rf\n", curr->domainFinish, x, frac);
			mpfr_fprintf(stdout, "  got : %Rf\n", frac);
			evaluateFaithful(fracSF, sf, x, getToolPrecision());
			mpfr_fix(fracSF, r.m_rangeWF);
			mpfr_fprintf(stdout, "  sollya : %Rf\n", fracSF);
			// TODO : Why does this not work?
			//assert(mpfr_equal_p(frac, fracSF));
			
			mpfr_set(res, frac, MPFR_RNDN);
			mpfr_add_d(res, res, 0.5, MPFR_RNDN);
			mpfr_mul_2si(res, res, eR, MPFR_RNDN);
			
			//mpfr_fprintf(stdout, "   result, res=%Rf, rangeFinish=%Rf\n", res, curr->rangeFinish);
			assert(mpfr_equal_p(res, curr->rangeFinish));
			
			assert(curr->isRangeFlat);
			assert(curr->isDomainFlat);
			
			++curr;
			mpfr_clears(x, res, frac, fracSF, (mpfr_ptr)0);
		}
		
		fprintf(stderr, "Starting RangePoly\n");
		RangePolys rp(r, degree);
		
		fprintf(stderr, "\n\nSplittin to error of 2^(-rangeWF)\n");
		rp.split_to_error(pow(2.0, -rangeWF));
		r.dump(stdout);
		int guard=0;
		
		fprintf(stderr, "\n\nCreating faithful with %d guard bits\n", guard);
		rp.calc_faithful_fixed_point(guard);
		r.dump(stdout);
		
		fprintf(stderr, "\n\nBuilding concrete with %d guard bits\n", guard);
		rp.build_concrete(guard);
		r.dump(stdout);
		
		::flopoco::verbose=FULL;
		
		flopoco::Virtex5 target;
		flopoco::PolynomialEvaluator *pe=rp.make_polynomial_evaluator(&target);
		
		StaticQuantiser *sq=rp.make_static_quantiser(&target);
		sq->outputVHDL(std::cout);
		
		
		mpfr_clears(domainStart, domainFinish, (mpfr_ptr)0);
		fprintf(stderr, "Exiting scope\n");		
	}catch(std::exception &e){
		fprintf(stderr, "Caught C++ exception : %s\n", e.what());
		return 1;
	}catch(...){
		fprintf(stderr, "Caught uknown exception.\n");
		return 1;
	}
		
	return 0;
}
