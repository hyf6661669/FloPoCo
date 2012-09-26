#ifndef flopoco_random_float_approx_range_polys_hpp
#define flopoco_random_float_approx_range_polys_hpp

#include "range.hpp"

#include <boost/smart_ptr.hpp>
#include <boost/any.hpp>
#include <boost/format.hpp>
#include <boost/utility.hpp>

#include "concrete_poly.hpp"


namespace flopoco
{
namespace random
{
namespace float_approx
{

/*
typedef boost::shared_ptr<sollya_node> sollya_node_ptr_t;

sollya_node_ptr_t make_shared_sollya(sollya_node_t n)
{ return sollya_node_ptr_t(n, free_memory); }
*/
	
typedef struct{
	sollya_node_t _x;
	
	sollya_node_t get()
	{ return _x; }
}sollya_node_ptr_t;
	
sollya_node_ptr_t make_shared_sollya(sollya_node_t n)
{ sollya_node_ptr_t x; x._x=n; return x; }


MPFRVec getPolyCoeffs(sollya_node_t poly)
{
	unsigned degree=getDegree(poly);
	MPFRVec res(degree+1, getToolPrecision());
	for(unsigned i=0;i<=degree;i++){
		sollya_node_t c=getIthCoefficient(poly, i);
		evaluateConstantExpression(res[i], c, getToolPrecision());
		//free_memory(c);
	}
	return res;
}

class RangePolys
{
public:	
	typedef Range::segment_it_t segment_it_t;
private:
	
	Range &m_range;
	unsigned m_degree;

	double infNorm(sollya_node_t func, sollya_node_t poly, mpfr_t a, mpfr_t b, unsigned points)
	{
		sollya_node_t diff=makeSub(copyTree(func), copyTree(poly));
		
		mpfr_t error;
		mpfr_init2(error, getToolPrecision());
		
		uncertifiedInfnorm(error, diff, a, b, points, getToolPrecision()); 
		// TODO : This causes corruption later on, but I don't see why, as both inputs are copied
		//free_memory(diff);
		
		double res=mpfr_get_d(error, MPFR_RNDD);	// round down
		mpfr_clear(error);
		
		return res;
	}

	void update_segment_minimax(segment_it_t curr)
	{
		if(curr->has_property("minimax") && curr->has_property("minimax_error"))
			return;
		
		sollya_node_ptr_t poly=make_shared_sollya(curr->minimax(m_degree));
		sollya_node_t func=curr->get_scaled_flat_function(); // doesn't need to be freed
		
		double error=infNorm(func, poly.get(), curr->domainStartFrac, curr->domainFinishFrac, 71);
		
		curr->properties["minimax"]=poly;
		curr->properties["minimax_error"]=error;
	}
public:
	RangePolys(Range &range, unsigned degree)
		: m_range(range)
		, m_degree(degree)
	{
		segment_it_t it=m_range.m_segments.begin();
		while(it!=m_range.m_segments.end()){
			update_segment_minimax(it);
			++it;
		}
	}
		
	void split_to_error(double error)
	{
		segment_it_t it=m_range.m_segments.begin();
		while(it!=m_range.m_segments.end()){
			update_segment_minimax(it);
			
			if( boost::any_cast<double>(it->properties["minimax_error"]) >= error){
				mpfr_fprintf(stderr, "Splitting [%Rf,%Rf]",
					it->domainStart, it->domainFinish
				);
				fprintf(stderr, " error=%lf\n", boost::any_cast<double>(it->properties["minimax_error"]));
				m_range.split_segment(it);
				update_segment_minimax(it);
			}			
			assert(it->has_property("minimax_error"));
			if(boost::any_cast<double>(it->properties["minimax_error"]) < error)
				++it;
		}
	}
	
	// Find fixed-point polynomial coefficients such that the result is faithful on rangeWF+guard bits
	void calc_faithful_fixed_point(segment_it_t curr, int guard)
	{
		std::string pnCoeffs=boost::str(boost::format("minimax_fixed%1%_coeffs")%guard);
		
		if(curr->has_property(pnCoeffs))
			return;
		
		update_segment_minimax(curr);
		
		double tol=pow(2.0, -m_range.m_rangeWF-guard-1);
		
		if(boost::any_cast<double>(curr->properties["minimax_error"]) > tol/2)
			throw std::runtime_error("calc_faithful_fixed_point - poly is not accurate enough to be faithful.");
		
		// For now this is just dumb, and uses too many bits
		std::vector<int> widths(m_degree+1, m_range.m_rangeWF+guard+1);
		
		sollya_node_t fp=curr->fpminimax(widths, boost::any_cast<sollya_node_ptr_t>(curr->properties["minimax"]).get() );
		
		double error=infNorm(curr->get_scaled_flat_function(), fp, curr->domainStartFrac, curr->domainFinishFrac, 71);
		
		MPFRVec coeffs=getPolyCoeffs(fp);
		
		// free_memory(fp)
		
		if(error > tol)
			throw std::runtime_error("calc_faithful_fixed_point - fixed-point poly was not faithful.");
		
		std::string pnWidths=boost::str(boost::format("minimax_fixed%1%_widths")%guard);
		std::string pnError=boost::str(boost::format("minimax_fixed%1%_error")%guard);
		
		curr->properties[pnCoeffs]=coeffs;
		curr->properties[pnWidths]=widths;
		curr->properties[pnError]=pnError;
	}
	
	void calc_faithful_fixed_point(int guard)
	{
		// We split it so that the minimax is twice as accurate as needed. This means we use too
		// many segments, and should be optimised.
		split_to_error(pow(2.0, -m_range.m_rangeWF-guard-2));
		
		segment_it_t curr=m_range.m_segments.begin();
		while(curr!=m_range.m_segments.end()){
			calc_faithful_fixed_point(curr, guard);
			++curr;
		}
	}
	
	MPFRVec m_concretePartition;
	std::vector<ConcretePoly> m_concretePolys;
	std::vector<int> m_concreteExp;
	
	void build_concrete(int guard)
	{
		std::string pnCoeffs=boost::str(boost::format("minimax_fixed%1%_coeffs")%guard);
		std::string pnWidths=boost::str(boost::format("minimax_fixed%1%_widths")%guard);
		
		fprintf(stderr, "segments = %d\n", m_range.m_segments.size());
		m_concretePartition.resize(m_range.m_segments.size(), m_range.m_domainWF);
		fprintf(stderr, "  len = %d\n", m_concretePolys.size());
		m_concretePolys.resize(0);
		m_concreteExp.resize(0);
		
		segment_it_t curr=m_range.m_segments.begin();
		for(unsigned i=0;i<m_range.m_segments.size();i++){
			assert(curr->has_property(pnCoeffs));
			MPFRVec coeffs=boost::any_cast<MPFRVec>(curr->properties[pnCoeffs]);
			std::vector<int> widths=boost::any_cast<std::vector<int> >(curr->properties[pnWidths]);
			
			MPFRVec sRange(2, m_range.m_domainWF);
			mpfr_set(sRange[0], curr->domainStartFrac, MPFR_RNDN);
			mpfr_set(sRange[1], curr->domainFinishFrac, MPFR_RNDN);
			
			ConcretePoly cp(coeffs, sRange, m_range.m_domainWF, widths, m_range.m_rangeWF+guard+1, m_range.m_rangeWF);
			
			fprintf(stderr, "  cp[%d], len=%d\n", i, m_concretePartition.size());
			mpfr_set(m_concretePartition[i], curr->domainFinish, MPFR_RNDN);
			m_concretePolys.push_back(cp);
			m_concreteExp.push_back(mpfr_get_exp(curr->rangeStart));
			
			++curr;
		}
		
		fprintf(stderr, "Done build concrete.\n");
	}
	
	unsigned eval_concrete(mpfr_t res, mpfr_t x)
	{
		int i;
		for(i=0;i<m_concretePartition.size();i++){
			if(mpfr_lessequal_p(x, m_concretePartition[i]))
				break;
		}
		if(i==m_concretePartition.size())
			throw std::invalid_argument("eval_concrete - out of range.");
		
		ConcretePoly &cp=m_concretePolys[i];
		
		mpfr_t rx, rres;
		
		mpfr_init2(rx, m_range.m_domainWF);
		mpfr_set(rx, x, MPFR_RNDN);
		mpfr_mul_2si(rx, rx, -mpfr_get_exp(x), MPFR_RNDN);
		mpfr_sub_d(rx, rx, 0.5, MPFR_RNDN);
		
		mpfr_init2(rres, m_range.m_rangeWF);
		cp.eval(rres, rx);
		
		mpfr_add_d(rres, rres, 0.5, MPFR_RNDN);
		mpfr_mul_2si(rres, rres, m_concreteExp[i], MPFR_RNDN);
		
		mpfr_set(res, rres, MPFR_RNDN);
		
		return i;
	}
	
	void dump_concrete(FILE *dst)
	{
		
		for(int i=0;i<m_concretePolys.size();i++){
			m_concretePolys[i].dump(dst, "  ");
		}
		fprintf(stderr, "Total polys=%u\n", m_concretePolys.size());
	}
	
	int ulps(mpfr_t got, mpfr_t want)
	{
		if(mpfr_equal_p(got, want))
			return 0;
		
		int dir=mpfr_greater_p(want, got) ? +1 : -1;
		
		mpfr_t tmp;
		mpfr_init_copy(tmp, got);
		int res=0;
		
		while(true){
			++res;
			if(dir==+1){
				mpfr_nextabove(tmp);
				if(mpfr_lessequal_p(want, tmp))
					break;
			}else{
				mpfr_nextbelow(tmp);
				if(mpfr_greaterequal_p(want, tmp))
					break;
			}
		}
			
		mpfr_clear(tmp);
		
		return res;
	}
	
	void exhaust_concrete()
	{
		mpfr_t x, got, want, actual, err, tol;
		
		mpfr_init2(got, m_range.m_rangeWF);
		mpfr_init2(want, m_range.m_rangeWF);
		mpfr_init2(actual, getToolPrecision());
		mpfr_init2(err, getToolPrecision());
		mpfr_init2(tol, getToolPrecision());
		mpfr_init_copy(x, m_range.m_segments.begin()->domainStart);
		
		mpfr_set_si_2exp(tol, 1, -m_range.m_rangeWF, MPFR_RNDN);
		
		int start=time(NULL);
		int delta_print=2;
		int next_print=start+delta_print;
		
		while(mpfr_lessequal_p(x, boost::prior(m_range.m_segments.end())->domainFinish)){
			unsigned seg=eval_concrete(got, x);
			m_range.m_function.eval(actual, x);
			mpfr_set(want, actual, MPFR_RNDN);
			
			mpfr_sub(err, got, actual, MPFR_RNDN);
			mpfr_div(err, err, actual, MPFR_RNDN);
			mpfr_abs(err, err, MPFR_RNDN);
			
			int u=ulps(got, want);
			
			if(u>1){
				mpfr_fprintf(stdout, "%u, %Rg, %Rg, %Rg, %Rg, %Rg, %s, %d\n", seg, x, got, want, actual, err, mpfr_greater_p(err, tol) ? "FAIL" : "pass", u);
			}
			/*if(!mpfr_equal_p(got, want)){
				if(mpfr_greater_p(err, tol)){
					fprintf(stderr, "Fail\n");
					exit(1);
				}
			}*/
			
			mpfr_nextabove(x);
			
			if(time(NULL) > next_print){
				mpfr_fprintf(stderr, "  curr=%Rg, range=[%Rg,%Rg]\n", x, m_range.m_segments.begin()->domainStart, boost::prior(m_range.m_segments.end())->domainFinish);
				delta_print=std::min(60, delta_print*2);
				next_print=time(NULL)+delta_print;
			}
		}
		
		mpfr_clears(x, got, want, (mpfr_ptr)0);
	}
};

}; // float_approx
};	// random
}; // flopco

#endif
