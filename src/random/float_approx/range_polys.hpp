#ifndef flopoco_random_float_approx_range_polys_hpp
#define flopoco_random_float_approx_range_polys_hpp

#include "range.hpp"

#include <boost/smart_ptr.hpp>
#include <boost/any.hpp>

namespace flopoco
{
namespace random
{
namespace float_approx
{

typedef boost::shared_ptr<sollya_node> sollya_node_ptr_t;

sollya_node_ptr_t make_shared_sollya(sollya_node_t n)
{ return sollya_node_ptr_t(n, free_memory); }

class RangePolys
{
public:	
	typedef Range::segment_it_t segment_it_t;
private:
	
	Range &m_range;
	unsigned m_degree;

	void update_segment(segment_it_t curr)
	{
		if(curr->has_property("minimax") && curr->has_property("minimax_error")){
			fprintf(stderr, "Skip update\n");
			return;
		}
		
		sollya_node_ptr_t poly=make_shared_sollya(curr->minimax(m_degree));
		sollya_node_t func=curr->get_scaled_flat_function(); // doesn't need to be freed
		
		sollya_node_t diff=makeSub(copyTree(poly.get()), copyTree(func));
		
		mpfr_t error;
		mpfr_init2(error, getToolPrecision());
		
		uncertifiedInfnorm(error, diff, curr->domainStartFrac, curr->domainFinishFrac, 71, getToolPrecision()); 
		free_memory(diff);
		
		fprintf(stderr, "Updated error = %lf\n", mpfr_get_d(error, MPFR_RNDN));
		
		curr->properties["minimax"]=poly;
		curr->properties["minimax_error"]=(double)mpfr_get_d(error, MPFR_RNDN);
		mpfr_clear(error);
	}
public:
	RangePolys(Range &range, unsigned degree)
		: m_range(range)
		, m_degree(degree)
	{
		segment_it_t it=m_range.m_segments.begin();
		while(it!=m_range.m_segments.end()){
			update_segment(it);
			++it;
		}
	}
		
	void split_to_error(double error)
	{
		segment_it_t it=m_range.m_segments.begin();
		while(it!=m_range.m_segments.end()){
			update_segment(it);
			
			if( boost::any_cast<double>(it->properties["minimax_error"]) >= error){
				mpfr_fprintf(stderr, "Splitting [%Rf,%Rf]",
					it->domainStart, it->domainFinish
				);
				fprintf(stderr, " error=%lf\n", boost::any_cast<double>(it->properties["minimax_error"]));
				m_range.split_segment(it);
				update_segment(it);
			}			
			assert(it->has_property("minimax_error"));
			if(boost::any_cast<double>(it->properties["minimax_error"]) < error)
				++it;
		}
	}
};

}; // float_approx
};	// random
}; // flopco

#endif
