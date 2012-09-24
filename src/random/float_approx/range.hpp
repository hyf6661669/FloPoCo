#ifndef flopoco_random_float_approx_range_hpp
#define flopoco_random_float_approx_range_hpp

#include <list>
#include <map>
#include <stdio.h>
#include <mpfr.h>
#include <assert.h>

#include <boost/any.hpp>

#include "FixedPointFunctions/Function.hpp"

namespace flopoco
{
namespace random
{
namespace float_approx
{

// initialise and copy same value (with same precision)
void mpfr_init_copy(mpfr_t dst, mpfr_t src);
	
// Fix to this many fractional bits
void mpfr_fix(mpfr_t x, int bits);
	
struct Range;

struct Segment
{
private:
	Segment();
public:
	
	Range *parent;
	
	// Start and finish points are inclusive
	mpfr_t domainStart, domainFinish;
	mpfr_t rangeStart, rangeFinish;

	mpfr_t domainStartFrac, domainFinishFrac;
	
	// Return true if the range/domain are all within one exponent value (i.e. one binade)
	bool isRangeFlat, isDomainFlat;

	typedef std::map<std::string,boost::any> property_map_t;
	property_map_t properties;
	
	void set_domain(mpfr_t _domainStart, mpfr_t _domainFinish);
	
	Segment(const Segment &o);
	Segment(Range *_parent, mpfr_t _domainStart, mpfr_t _domainFinish);
	
	Segment &operator=(const Segment &o);
	
	~Segment();
	
	void dump(FILE *dst) const;

	// Compute minimax with "real" coefficients
	sollya_node_t minimax(unsigned degree);
	
	// Compute minimax with fixed-point coefficients
	sollya_node_t fpminimax(sollya_chain_t monomials, sollya_chain_t formats);
	
	// The returned node should not be freed (it is cached)
	sollya_node_t get_scaled_flat_function();
	
	bool has_property(std::string name) const;
};

struct Range
{
	Function &m_function;
	int m_domainWF;
	int m_rangeWF;
	
	// upper (exclusive) bound for integer value of the domain and range, not including the implicit bit
	mpfr_t m_domainFractionEnd, m_rangeFractionEnd;
	
	typedef std::list<Segment> segment_list_t;
	typedef segment_list_t::iterator segment_it_t;
	segment_list_t m_segments;
	
	typedef std::map<std::pair<int,int>,sollya_node_t> flat_function_cache_t;
	flat_function_cache_t m_flatFunctions;
	
	Range(Function &f, int domainWF, int rangeWF, mpfr_t domainStart, mpfr_t domainFinish);
	
	~Range();
	
	/* Split the segment, reducing the range of the original and creating a new one.
		If the original range was  [start,finish], then its range will change to [start,newFinish],
		and the next segment will have [next(newFinish),finish]
	*/
	void split_segment(segment_it_t src, mpfr_t newFinish);
	
	// Splits the segment into two equal parts
	void split_segment(segment_it_t src);
	
	// Keep splitting in domain until all segments are flat
	void flatten_domain();
	
	/* Requires: binade(f(start)) < binade(f(finish))
		Ensures: start<=newFinish < finish  and  binade(f(start))==binade(f(newFinish))  and  binade(f(newFinish))<binade(f(next(newFinish)))*/
	void find_range_jump(mpfr_t newFinish, mpfr_t start, mpfr_t finish);
	
	// Keep splitting in range until all segments are flat
	void flatten_range();
	
	// Evaluate function into previously initialised variable
	void eval(mpfr_t res, mpfr_t x);
	
	// Evaluate for fixed-point fraction in [0,0.5) (i.e. excluding implicit bit) over given input and output binade
	void eval_scaled_flat_function(mpfr_t res, mpfr_t x, int eD, int eR);
	
	// Initialise to the range precision, and evaluate
	void init_eval(mpfr_t res, mpfr_t x);
	
	void dump(FILE *dst);
	
	// Get the function scaled to fixed point for input binade eD and output eR
	// The returned node should not be freed (it is cached)
	sollya_node_t get_scaled_flat_function(int eD, int eR);
};

}; // float_approx
}; // random
}; // flopoco

#endif
