#ifndef flopoco_random_float_approx_range_hpp
#define flopoco_random_float_approx_range_hpp

#include <list>
#include <vector>
#include <map>
#include <stdio.h>
#include <mpfr.h>
#include <assert.h>

#include <boost/any.hpp>
#include <boost/smart_ptr.hpp>

#include "FixFunctions/Function.hpp"

namespace flopoco
{
namespace random
{
namespace float_approx
{
	
inline void freeMpfrPtr(void *p)
{
	mpfr_clear(*(mpfr_t*)p);
	free(p);
}

	
typedef boost::shared_ptr<sollya_node> sollya_node_ptr_t;

inline sollya_node_ptr_t make_shared_sollya(sollya_node_t n)
{ return sollya_node_ptr_t(n, free_memory); }

	
/*typedef struct{
	sollya_node_t _x;
	
	sollya_node_t get()
	{ return _x; }
}sollya_node_ptr_t;
	
sollya_node_ptr_t make_shared_sollya(sollya_node_t n)
{ sollya_node_ptr_t x; x._x=n; return x; }
*/

	
void unblockSignals();

// initialise and copy same value (with same precision)
void mpfr_init_copy(mpfr_t dst, mpfr_t src);
	
// Fix to this many fractional bits
void mpfr_fix(mpfr_t x, int bits, mpfr_rnd_t rnd=MPFR_RNDN);
	
struct Range;

struct Segment
{
private:
	Segment();

	boost::shared_ptr<sollya_node> flatFunction;
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
	/* coeffFormats is a list of integers, such that coefficient i will be representable as x/2^coeffFormats[i],
		i.e. it has that many fractional bits :) */
	sollya_node_t fpminimax(const std::vector<int> &coeffFormats, sollya_node_t minimax=NULL);
	
	// Compute minimax with fixed-point coefficients all with the given precision
	sollya_node_t fpminimax(unsigned degree, int coeffFormat, sollya_node_t minimax=NULL);
	
	// The returned node should not be freed (it is cached)
	sollya_node_t get_scaled_flat_function();
	
	bool has_property(std::string name) const;
};

struct Range
{
	const Function &m_function;
	int m_domainWF;
	int m_rangeWF;
	bool m_offsetPolyInputs;
	
	// upper (exclusive) bound for integer value of the domain and range, not including the implicit bit
	mpfr_t m_domainFractionEnd, m_rangeFractionEnd;
	
	typedef std::list<Segment> segment_list_t;
	typedef segment_list_t::iterator segment_it_t;
	segment_list_t m_segments;
	
	typedef std::map<std::pair<int,int>,sollya_node_t> flat_function_cache_t;
	flat_function_cache_t m_flatFunctions;
	
	typedef std::map<std::string,boost::any> property_map_t;
	property_map_t properties;
	
	Range(const Function &f, int domainWF, int rangeWF, mpfr_t domainStart, mpfr_t domainFinish);
	
	// Return a (wE,wF) pair that can represent all values in the domain
	std::pair<int,int> GetFloatTypeEnclosingDomain() const;
	
	// Return a (wE,wF) pair that can represent all values in the range
	std::pair<int,int> GetFloatTypeEnclosingRange() const;
	
	~Range();
	
	/* Split the segment, reducing the range of the original and creating a new one.
		If the original range was  [start,finish], then its range will change to [start,newFinish],
		and the next segment will have [next(newFinish),finish]
	*/
	void split_segment(segment_it_t src, mpfr_t newFinish);
	
	// Splits the segment into two equal parts
	void split_segment(segment_it_t src);
	
	// Split up ranges to ensure that each segment is either monotonically increasing, decreasing, or flat in the range (not nesc. in domain though)
	void make_monotonic_or_range_flat();
	
	// Keep splitting in domain until all segments are flat
	void flatten_domain();
	
	/* Requires: binade(f(start)) != binade(f(finish))
		Ensures: start<=newFinish < finish  and  binade(f(start))==binade(f(newFinish))  and  binade(f(newFinish))!=binade(f(next(newFinish)))*/
	void find_range_jump(mpfr_t newFinish, mpfr_t start, mpfr_t finish);
	
	// Keep splitting in range until all segments are flat
	void flatten_range(bool domainAlreadyFlat=false);
	
	// Evaluate function into previously initialised variable
	void eval(mpfr_t res, mpfr_t x);
	
	// Evaluate for fixed-point fraction in [0,0.5) (i.e. excluding implicit bit) over given input and output binade
	void eval_scaled_flat_function(mpfr_t res, mpfr_t x, int eD, int eR);
	
	// Initialise to the range precision, and evaluate
	void init_eval(mpfr_t res, mpfr_t x);
	
	void dump(FILE *dst);
	
	// Return true if this point is included in the function domain (is in a segment and is representable)
	bool is_in_domain(mpfr_t x);
	
	// Get the function scaled to fixed point for input binade eD and output eR
	// The returned node should not be freed (it is cached)
	sollya_node_t get_scaled_flat_function(int eD, int eR);
	
	segment_it_t find_segment(mpfr_t x);
};

}; // float_approx
}; // random
}; // flopoco

#endif
