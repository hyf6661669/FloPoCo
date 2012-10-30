#include "range.hpp"
#include <stdexcept>
#include <vector>
#include <stdlib.h>
#include <boost/utility.hpp>

#include "Operator.hpp"

#include "mpreal.h"

namespace flopoco
{
namespace random
{
namespace float_approx
{
	
// Try to undo sollya's "helpful" stuff
void unblockSignals()
{
  sigset_t mask;

  sigemptyset(&mask);
  sigaddset(&mask,SIGINT);
  sigaddset(&mask,SIGSEGV);
  sigaddset(&mask,SIGBUS);
  sigaddset(&mask,SIGFPE);
  sigaddset(&mask,SIGPIPE);
  sigprocmask(SIG_UNBLOCK, &mask, NULL);
  signal(SIGINT,SIG_DFL);
  signal(SIGSEGV,SIG_DFL);
  signal(SIGBUS,SIG_DFL);
  signal(SIGFPE,SIG_DFL);
  signal(SIGPIPE,SIG_DFL);
}

// initialise and copy same value (with same precision)
void mpfr_init_copy(mpfr_t dst, mpfr_t src)
{
	mpfr_init2(dst, mpfr_get_prec(src));
	mpfr_set(dst, src, MPFR_RNDN);
}

// Fix to this many fractional bits
void mpfr_fix(mpfr_t x, int bits, mpfr_rnd_t rnd)
{
	mpfr_mul_2si(x, x, bits, MPFR_RNDN);
	mpfr_rint(x,x, rnd);
	mpfr_mul_2si(x, x, -bits, MPFR_RNDN);
}

void Segment::set_domain(mpfr_t _domainStart, mpfr_t _domainFinish)
{
	flatFunction.reset();
	
	assert(mpfr_sgn(_domainStart) > 0);
	assert(mpfr_lessequal_p (_domainStart, _domainFinish));
	
	mpfr_set(domainStart, _domainStart, MPFR_RNDN);
	mpfr_set(domainFinish, _domainFinish, MPFR_RNDN);
	
	mpfr_mul_2si(domainStartFrac, domainStart, -mpfr_get_exp(domainStart), MPFR_RNDN);
	mpfr_sub_d(domainStartFrac, domainStartFrac, 0.5, MPFR_RNDN);
	mpfr_fix(domainStartFrac, parent->m_domainWF);
	
	mpfr_mul_2si(domainFinishFrac, domainFinish, -mpfr_get_exp(domainFinish), MPFR_RNDN);
	mpfr_sub_d(domainFinishFrac, domainFinishFrac, 0.5, MPFR_RNDN);
	mpfr_fix(domainFinishFrac, parent->m_domainWF);
	
	parent->eval(rangeStart, _domainStart);
	parent->eval(rangeFinish, _domainFinish);
	
	assert(mpfr_sgn(rangeStart) > 0);
	
	isRangeFlat = mpfr_get_exp(rangeStart) == mpfr_get_exp(rangeFinish);
	isDomainFlat = mpfr_get_exp(domainStart) == mpfr_get_exp(domainFinish);
}

Segment::Segment(Range *_parent, mpfr_t _domainStart, mpfr_t _domainFinish)
	: parent(_parent)
{
	mpfr_inits2(parent->m_domainWF, domainStart, domainFinish, domainStartFrac, domainFinishFrac, (mpfr_ptr)0);
	mpfr_inits2(parent->m_rangeWF, rangeStart, rangeFinish, (mpfr_ptr)0);
	
	set_domain(_domainStart, _domainFinish);
}

Segment::Segment(const Segment &o)
	: parent(o.parent)
{
	mpfr_inits2(parent->m_domainWF, domainStart, domainFinish, domainStartFrac, domainFinishFrac, (mpfr_ptr)0);
	mpfr_inits2(parent->m_rangeWF, rangeStart, rangeFinish, (mpfr_ptr)0);
	
	mpfr_set(domainStart, o.domainStart, MPFR_RNDN);
	mpfr_set(domainFinish, o.domainFinish, MPFR_RNDN);
	mpfr_set(domainStartFrac, o.domainStartFrac, MPFR_RNDN);
	mpfr_set(domainFinishFrac, o.domainFinishFrac, MPFR_RNDN);
	mpfr_set(rangeStart, o.rangeStart, MPFR_RNDN);
	mpfr_set(rangeFinish, o.rangeFinish, MPFR_RNDN);
	
	flatFunction=o.flatFunction;
	isRangeFlat=o.isRangeFlat;
	isDomainFlat=o.isDomainFlat;	
	properties=o.properties;
}

Segment &Segment::operator=(const Segment &o)
{
	assert(parent==o.parent);
	if(this!=&o){
		mpfr_set(domainStart, o.domainStart, MPFR_RNDN);
		mpfr_set(domainFinish, o.domainFinish, MPFR_RNDN);
		mpfr_set(domainStartFrac, o.domainStartFrac, MPFR_RNDN);
		mpfr_set(domainFinishFrac, o.domainFinishFrac, MPFR_RNDN);
		mpfr_set(rangeStart, o.rangeStart, MPFR_RNDN);
		mpfr_set(rangeFinish, o.rangeFinish, MPFR_RNDN);
		flatFunction=o.flatFunction;
		isRangeFlat=o.isRangeFlat;
		isDomainFlat=o.isDomainFlat;
		properties=o.properties;
	}
	return *this;
}

Segment::~Segment()
{
	mpfr_clears(domainStart, domainFinish,rangeStart, rangeFinish, domainStartFrac, domainFinishFrac, (mpfr_ptr)0);
}

void Segment::dump(FILE *dst) const
{
	mpfr_fprintf(dst, "%s[%Rg,%Rg] -> %s[%Rg,%Rg]",
		isDomainFlat?"flat":"split", domainStart, domainFinish,
		isRangeFlat?"flat":"split", rangeStart, rangeFinish
	);
}

bool Segment::has_property(std::string name) const
{
	return properties.find(name)!=properties.end();
}



/* This function is flat, with eD=binade(domainStart), eR=binade(f(domainStart))
	That means that 	2^(eD-1) <= domainStart <= domainFinish < 2^eD
	and 2^(eR-1) <= f(domainStart) <= f(domainFinish) < 2^eR
	we want to convert to fixed-point, so we have inputs in [0,1) with domainWF-1
	fractional bits, and outputs in [0,1) with rangeWF-1 fractional bits.
	
	To move a real x to fixed-point we have:
	0.5*2^eD <= x < 1.0*2^eD
	0.5 <= x * 2^-eD < 1.0
	0 <= x * 2^-eD - 0.5 < 0.5
	0 <= 2*(x * 2^-eD - 0.5) < 1.0
	0 <= x*2^(-eD+1) -1.0 < 1.0

	To convert a fixed-point r to real we have
	0 <= r < 1.0
	0 <= (r+1.0) * 2^(eR-1) < 1.0
*/
sollya_node_t Range::get_scaled_flat_function(int eD, int eR)
{	
	if(m_flatFunctions.find(std::make_pair(eD,eR))==m_flatFunctions.end()){
		// I'm sure this must be leaking sollya objects. The caching is a hedge against
		// that, rather than about efficiency
		
		sollya_node_t f=m_function.getSollyaNode();
		
		sollya_node_t tx=makeAdd(makeVariable(), makeConstantDouble(0.5));
		tx=makeMul(tx, makeConstantDouble(pow(2.0, eD)));
		
		tx=substitute(f, tx);
		
		tx=makeMul(tx, makeConstantDouble(pow(2.0, -eR)));
		tx=makeSub(tx, makeConstantDouble(0.5));
		
		//fprintTree(stderr, tx);
		
		m_flatFunctions[std::make_pair(eD,eR)]=tx;
		return tx;
	}else{
		return m_flatFunctions[std::make_pair(eD,eR)];
	}
}

sollya_node_t Segment::get_scaled_flat_function()
{
	assert(isRangeFlat && isDomainFlat);
	if(parent->m_offsetPolyInputs){
		if(!flatFunction){
			// We want a function where the input starts at 0, so we add on domainFracStart beforehand
			sollya_node_t offset=makeConstant(domainStartFrac);
			offset=makeAdd(makeVariable(), offset);
			
			sollya_node_t tree=substitute(parent->get_scaled_flat_function(mpfr_get_exp(domainStart), mpfr_get_exp(rangeStart)), offset);
			flatFunction.reset(tree, free_memory);
		}
		return flatFunction.get();
	}else{
		return parent->get_scaled_flat_function(mpfr_get_exp(domainStart), mpfr_get_exp(rangeStart));
	}
}

static void real_to_frac(mpfr_t frac, mpfr_t real)
{
	// This paranoia is because I seemed to be getting non lossless transforms at some point
	
	mpfr_t tmp;
	mpfr_init_copy(tmp, frac);
	
	// Go forwards
	mpfr_set(tmp, real, MPFR_RNDN);
	mpfr_mul_2si(tmp, tmp, -mpfr_get_exp(real), MPFR_RNDN);
	mpfr_sub_d(tmp, tmp, 0.5, MPFR_RNDN);
	
	// Go backwards
	mpfr_add_d(tmp, tmp, 0.5, MPFR_RNDN);
	mpfr_mul_2si(tmp, tmp, mpfr_get_exp(real), MPFR_RNDN);
	assert(mpfr_equal_p(tmp, real));
	
	// Go forwards again
	mpfr_mul_2si(tmp, tmp, -mpfr_get_exp(real), MPFR_RNDN);
	mpfr_sub_d(tmp, tmp, 0.5, MPFR_RNDN);
	
	mpfr_set(frac, tmp, MPFR_RNDN);
	mpfr_clear(tmp);
}

static void frac_to_real(mpfr_t real, mpfr_t frac, int e)
{
	// This paranoia is because I seemed to be getting non lossless transforms at some point
	
	mpfr_t tmp;
	mpfr_init_copy(tmp, real);
	
	// Go forwards
	mpfr_set(tmp, frac, MPFR_RNDN);
	mpfr_add_d(tmp, tmp, 0.5, MPFR_RNDN);
	mpfr_mul_2si(tmp, tmp, e, MPFR_RNDN);
	
	// Go backwards
	mpfr_mul_2si(tmp, tmp, -e, MPFR_RNDN);
	mpfr_sub_d(tmp, tmp, 0.5, MPFR_RNDN);
	assert(mpfr_equal_p(tmp, frac));
	
	// Go forwards again
	mpfr_add_d(tmp, tmp, 0.5, MPFR_RNDN);
	mpfr_mul_2si(tmp, tmp, e, MPFR_RNDN);
	
	mpfr_set(real, tmp, MPFR_RNDN);
	mpfr_clear(tmp);
}

// Return a (wE,wF) pair that can represent all values in the domain
std::pair<int,int> Range::GetFloatTypeEnclosingDomain() const
{
	int minE=mpfr_get_exp(m_segments.begin()->domainStart);
	int maxE=mpfr_get_exp(boost::prior(m_segments.begin())->domainFinish);
	int wE=1;
	while( (minE<(1<<(wE-1))) || ((1<<(wE-1)) < maxE)){
		wE++;
	}
	return std::make_pair(wE, m_domainWF-1);
}

/*
// Return a (wE,wF) pair that can represent all values in the range
std::pair<int,int> Range::GetFloatTypeEnclosingRange() const
{
	
}
*/

void Range::eval_scaled_flat_function(mpfr_t res, mpfr_t x, int eD, int eR)
{
	throw std::string("Here");	// Is this function used?
	
	mpfr_t tx;
	mpfr_init_copy(tx, x);
	
	assert(mpfr_sgn(tx) >= 0);
	assert(mpfr_less_p(tx,m_domainFractionEnd));
	
	frac_to_real(tx, tx, eD);
	
	assert(mpfr_get_exp(tx)==eD);
	eval(res, tx);
	assert(mpfr_get_exp(res)==eR);
	
	real_to_frac(res, res);
	
	assert(mpfr_sgn(res) >= 0);
	assert(mpfr_less_p(res, m_rangeFractionEnd));
	
	mpfr_clear(tx);
}

sollya_node_t Segment::minimax(unsigned degree)
{
	assert(isRangeFlat && isDomainFlat);
	
	sollya_node_t res=NULL;
	
	mpfr_t fracStart, fracFinish;
	mpfr_init2(fracStart, parent->m_domainWF);
	mpfr_init2(fracFinish, parent->m_domainWF);
	
	
	real_to_frac(fracStart, domainStart);
	real_to_frac(fracFinish, domainFinish);
	
	//mpfr_fprintf(stderr, " minimax, dRange=[%Rf,%Rf] -> [%Rf,%Rf]\n", domainStart, domainFinish, rangeStart, rangeFinish);
	
	sollya_node_t func=get_scaled_flat_function();
	sollya_node_t weight=makeConstantDouble(1.0);
	sollya_chain_t monom=makeIntPtrChainFromTo(0, degree);
	
	mp_prec_t prec=getToolPrecision();
	mpfr_t requestedQuality;
	mpfr_init2(requestedQuality, prec);
	mpfr_set_d(requestedQuality, pow(2.0, -parent->m_rangeWF-10), MPFR_RNDN);
	
	mpfr::mpreal zero;

	if(mpfr_equal_p(domainStart, domainFinish)){
		mpfr_t *coeffs=(mpfr_t*)malloc(sizeof(mpfr_t)*(degree+1));
		for(unsigned i=0;i<=degree;i++){
			mpfr_init2(coeffs[i], getToolPrecision());
			mpfr_set_d(coeffs[i], 0.0, MPFR_RNDN);
		}
		if(!parent->m_offsetPolyInputs){
			evaluateFaithful(coeffs[0], func, fracStart, prec);
		}else{
			evaluateFaithful(coeffs[0], func, zero.mpfr_ptr(), prec);
		}
		//mpfr_fix(coeffs[0], parent->m_rangeWF);
		res=makePolynomial(&coeffs[0], degree);
		for(unsigned i=0;i<=degree;i++){
			mpfr_clear(coeffs[i]);
		}
		free(coeffs);
	}else{	
		if(!parent->m_offsetPolyInputs){
			res=remez(func, weight, monom, fracStart, fracFinish, &requestedQuality, prec);
		}else{
			mpfr::mpreal length(fracFinish);
			length=length-fracStart;
			res=remez(func, weight, monom, zero.mpfr_ptr(), length.mpfr_ptr(), &requestedQuality, prec);
		}
	}
		
	// func doesn't need freeing
	freeChain(monom,freeIntPtr);
	free_memory(weight);
	mpfr_clears(fracStart, fracFinish, (mpfr_ptr)0);
	
	return res;
}

sollya_chain_t makeIntPtrChain(const std::vector<int> &values)
{
	sollya_chain_t res=NULL;
	for(int j=values.size()-1;j>=0;--j){
		int *elem=(int*)malloc(sizeof(int));
		*elem=values[j];
		res=addElement(res, elem);
	}
	return res;
}

sollya_node_t Segment::fpminimax(const std::vector<int> &coeffWidths, sollya_node_t minimax)
{
	unsigned degree=coeffWidths.size()-1;
	
	assert(isRangeFlat && isDomainFlat);
	if(minimax){
		assert(isPolynomial(minimax));
		assert(getDegree(minimax)==(int)degree);
	}
	
	sollya_node_t res=NULL;
	
	//mpfr_fprintf(stderr, " fpminimax, dRange=[%Rf,%Rf] -> [%Rf,%Rf]\n", domainStart, domainFinish, rangeStart, rangeFinish);
	
	sollya_node_t func=get_scaled_flat_function();
	
	mp_prec_t prec=getToolPrecision();
	mpfr::mpreal zero;
	
	if(mpfr_equal_p(domainStart, domainFinish)){
		mpfr_t *coeffs=(mpfr_t*)malloc(sizeof(mpfr_t)*(degree+1));
		for(unsigned i=0;i<=degree;i++){
			mpfr_init2(coeffs[i], getToolPrecision());
			mpfr_set_d(coeffs[i], 0.0, MPFR_RNDN);
		}
		if(!parent->m_offsetPolyInputs){
			evaluateFaithful(coeffs[0], func, domainStartFrac, prec);
		}else{
			evaluateFaithful(coeffs[0], func, zero.mpfr_ptr(), prec);
		}
		
		// We need to round to whatever the format of coeff zero should be
		mpfr_fix(coeffs[0], coeffWidths[0]);
		
		res=makePolynomial(&coeffs[0], degree);
		for(unsigned i=0;i<=degree;i++){
			mpfr_clear(coeffs[i]);
		}
		free(coeffs);
	}else{
		sollya_chain_t monom=makeIntPtrChainFromTo(0, degree);
		sollya_chain_t formats=makeIntPtrChain(coeffWidths);
		sollya_node_t constantPart=makeConstantDouble(0.0);
		
		if(!parent->m_offsetPolyInputs){
			res=FPminimax(func, monom, formats, /*points*/NULL, domainStartFrac, domainFinishFrac, FIXED, ABSOLUTESYM, constantPart,minimax);
		}else{
			mpfr::mpreal length(domainFinishFrac);
			length=length-domainStartFrac;
			res=FPminimax(func, monom, formats, /*points*/NULL, zero.mpfr_ptr(), length.mpfr_ptr(), FIXED, ABSOLUTESYM, constantPart,minimax);
		}
		
		free_memory(constantPart);
		freeChain(monom,freeIntPtr);
		freeChain(formats,freeIntPtr);
	}
	
	return res;
}

sollya_node_t Segment::fpminimax(unsigned degree, int coeffWidths, sollya_node_t minimax)
{
	std::vector<int> coeffWidthsVec(degree+1, coeffWidths);
	
	return fpminimax(coeffWidthsVec, minimax);
}

Range::Range(const Function &f, int domainWF, int rangeWF, mpfr_t domainStart, mpfr_t domainFinish)
	: m_function(f)
	, m_domainWF(domainWF)
	, m_rangeWF(rangeWF)
	, m_offsetPolyInputs(false)
{
	if(mpfr_sgn(domainStart)<0)
		throw std::invalid_argument("Domain must be strictly positive.");
	if(mpfr_greater_p(domainStart, domainFinish))
		throw std::invalid_argument("Domain must be ordered.");
	
	mpfr_init2(m_domainFractionEnd, domainWF);
	mpfr_init2(m_rangeFractionEnd, rangeWF);
	
	mpfr_set_d(m_domainFractionEnd, 0.5, MPFR_RNDN);
	mpfr_set_d(m_rangeFractionEnd, 0.5, MPFR_RNDN);
	
	Segment base(this, domainStart, domainFinish);
	m_segments.push_back(base);
}

Range::~Range()
{
	mpfr_clears(m_domainFractionEnd, m_rangeFractionEnd, (mpfr_ptr)0);
};

/* Split the segment, reducing the range of the original and creating a new one.
	If the original range was  [start,finish], then its range will change to [start,newFinish],
	and the next segment will have [next(newFinish),finish]
*/
void Range::split_segment(segment_it_t src, mpfr_t newFinish)
{
	if(::flopoco::verbose>=DEBUG){
		mpfr_fprintf(stderr, "split_segment, [%Rf,%Rf] < %Rf,  ", src->domainStart, src->domainFinish, newFinish);
	}
	
	// require:   start<=newFinish   and   newFinish < finish
	assert(mpfr_lessequal_p(src->domainStart, newFinish));
	assert(mpfr_less_p(newFinish, src->domainFinish));
	
	mpfr_t newStart;	// starting point for new segment
	mpfr_init_copy(newStart, newFinish);
	mpfr_nextabove(newStart);
	
	Segment left(this, src->domainStart, newFinish), right(this, newStart, src->domainFinish);
	if(::flopoco::verbose>=DEBUG){
		fprintf(stderr, "  left=");
		left.dump(stderr);
		fprintf(stderr, ", right=");
		right.dump(stderr);
		fprintf(stderr, "\n");
	}
	
	*src=left;
	segment_it_t tmp=src;
	++tmp;
	m_segments.insert(tmp, right);
}

void Range::split_segment(segment_it_t src)
{
	// Segment can't be a point
	assert(mpfr_less_p(src->domainStart,src->domainFinish));
	
	mpfr_t mid;
	mpfr_init2(mid, m_domainWF);
	mpfr_add(mid, src->domainStart, src->domainFinish, MPFR_RNDD); // note that we round down here
	mpfr_mul_2si(mid, mid, -1, MPFR_RNDN);
	
	split_segment(src, mid);

	mpfr_clear(mid);
}

void Range::make_monotonic_or_range_flat()
{
	sollya_node_ptr_t deriv=make_shared_sollya(differentiate(m_function.getSollyaNode()));
	
	fprintf(stderr, "Deriv = ");
	fprintTree(stderr, deriv.get());
	fprintf(stderr, "\n");
	
	sollya_range_t rr={ &(m_segments.begin()->domainStart), &(boost::prior(m_segments.end())->domainFinish) };
	sollya_chain_t zeros=fpFindZerosFunction(deriv.get(), rr, getToolPrecision());
	
	mpfr_t extremaLow, extremaHigh, extremaLowRange, extremaHighRange;
	mpfr_inits2(m_domainWF, extremaLow, extremaHigh, (mpfr_ptr)0);
	mpfr_inits2(m_rangeWF, extremaLowRange, extremaHighRange, (mpfr_ptr)0);
	
	sollya_chain_t tmp=zeros;
	while (tmp != NULL) {
		mpfr_t &extrema = *(mpfr_t*)(tmp->value);	// some maximum or minimum
		
		mpfr_fprintf(stderr, "  extrema at %Rg\n", extrema);
		
		// The extrema might be at a representable number...
		if(0==mpfr_set(extremaLow, extrema, MPFR_RNDD)){
			// in which case it must be the first or last element of a segment
			segment_it_t seg=Range::find_segment(extremaLow);
			
			if(mpfr_equal_p(seg->domainStart,extremaLow) || mpfr_equal_p(seg->domainFinish,extremaHigh)){
				// This is fine, one segment already finishes on this boundary, so the function must
				// be curving up (down) from the minima (maxima)
				fprintf(stderr, "    extrema already on function boundary.\n");
			}else{
				// otherwise the maxima is located within this segment.
				if((mpfr_get_exp(seg->domainStart) == mpfr_get_exp(extremaLow)) &&(mpfr_get_exp(seg->domainFinish) == mpfr_get_exp(extremaLow))){
					// We don't need to do anything, as this segment is already flat in the range
					// TODO : There is a slight risk that this might result in functions that are very difficult
					// to approximate, and it would be better to split at the extrema. However, this should
					// be handled later on by the adaptive splitting if it occurs.
					fprintf(stderr, "    extrema already in flat range segment.\n");
				}else{
					// Split the range at the (exact) extrema
					// Each segment should now be curving up or down away from the boundary
					split_segment(seg, extremaLow);
				}
			}			
		}else{
			mpfr_set(extremaHigh, extremaLow, MPFR_RNDN);
			mpfr_nextabove(extremaHigh);
			eval(extremaLowRange, extremaLow);
			eval(extremaHighRange, extremaHigh);
			mpfr_fprintf(stderr, "    extremaLow=%Rg, extremaHigh=%Rg\n", extremaLow, extremaHigh);
			fprintf(stderr,"    lowInDomain=%d, highInDomain=%d\n", is_in_domain(extremaLow)?1:0, is_in_domain(extremaHigh)?1:0);
			
			// The function has an extrema at extremaLow < extrema < extremaHigh
			if((!is_in_domain(extremaLow)) || (!is_in_domain(extremaHigh))){
				// We don't have to worry, as one or both of the extrema turns out
				// not to be representable. If one is representable, then they other
				// must be curving up or down quite happily
				fprintf(stderr, "    one side of extrema not in domain.\n");
			}else{
				segment_it_t seg=Range::find_segment(extremaLow);
				mpfr_fprintf(stderr, "    [%Rg,%Rg]\n", seg->domainStart, seg->domainFinish);
				
				if(mpfr_equal_p(seg->domainFinish, extremaLow)){
					// The lower bound is at the end of a segment, we don't need to
					// do anything
					fprintf(stderr, "    extremaLow already at finish of segment.\n");
				}else if((mpfr_get_exp(seg->rangeStart)==mpfr_get_exp(extremaLowRange)) && (mpfr_get_exp(seg->rangeFinish)==mpfr_get_exp(extremaLowRange))){
					// The segment is flat, we don't need to do anything
					fprintf(stderr, "    extremaLow already in a flat segment.\n");
				}else{
					// We are splitting
					fprintf(stderr, "    extremaLow forcing split.\n");
					split_segment(seg, extremaLow);
				}
				
				// This is the same again. I'm pretty sure it won't do anything.
				seg=Range::find_segment(extremaHigh);
				if(mpfr_equal_p(seg->domainFinish, extremaHigh)){
					fprintf(stderr, "    extremaHigh already at finish of segment.\n");
				}else if((mpfr_get_exp(seg->rangeStart)==mpfr_get_exp(extremaHighRange)) && (mpfr_get_exp(seg->rangeFinish)==mpfr_get_exp(extremaHighRange))){
					fprintf(stderr, "    extremaHigh already in a flat segment.\n");
				}else{
					fprintf(stderr, "    extremaHigh forcing split.\n");
					split_segment(seg, extremaHigh);
				}
			}
		}
			
		
		tmp=tmp->next;
	}
	
	mpfr_clears(extremaLow, extremaHigh, extremaLowRange, extremaHighRange, (mpfr_ptr)0);
	freeChain(zeros, freeMpfrPtr);
}

// Keep splitting range until all segments are flat
// This could be done more efficiently at construction time
void Range::flatten_domain()
{
	mpfr_t newFinish;
	mpfr_init2(newFinish, m_domainWF);
	
	segment_it_t curr=m_segments.begin();
	while(curr!=m_segments.end()){
		if(!curr->isDomainFlat){
			// So currently we have [start,finish] with binade(start) < binade(finish)
			// We need to find newFinish where  start<=newFinish and binade(start)==binade(newFinish)
			int eStart=mpfr_get_exp(curr->domainStart);
			// To get the endpoint, we can just do prev(2^(e+1))
			
			mpfr_set_si_2exp(newFinish, 1, eStart, MPFR_RNDN);
			mpfr_nextbelow(newFinish);
			
			split_segment(curr, newFinish);
			
			assert(curr->isDomainFlat);
		}
		++curr;
	}
	
	mpfr_clear(newFinish);
}

void Range::find_range_jump(mpfr_t newFinish, mpfr_t start, mpfr_t finish)
{
	mpfr_t low, mid, high;
	mpfr_t val;
	int lowE, highE, midE;
	
	mpfr_init_copy(low, start);
	mpfr_init_copy(high, finish);
	mpfr_init2(mid, m_domainWF);
	mpfr_init2(val, m_rangeWF);
	
	eval(val, low);
	lowE=mpfr_get_exp(val);
	eval(val, high);
	highE=mpfr_get_exp(val);
	
	assert(lowE!=highE);
	
	while(true){
		//mpfr_fprintf(stderr, "  low=%Rf, high=%Rf, lowE=%d, highE=%d\n", low, high, lowE, highE);
		
		int finish=false;
		mpfr_nextabove(low);
		finish = mpfr_equal_p(low, high) != 0;
		mpfr_nextbelow(low);
		
		if(finish)
			break;
		
		mpfr_add(mid, low, high, MPFR_RNDN);
		mpfr_div_2ui (mid, mid, 1, MPFR_RNDN);
		
		eval(val, mid);
		midE=mpfr_get_exp(val);
		
		if(lowE==midE){
			mpfr_set(low, mid, MPFR_RNDN);
			lowE=midE;
		}else{
			mpfr_set(high, mid, MPFR_RNDN);
			highE=midE;
		}
	}
	
	mpfr_set(newFinish, low, MPFR_RNDN);
	
	//mpfr_printf("  newFinish=%Rf\n", newFinish);
	
	mpfr_clears(low, mid, high, val, (mpfr_ptr)0);
}

void Range::flatten_range(bool domainAlreadyFlat)
{
	mpfr_t newFinish;
	mpfr_init2(newFinish, m_domainWF);
	
	segment_it_t curr=m_segments.begin();
	while(curr!=m_segments.end()){
		if(!curr->isRangeFlat){
			// So currently we have [start,finish] with binade(f(start)) != binade(f(finish))
			// We need to find newFinish where  start<=newFinish<finish and binade(f(newFinish)) != binade(f(next(newFinish)))
		
			find_range_jump(newFinish, curr->domainStart, curr->domainFinish);
		
			//mpfr_fprintf(stderr, "  newFinish=%Rf\n", newFinish);
			split_segment(curr, newFinish);
			
			assert(curr->isRangeFlat);
			assert(domainAlreadyFlat ? curr->isRangeFlat : true);
		}
		if(curr->isRangeFlat)
			++curr;
	}
	
	mpfr_clear(newFinish);
}

void Range::eval(mpfr_t res, mpfr_t x)
{
	m_function.eval(res, x);
}

void Range::init_eval(mpfr_t res, mpfr_t x)
{
	mpfr_init2(res, m_rangeWF);
	m_function.eval(res, x);
}

void Range::dump(FILE *dst)
{
	segment_it_t curr=m_segments.begin();
	while(curr!=m_segments.end()){
		curr->dump(dst);
		fprintf(dst, "\n");
		++curr;
	}
	fprintf(dst, "Total segments = %u\n", m_segments.size());
}

// TODO : This is O(n)... fix it
bool Range::is_in_domain(mpfr_t x)
{
	bool res=false;
	
	mpfr_t tmp;
	mpfr_init2(tmp, m_domainWF);
	if(0==mpfr_set(tmp, x, MPFR_RNDN)){
		segment_it_t curr=m_segments.begin();
		while(curr!=m_segments.end()){
			if(mpfr_lessequal_p(curr->domainStart, tmp) && mpfr_lessequal_p(tmp,curr->domainFinish)){
				res=true;
				break;
			}				
			++curr;
		}
	}
	mpfr_clear(tmp);
	return res;
}

// TODO : This is O(n)... fix it
Range::segment_it_t Range::find_segment(mpfr_t x)
{
	segment_it_t curr=m_segments.begin();
	while(curr!=m_segments.end()){
		if(mpfr_lessequal_p(curr->domainStart, x) && mpfr_lessequal_p(x,curr->domainFinish))
			return curr;
		++curr;
	}
	
	throw std::invalid_argument("find_segment - value out of range.");
}

}; // float_approx
}; // random
}; // flopoco

