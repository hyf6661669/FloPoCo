#include "range.hpp"
#include <stdexcept>
#include <vector>
#include <stdlib.h>

namespace flopoco
{
namespace random
{
namespace float_approx
{

// initialise and copy same value (with same precision)
void mpfr_init_copy(mpfr_t dst, mpfr_t src)
{
	mpfr_init2(dst, mpfr_get_prec(src));
	mpfr_set(dst, src, MPFR_RNDN);
}

// Fix to this many fractional bits
void mpfr_fix(mpfr_t x, int bits)
{
	mpfr_mul_2si(x, x, bits, MPFR_RNDN);
	mpfr_round(x,x);
	mpfr_mul_2si(x, x, -bits, MPFR_RNDN);
}

void Segment::set_domain(mpfr_t _domainStart, mpfr_t _domainFinish)
{
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
	
	mpfr_fprintf(stderr, "Pre eval\n");
	parent->eval(rangeStart, _domainStart);
	parent->eval(rangeFinish, _domainFinish);
	mpfr_fprintf(stderr, "Post eval : (%Rf, %Rf)\n", rangeStart, rangeFinish);
	
	assert(mpfr_sgn(rangeStart) > 0);
	assert(mpfr_lessequal_p(rangeStart, rangeFinish));
	
	isRangeFlat = mpfr_get_exp(rangeStart) == mpfr_get_exp(rangeFinish);
	isDomainFlat = mpfr_get_exp(domainStart) == mpfr_get_exp(domainFinish);
}

Segment::Segment(Range *_parent, mpfr_t _domainStart, mpfr_t _domainFinish)
	: parent(_parent)
{
	mpfr_inits2(parent->m_domainWF, domainStart, domainFinish, domainStartFrac, domainFinishFrac, (mpfr_ptr)0);
	mpfr_inits2(parent->m_rangeWF, rangeStart, rangeFinish, (mpfr_ptr)0);
	
	set_domain(_domainStart, _domainFinish);
	
	fprintf(stderr, "Finish Segment()\n");
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
	
	isRangeFlat=o.isRangeFlat;
	isDomainFlat=o.isDomainFlat;	
	properties=o.properties;
}

Segment &Segment::operator=(const Segment &o)
{
	fprintf(stderr, "Begin copy\n");
	assert(parent==o.parent);
	if(this!=&o){
		mpfr_set(domainStart, o.domainStart, MPFR_RNDN);
		mpfr_set(domainFinish, o.domainFinish, MPFR_RNDN);
		mpfr_set(domainStartFrac, o.domainStartFrac, MPFR_RNDN);
		mpfr_set(domainFinishFrac, o.domainFinishFrac, MPFR_RNDN);
		mpfr_set(rangeStart, o.rangeStart, MPFR_RNDN);
		mpfr_set(rangeFinish, o.rangeFinish, MPFR_RNDN);
		isRangeFlat=o.isRangeFlat;
		isDomainFlat=o.isDomainFlat;
		properties=o.properties;
	}
	fprintf(stderr, "End copy\n");
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
		
		fprintf(stderr, "Created tree : ");
		fprintTree(stderr, tx);
		fprintf(stderr, "\n");
		
		m_flatFunctions[std::make_pair(eD,eR)]=tx;
		return tx;
	}else{
		return m_flatFunctions[std::make_pair(eD,eR)];
	}
}

sollya_node_t Segment::get_scaled_flat_function()
{
	assert(isRangeFlat && isDomainFlat);
	return parent->get_scaled_flat_function(mpfr_get_exp(domainStart), mpfr_get_exp(rangeStart));
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

void Range::eval_scaled_flat_function(mpfr_t res, mpfr_t x, int eD, int eR)
{
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
	
	mpfr_fprintf(stderr, " minimax, dRange=[%Rf,%Rf] -> [%Rf,%Rf]\n", domainStart, domainFinish, fracStart, fracFinish);
	mpfr_fprintf(stderr, "  rRange=[%Rf,%Rf]\n", rangeStart, rangeFinish);
	
	sollya_node_t func=parent->get_scaled_flat_function(mpfr_get_exp(domainStart), mpfr_get_exp(rangeStart));
	sollya_node_t weight=makeConstantDouble(1.0);
	sollya_chain_t monom=makeIntPtrChainFromTo(0, degree);
	
	mp_prec_t prec=getToolPrecision();
	mpfr_t requestedQuality;
	mpfr_init2(requestedQuality, prec);
	mpfr_set_d(requestedQuality, pow(2.0, -parent->m_rangeWF-2), MPFR_RNDN);

	blockSignals();
	
	fprintf(stderr, "Minimax tree : ");
	fprintTree(stderr, func);
	fprintf(stderr, "\n");
	
	if(mpfr_equal_p(domainStart, domainFinish)){
		mpfr_t *coeffs=(mpfr_t*)malloc(sizeof(mpfr_t)*(degree+1));
		for(unsigned i=0;i<=degree;i++){
			mpfr_init2(coeffs[i], getToolPrecision());
			mpfr_set_d(coeffs[i], 0.0, MPFR_RNDN);
		}
		evaluateFaithful(coeffs[0], func, fracStart, prec);
		//mpfr_fix(coeffs[0], parent->m_rangeWF);
		mpfr_fprintf(stderr, "  creating constant : %Rf -> %Rf\n", fracStart, coeffs[0]);
		res=makePolynomial(&coeffs[0], degree);
		for(unsigned i=0;i<=degree;i++){
			mpfr_clear(coeffs[i]);
		}
		free(coeffs);
	}else{	
		res=remez(func, weight, monom, fracStart, fracFinish, &requestedQuality, prec);
	}
		
	// func doesn't need freeing
	freeChain(monom,freeIntPtr);
	free_memory(weight);
	mpfr_clears(fracStart, fracFinish, (mpfr_ptr)0);
	
	return res;
}

sollya_node_t Segment::fpminimax(sollya_chain_t monomials, sollya_chain_t coefficients)
{
	assert(isRangeFlat && isDomainFlat);
	
	sollya_node_t res=NULL;
	
	mpfr_t fracStart, fracFinish;
	mpfr_init2(fracStart, parent->m_domainWF);
	mpfr_init2(fracFinish, parent->m_domainWF);
	
	
	real_to_frac(fracStart, domainStart);
	real_to_frac(fracFinish, domainFinish);
	
	mpfr_fprintf(stderr, " minimax, dRange=[%Rf,%Rf] -> [%Rf,%Rf]\n", domainStart, domainFinish, fracStart, fracFinish);
	mpfr_fprintf(stderr, "  rRange=[%Rf,%Rf]\n", rangeStart, rangeFinish);
	
	sollya_node_t func=parent->get_scaled_flat_function(mpfr_get_exp(domainStart), mpfr_get_exp(rangeStart));
	sollya_node_t weight=makeConstantDouble(1.0);
	sollya_chain_t monom=makeIntPtrChainFromTo(0, degree);
	
	mp_prec_t prec=getToolPrecision();
	mpfr_t requestedQuality;
	mpfr_init2(requestedQuality, prec);
	mpfr_set_d(requestedQuality, pow(2.0, -parent->m_rangeWF-2), MPFR_RNDN);

	blockSignals();
	
	fprintf(stderr, "Minimax tree : ");
	fprintTree(stderr, func);
	fprintf(stderr, "\n");
	
	if(mpfr_equal_p(domainStart, domainFinish)){
		mpfr_t *coeffs=(mpfr_t*)malloc(sizeof(mpfr_t)*(degree+1));
		for(unsigned i=0;i<=degree;i++){
			mpfr_init2(coeffs[i], getToolPrecision());
			mpfr_set_d(coeffs[i], 0.0, MPFR_RNDN);
		}
		evaluateFaithful(coeffs[0], func, fracStart, prec);
		
		// We need to round to whatever the format of coeff zero should be
		mpfr_fix(coeffs[0], parent->m_rangeWF);
		
		mpfr_fprintf(stderr, "  creating constant : %Rf -> %Rf\n", fracStart, coeffs[0]);
		res=makePolynomial(&coeffs[0], degree);
		for(unsigned i=0;i<=degree;i++){
			mpfr_clear(coeffs[i]);
		}
		free(coeffs);
	}else{	
		res=FPminimax(sY, tempChain ,tempChain2, NULL, zero, bi, FIXED, ABSOLUTESYM, tempNode2,NULL);
	}
		
	free_memory(weight);
	mpfr_clears(fracStart, fracFinish, (mpfr_ptr)0);
	
	return res;
}

Range::Range(Function &f, int domainWF, int rangeWF, mpfr_t domainStart, mpfr_t domainFinish)
	: m_function(f)
	, m_domainWF(domainWF)
	, m_rangeWF(rangeWF)
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
	fprintf(stderr, "Pre push_back\n");
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
	mpfr_fprintf(stderr, "split_segment, [%Rf,%Rf] < %Rf\n", src->domainStart, src->domainFinish, newFinish);
	
	// require:   start<=newFinish   and   newFinish < finish
	assert(mpfr_lessequal_p(src->domainStart, newFinish));
	assert(mpfr_less_p(newFinish, src->domainFinish));
	
	mpfr_t newStart;	// starting point for new segment
	mpfr_init_copy(newStart, newFinish);
	mpfr_nextabove(newStart);
	
	Segment left(this, src->domainStart, newFinish), right(this, newStart, src->domainFinish);
	fprintf(stderr, "  left=");
	left.dump(stderr);
	fprintf(stderr, "\n  right=");
	right.dump(stderr);
	fprintf(stderr, "\n");
	
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
	assert(mpfr_less_p(start, finish));
	
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
	
	assert(lowE<highE);
	
	while(true){
		mpfr_fprintf(stderr, "  low=%Rf, high=%Rf, lowE=%d, highE=%d\n", low, high, lowE, highE);
		
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
	
	mpfr_printf("  newFinish=%Rf\n", newFinish);
	
	mpfr_clears(low, mid, high, val, (mpfr_ptr)0);
}

void Range::flatten_range()
{
	mpfr_t newFinish;
	mpfr_init2(newFinish, m_domainWF);
	
	segment_it_t curr=m_segments.begin();
	while(curr!=m_segments.end()){
		if(!curr->isRangeFlat){
			// So currently we have [start,finish] with binade(f(start)) < binade(f(finish))
			// We need to find newFinish where  start<=newFinish<finish and binade(f(newFinish)) < binade(f(next(newFinish)))
		
			find_range_jump(newFinish, curr->domainStart, curr->domainFinish);
		
			mpfr_fprintf(stderr, "  newFinish=%Rf\n", newFinish);
			split_segment(curr, newFinish);
			
			assert(curr->isRangeFlat);
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

}; // float_approx
}; // random
}; // flopoco

