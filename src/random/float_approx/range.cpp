#include "range.hpp"
#include <stdexcept>
#include <vector>
#include <stdlib.h>
#include <signal.h>
#include <boost/utility.hpp>

#include "dirty_fixed_minimax.hpp"

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
#ifndef WIN32
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
#endif
	    }

	    // 2^binade(x) <= x < 2^(binade(x+1))
	    std::pair<NumberClass,int> binade(mpfr_t x)
	    {
		if(mpfr_inf_p(x)){
		    return std::make_pair(mpfr_sgn(x)<0 ? NegInf : PosInf,0);
		}else if(mpfr_zero_p(x)){
		    return std::make_pair(mpfr_sgn(x)<0 ? NegZero : PosZero,0);
		}else if(mpfr_nan_p(x)){
		    return std::make_pair(NaN,0);
		}

		int e=mpfr_get_exp(x)-1;
		if(mpfr_sgn(x)>0){
		    assert(mpfr_cmp_si_2exp(x,1,e)>=0);
		    assert(mpfr_cmp_si_2exp(x,1,e+1)<0);
		}else{
		    assert(mpfr_cmp_si_2exp(x,-1,e+1)>0);
		    assert(mpfr_cmp_si_2exp(x,-1,e)<=0);
		}
		return std::make_pair(mpfr_sgn(x)>0 ? PosNormal : NegNormal, e);
	    }

	    /*! This returns the positive fractional part of the
	      number, not including the implicit bit. */
	    static void real_to_frac(mpfr_t frac, mpfr_t real)
	    {
		// This paranoia is because I seemed to be getting non lossless transforms at some point
	
		if(mpfr_zero_p(real)){
		    mpfr_set_d(frac, 0.0, MPFR_RNDN);
		    return;
		}

		mpfr_t res, tmp;
		mpfr_init2(tmp, mpfr_get_prec(real));
		mpfr_init2(res, mpfr_get_prec(real));
	
		// Go forwards
		mpfr_set(tmp, real, MPFR_RNDN);
		mpfr_abs(tmp, tmp, MPFR_RNDN);
		mpfr_mul_2si(tmp, tmp, -mpfr_get_exp(real)+1, MPFR_RNDN);
		mpfr_sub_d(res, tmp, 1.0, MPFR_RNDN);
	
		// Go backwards
		mpfr_add_d(tmp, res, 1.0, MPFR_RNDN);
		mpfr_mul_2si(tmp, tmp, mpfr_get_exp(real)-1, MPFR_RNDN);
		if(mpfr_sgn(real)<0)
		    mpfr_neg(tmp,tmp,MPFR_RNDN);
		if(!mpfr_equal_p(tmp, real)){
		    std::vector<char> tmp(8192,0);
		    mpfr_snprintf(&tmp[0], tmp.size()-1, "real_to_frac(%Re)=%Re - number does not round-trip.", real, frac);
		    throw std::string(&tmp[0]);
		}
	
		if(0!=mpfr_set(frac, res, MPFR_RNDN)){
		    mpfr_fprintf(stderr, " Assiging %Re (prec=%d) to (prec=%d)\n",
				 res, mpfr_get_prec(res), mpfr_get_prec(frac));
		    throw std::string("real_to_frac - Target does not have enough bits.");
		}
		mpfr_clear(tmp);
		mpfr_clear(res);
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

	    /*! This is now extended across the entire
	     * 	domain of floats, including infinities and nans
	     * 	The equivalent flat binades are:
	     * 	-inf
	     * 	-normal
	     * 	-zero
	     * 	+zero
	     * 	+normal
	     * 	+inf
	     * 	nan (all types)
	     */
	    bool isFlat(mpfr_t start, mpfr_t finish)
	    {
		return binade(start)==binade(finish);
	    }


	    std::ostream &operator<<(std::ostream &dst, const binade_t &b)
	    {
		switch(b.first){
		case NaN:
		    dst<<"NaN";break;
		case NegInf:
		    dst<<"NegInf"; break;
		case NegNormal:
		    dst<<"NegNormal<e="<<b.second<<">"; break;
		case NegZero:
		    dst<<"NegZero"; break;
		case PosZero:
		    dst<<"PosZero"; break;
		case PosNormal:
		    dst<<"PosNormal<u="<<b.second<<">"; break;
		case PosInf:
		    dst<<"PosInf"; break;
		default:
		    throw std::invalid_argument("Invalid binade.");
		}
		return dst;
	    }


	    void Segment::set_domain(mpfr_t _domainStart, mpfr_t _domainFinish)
	    {

		assert(mpfr_lessequal_p (_domainStart, _domainFinish));
	
		mpfr_set(domainStart, _domainStart, MPFR_RNDN);
		mpfr_set(domainFinish, _domainFinish, MPFR_RNDN);
	
		if(mpfr_regular_p(domainStart)){
		    if(mpfr_sgn(domainStart)>0){
			real_to_frac(domainStartFrac, domainStart);
			real_to_frac(domainFinishFrac, domainFinish);
		    }else{
			real_to_frac(domainStartFrac, domainFinish);
			real_to_frac(domainFinishFrac, domainStart);
	    
		    }
		}else{
		    mpfr_set_d(domainStartFrac, 0.0, MPFR_RNDN);
		    mpfr_set_d(domainFinishFrac, 0.0, MPFR_RNDN);
		}

		parent->eval_clamp(rangeStart, _domainStart);
		parent->eval_clamp(rangeFinish, _domainFinish);

		if(mpfr_regular_p(rangeStart)){
		    real_to_frac(rangeStartFrac, rangeStart);
		}else{
		    mpfr_set_d(rangeStartFrac, 0.0, MPFR_RNDN);
		}
		if(mpfr_regular_p(rangeFinish)){
		    real_to_frac(rangeFinishFrac, rangeFinish);
		}else{
		    mpfr_set_d(rangeFinishFrac, 0.0, MPFR_RNDN);
		}
	
		isRangeFlat = isFlat(rangeStart,rangeFinish);
		isDomainFlat = isFlat(domainStart, domainFinish);

		if(parent->m_isDomainScaled){
		    mpfr_set(domainOffset, domainStartFrac, MPFR_RNDN);
		    domainScale=0;
		    mpfr_set_zero(domainStartScaled,+1);
		    mpfr_sub(domainFinishScaled, domainFinishFrac, domainStartFrac, MPFR_RNDN);

		    mpfr_fprintf(stderr, "  [%Re,%Re], %Re\n", domainStartFrac, domainFinishFrac, domainFinishScaled);

		    if(isDomainFlat && mpfr_regular_p(rangeStart) && mpfr_regular_p(domainFinishScaled) && !mpfr_equal_p(domainStart,domainFinish)){
			assert(mpfr_cmp_d(domainFinishScaled,0.0)>=0);
			assert(mpfr_cmp_d(domainFinishScaled,1.0)<0);

			while(mpfr_cmp_d(domainFinishScaled,0.5)<0){
			    ++domainScale;
			    mpfr_mul_d(domainFinishScaled,domainFinishScaled,2.0, MPFR_RNDN);
			}
		    }						     
		}else{
		    mpfr_set_zero(domainOffset, +1);
		    domainScale=0;
		    mpfr_set(domainStartScaled, domainStartFrac, MPFR_RNDN);
		    mpfr_set(domainFinishScaled, domainFinishFrac, MPFR_RNDN);
		}

		// Check for constant functions
		isRangeConstant=false;
		if( mpfr_cmp(domainStart, domainFinish)==0 ){
		    // If the domain is degenerate, range must be constant
		    isRangeConstant=true;
		}else if(isRangeFlat){
		    // Maybe the endpoints are equal?
		    // TODO - This is pessimistic, as it ignores two endpoints that
		    // could be rounded faithfully to a constant.
		    if(mpfr_cmp(rangeStart, rangeFinish)==0){

			// Ok, so they are the same at the ends, part there might be
			// a difference in the middle
			mpfr_t rangeMid, domainMid;
			mpfr_init2(domainMid,parent->m_domainWF+1);
			mpfr_init2(rangeMid,parent->m_rangeWF+1);

			parent->eval_clamp(rangeMid, domainMid);

			if(mpfr_cmp(rangeMid,rangeStart)==0){
			    // Ok, we definitely have a line, just double-check that
			    // there isn't a curve... (Really? How likely)

			    sollya_node_ptr_t deriv=make_shared_sollya(differentiate(parent->m_function.getSollyaNode()));

			    sollya_range_t rr={ &domainStart, &domainFinish };
			    sollya_chain_t zeros=fpFindZerosFunction(deriv.get(), rr, getToolPrecision());

			    if(zeros==NULL){
				isRangeConstant=true;
			    }else{
				fprintf(stderr, "Found flat segment that looks constant, but no zeros.");
				// Needs to be investigated.
				exit(1);
			    }

			    freeChain(zeros, freeMpfrPtr);
			}

			mpfr_clears(domainMid, rangeMid, NULL);
		    }
		}
	    }

	    Segment::Segment(Range *_parent, mpfr_t _domainStart, mpfr_t _domainFinish)
		: parent(_parent)
	    {
		mpfr_inits2(parent->m_domainWF+1, domainStart, domainFinish, domainStartFrac, domainFinishFrac, (mpfr_ptr)0);
		mpfr_inits2(parent->m_rangeWF+1, rangeStart, rangeFinish, rangeStartFrac, rangeFinishFrac, (mpfr_ptr)0);

		mpfr_inits2(parent->m_domainWF+1, domainOffset, domainStartScaled,  domainFinishScaled, (mpfr_ptr)0);
		domainScale=0;
	
		set_domain(_domainStart, _domainFinish);
	    }

	    Segment::Segment(const Segment &o)
		: parent(o.parent)
	    {
		mpfr_inits2(parent->m_domainWF+1, domainStart, domainFinish, (mpfr_ptr)0);
		mpfr_inits2(parent->m_domainWF+1, domainStartFrac, domainFinishFrac, (mpfr_ptr)0);
		mpfr_inits2(parent->m_rangeWF+1, rangeStart, rangeFinish, rangeStartFrac, rangeFinishFrac, (mpfr_ptr)0);
		mpfr_inits2(parent->m_domainWF+1, domainOffset, domainStartScaled,  domainFinishScaled, (mpfr_ptr)0);
	
		mpfr_set(domainStart, o.domainStart, MPFR_RNDN);
		mpfr_set(domainFinish, o.domainFinish, MPFR_RNDN);
		mpfr_set(domainStartFrac, o.domainStartFrac, MPFR_RNDN);
		mpfr_set(domainFinishFrac, o.domainFinishFrac, MPFR_RNDN);
		mpfr_set(rangeStart, o.rangeStart, MPFR_RNDN);
		mpfr_set(rangeFinish, o.rangeFinish, MPFR_RNDN);
		mpfr_set(rangeStartFrac, o.rangeStartFrac, MPFR_RNDN);
		mpfr_set(rangeFinishFrac, o.rangeFinishFrac, MPFR_RNDN);

		mpfr_set(domainOffset, o.domainOffset, MPFR_RNDN);
		domainScale=o.domainScale;
		mpfr_set(domainStartScaled, o.domainStartScaled, MPFR_RNDN);
		mpfr_set(domainFinishScaled, o.domainFinishScaled, MPFR_RNDN);
	
		isRangeFlat=o.isRangeFlat;
		isDomainFlat=o.isDomainFlat;
		isRangeConstant=o.isRangeConstant;
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
		    mpfr_set(rangeStartFrac, o.rangeStartFrac, MPFR_RNDN);
		    mpfr_set(rangeFinishFrac, o.rangeFinishFrac, MPFR_RNDN);
		    mpfr_set(domainOffset, o.domainOffset, MPFR_RNDN);
		    domainScale=o.domainScale;
		    mpfr_set(domainStartScaled, o.domainStartScaled, MPFR_RNDN);
		    mpfr_set(domainFinishScaled, o.domainFinishScaled, MPFR_RNDN);
		    isRangeFlat=o.isRangeFlat;
		    isDomainFlat=o.isDomainFlat;
		    isRangeConstant=o.isRangeConstant;
		    properties=o.properties;
		}
		return *this;
	    }

	    Segment::~Segment()
	    {
		mpfr_clears(domainStart, domainFinish,rangeStart, rangeFinish, domainStartFrac, domainFinishFrac, rangeStartFrac, rangeFinishFrac, domainOffset, domainStartScaled, domainFinishScaled, (mpfr_ptr)0);
	    }

	    void Segment::dump(FILE *dst) const
	    {
		mpfr_fprintf(dst, "%s[%Re,%Re] -> %s[%Re,%Re], off=%Re, scale=%d",
			     isDomainFlat?"flat":"split", domainStart, domainFinish,
			     isRangeFlat?"flat":"split", rangeStart, rangeFinish,
			     domainOffset, domainScale
			     );
	    }

	    bool Segment::has_property(std::string name) const
	    {
		return properties.find(name)!=properties.end();
	    }

	    /* This function is flat, with eD=binade(domainStart), eR=binade(f(domainStart))
	       That means that 	2^(eD-1) <= domainStart <= domainFinish < 2^eD
	       and 2^(eR-1) <= f(domainStart) <= f(domainFinish) < 2^eR
	       we want to convert to fixed-point, so we have inputs in [0,1) with domainWF
	       fractional bits, and outputs in [0,1) with rangeWF fractional bits.
	
	       To move a real x to fixed-point we have:
	       1.0*2^eD <= x < 2.0*2^eD
	       1.0 <= x * 2^-eD < 2.0
	       0 <= x * 2^-eD - 1.0 < 1.0
	       0 <= 2*(x * 2^-eD - 1.0) < 2.0
	       0 <= x*2^(-eD+1) -2.0 < 2.0
	       0 <= x*2^(-eD) - 1.0 < 1.0
	    */
	    sollya_node_t Range::get_scaled_flat_function(binade_t bDom, binade_t bRan)
	    {	
		sollya_node_t f=m_function.getSollyaNode();
		sollya_node_t tx=NULL;

		// I have no brain!!!!

	    

		if(bDom.first==NegNormal){
		    if(bRan.first==PosInf){
			tx=makeConstantDouble(INFINITY);
		    }else if(bRan.first==PosNormal){
			tx=makeAdd(makeVariable(), makeConstantDouble(1.0));
			tx=makeMul(tx, makeConstantDouble(pow(2.0, bDom.second)));
		
			tx=makeNeg(tx);
			tx=substitute(f, tx);
		
			tx=makeMul(tx, makeConstantDouble(pow(2.0, -bRan.second)));
			tx=makeSub(tx, makeConstantDouble(1.0));
		    }else if(bRan.first==PosZero || bRan.first==NegZero){
			tx=makeConstantDouble(bRan.first==PosZero ? 0.0 : -0.0);
		    }else if(bRan.first==NegNormal){
			tx=makeAdd(makeVariable(), makeConstantDouble(1.0));
			tx=makeMul(tx, makeConstantDouble(pow(2.0, bDom.second)));

			tx=makeNeg(tx);
			tx=substitute(f, tx);
			tx=makeNeg(tx);

			tx=makeMul(tx, makeConstantDouble(pow(2.0, -bRan.second)));
			tx=makeSub(tx, makeConstantDouble(1.0));
		    }else if(bRan.first==NegInf){
			tx=makeConstantDouble(-INFINITY);
		    }else if(bRan.first==NaN){
			tx=makeConstantDouble(nan(""));
		    }
		}else if((bDom.first==PosZero) || (bDom.first==NegZero)){
		    if(bRan.first==PosInf){
			tx=makeConstantDouble(INFINITY);
		    }else if((bRan.first==PosNormal)){
			tx=makeMul(makeConstantDouble(0.0),makeVariable());
			tx=substitute(f, tx);
			tx=makeMul(tx, makeConstantDouble(pow(2.0, -bRan.second)));
			tx=makeSub(tx, makeConstantDouble(1.0));
		    }else if(bRan.first==PosZero){
			tx=makeConstantDouble(0.0);
		    }else if(bRan.first==NegZero){
			tx=makeConstantDouble(-0.0);
		    }else if(bRan.first==NegInf){
			tx=makeConstantDouble(-INFINITY);
		    }else if(bRan.first==NaN){
			tx=makeConstantDouble(nan(""));
		    }
		}else if(bDom.first==PosNormal){
		    if(bRan.first==PosInf){
			tx=makeConstantDouble(INFINITY);
		    }else if(bRan.first==PosNormal){
			tx=makeAdd(makeVariable(), makeConstantDouble(1.0));
			tx=makeMul(tx, makeConstantDouble(pow(2.0, bDom.second)));
		
			tx=substitute(f, tx);
		
			tx=makeMul(tx, makeConstantDouble(pow(2.0, -bRan.second)));
			tx=makeSub(tx, makeConstantDouble(1.0));
		    }else if(bRan.first==PosZero || bRan.first==NegZero){
			tx=makeConstantDouble(bRan.first==PosZero ? 0.0 : -0.0);
		    }else if(bRan.first==NegNormal){
			tx=makeAdd(makeVariable(), makeConstantDouble(1.0));
			tx=makeMul(tx, makeConstantDouble(pow(2.0, bDom.second)));

			tx=substitute(f, tx);
			tx=makeNeg(tx);

			tx=makeMul(tx, makeConstantDouble(pow(2.0, -bRan.second)));
			tx=makeSub(tx, makeConstantDouble(1.0));
		    }else if(bRan.first==PosInf){
			tx=makeConstantDouble(INFINITY);
		    }else if(bRan.first==NaN){
			tx=makeConstantDouble(nan(""));
		    }
		}

		if(tx==NULL){
		    std::stringstream acc;
		    acc<<"get_scaled_flat_function("<<bDom<<","<<bRan<<") - Not implemented.";
		    throw std::invalid_argument(acc.str());
		}

	

		return tx;
	    }

	    sollya_node_t Segment::get_scaled_flat_function()
	    {
		assert(isRangeFlat && isDomainFlat);
		auto bDom=binade(domainStart);
		auto bRan=binade(rangeStart);

		if(parent->m_isDomainScaled){
		    // Gives us a function over [domainStartFrac,domainFinishFrac]
		    sollya_node_t raw=parent->get_scaled_flat_function(bDom,bRan);
		    // Initial scaling that moves it into that range
		    sollya_node_t pre=makeMul(makeVariable(), makeConstantDouble(ldexp(1.0, -domainScale)));
		    pre=makeAdd(pre, makeConstant(domainOffset));
		    return substitute(raw, pre);
		}else{
		    return parent->get_scaled_flat_function(bDom, bRan);
		}
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


	    // Provided by sollya
	    extern "C" int getVerbosity();
	    extern "C" int setVerbosity(int);

	    void infNorm(mpfr_t error, sollya_node_t func, sollya_node_t poly, mpfr_t a, mpfr_t b, unsigned points)
	    {
		sollya_node_t diff=makeSub(copyTree(func), copyTree(poly));

		uncertifiedInfnorm(error, diff, a, b, points, getToolPrecision());
		free_memory(diff);
	    }

	    bool Segment::check_poly(sollya_node_t func, sollya_node_t poly)
	    {
		mpfr::mpreal error=mpfr::create_zero(getToolPrecision());
		infNorm(get_mpfr_ptr(error), func, poly, domainStartScaled, domainFinishScaled, 71);
		if(mpfr_cmp_d(get_mpfr_ptr(error), ldexp(1.0, -(int)parent->m_rangeWF-10 ))<0){
		    return true;
		}else{
		    return false;
		}
	    }

	    sollya_node_t Segment::make_linear_poly(sollya_node_t func, unsigned degree)
	    {
		unsigned prec=getToolPrecision();

		mpfr::mpreal x0(domainStartScaled), x1(domainFinishScaled);

		x0.backend().precision(prec);
		x1.backend().precision(prec);

		mpfr::mpreal y0=mpfr::create_zero(prec), y1=mpfr::create_zero(prec);
		evaluateFaithful(get_mpfr_ptr(y0), func, get_mpfr_ptr(x0), prec);
		evaluateFaithful(get_mpfr_ptr(y1), func, get_mpfr_ptr(x1), prec);

		mpfr::mpreal m=(y0 - y1)/(x0 - x1);
		mpfr::mpreal c=(x0*y1 - x1*y0)/(x0 - x1);

		mpfr_t *coeffs=(mpfr_t*)malloc(sizeof(mpfr_t)*(degree+1));
		for(unsigned i=0;i<=degree;i++){
		    mpfr_init2(coeffs[i], prec);
		    mpfr_set_d(coeffs[i], 0.0, MPFR_RNDN);
		}

		//mpfr_fprintf(stderr, "x0=%Re, x1=%Re, y0=%Re, y1=%Re\n", get_mpfr_ptr(x0), get_mpfr_ptr(x1), get_mpfr_ptr(y0), get_mpfr_ptr(y1));

		mpfr_set(coeffs[0], get_mpfr_ptr(c), MPFR_RNDN);
		mpfr_set(coeffs[1], get_mpfr_ptr(m), MPFR_RNDN);

		sollya_node_t res=makePolynomial(&coeffs[0], degree);
		for(unsigned i=0;i<=degree;i++){
		    mpfr_clear(coeffs[i]);
		}
		free(coeffs);

		return res;
	    }

	    sollya_node_t Segment::minimax(unsigned degree)
	    {
		assert(isRangeFlat && isDomainFlat);
	
		sollya_node_t res=NULL;

		mpfr_fprintf(stderr, " minimax, dom=[%Re,%Re], fDom=[%Re,%Re] -> ran=[%Re,%Re], fRan=[%Re,%Re], domFlat=%s, ranFlat=%s, ranConst=%s\n", domainStart, domainFinish, domainStartFrac, domainFinishFrac, rangeStart, rangeFinish, rangeStartFrac, rangeFinishFrac, isDomainFlat?"Y":"N", isRangeFlat?"Y":"N", isRangeConstant?"Y":"N");
	
		sollya_node_t func=get_scaled_flat_function();
		fprintTree(stderr, func);
		fprintf(stderr, "\n");
		sollya_node_t weight=makeConstantDouble(1.0);
	
		mp_prec_t prec=getToolPrecision();
		mpfr_t requestedQuality;
		mpfr_init2(requestedQuality, prec);
		mpfr_set_d(requestedQuality, pow(2.0, -parent->m_rangeWF-10), MPFR_RNDN);
	
		mpfr::mpreal zero;

		if(mpfr_equal_p(domainStart, domainFinish) || isRangeConstant){
		    mpfr_t *coeffs=(mpfr_t*)malloc(sizeof(mpfr_t)*(degree+1));
		    for(unsigned i=0;i<=degree;i++){
			mpfr_init2(coeffs[i], getToolPrecision());
			mpfr_set_d(coeffs[i], 0.0, MPFR_RNDN);
		    }
		    if(mpfr_regular_p(rangeStart)){
			evaluateFaithful(coeffs[0], func, domainStartScaled, prec);
		    }else{
			mpfr_set(coeffs[0], rangeStartFrac, MPFR_RNDN);
		    }
		    mpfr_fprintf(stderr, "domainStartScaled=%Re, coeff[0]=%Re\n", domainStartScaled, coeffs[0]);

		    //mpfr_fix(coeffs[0], parent->m_rangeWF);
		    res=makePolynomial(&coeffs[0], degree);
		    for(unsigned i=0;i<=degree;i++){
			mpfr_clear(coeffs[i]);
		    }
		    free(coeffs);
		}else{	
		    unsigned degreeCurr=degree;

		    sollya_chain_t monom=makeIntPtrChainFromTo(0, degreeCurr);

		    sollya_node_t poly=make_linear_poly(func, degree);

		    if(check_poly(func,poly)){
			res=poly;
		    }else{
			bool failedRemez=false;
			jmp_buf env;
			if(setjmp (env)){
			    failedRemez=true;
			    fprintf(stderr, "Failed");
			    exit(1);

			}else{
			    setRecoverEnvironment(&env);

			    res=remez(func, weight, monom, domainStartScaled,domainFinishScaled, &requestedQuality, prec);
			}

			invalidateRecoverEnvironment();

			if(failedRemez){
			    res=poly;
			}else{
			    free_memory(poly);
			}
		    }


		    freeChain(monom,freeIntPtr);
		}
		
		free_memory(func);
		free_memory(weight);
	
		assert(getDegree(res)<=(int)degree);
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
	
		mpfr_fprintf(stderr, " fpminimax, dom=[%Re,%Re] -> range [%Rf,%Re], fDom=[%Rf,%Rf], fRange=[%Rf,%Rf]\n", domainStart, domainFinish, rangeStart, rangeFinish, domainStartFrac, domainFinishFrac, rangeStartFrac, rangeFinishFrac);
	
		sollya_node_t func=get_scaled_flat_function();
	
		mp_prec_t prec=getToolPrecision();
	
		if(mpfr_equal_p(domainStart, domainFinish) || isRangeConstant){
		    mpfr_t *coeffs=(mpfr_t*)malloc(sizeof(mpfr_t)*(degree+1));
		    for(unsigned i=0;i<=degree;i++){
			mpfr_init2(coeffs[i], getToolPrecision());
			mpfr_set_d(coeffs[i], 0.0, MPFR_RNDN);
		    }

		    evaluateFaithful(coeffs[0], func, domainStartScaled, prec);
		    	  
		    // We need to round to whatever the format of coeff zero should be
		    mpfr_fix(coeffs[0], coeffWidths[0]);

		    //	  mpfr_set(coeffs[0], rangeStartFrac, MPFR_RNDN);
		
		    res=makePolynomial(&coeffs[0], degree);
		    for(unsigned i=0;i<=degree;i++){
			mpfr_clear(coeffs[i]);
		    }
		    free(coeffs);
	  
		}else{

		    sollya_node_t constantPart=makeConstantDouble(0.0);

		    /* This is now updated to walk up from low degree to
		       high degree, to try to avoid the creation of exceptionally
		       precise polynomials which have very large constant for
		       high coefficients. As soon as we hit something that has
		       more than two spare bits we quit. This helps us get quadratics
		       with reasonable constants, rather than quintics with crazy
		       ones.
		    */
		

		    bool success=false;
		    for(unsigned currDegree=1;currDegree<=degree;currDegree++){

			// The precision has to be wacked _way_ up here, as otherwise this just
			// hangs on some polynomials. This is not a great way of doing things...
			setToolPrecision(max((int)getToolPrecision(),1024));

			sollya_chain_t monom=makeIntPtrChainFromTo(0, currDegree);
			std::vector<int> localCoeffWidths(coeffWidths.begin(), coeffWidths.begin()+currDegree+1);
			sollya_chain_t formats=makeIntPtrChain(localCoeffWidths);

			// Move bounds to higher precision, does it help?
			mpfr_t a, b;
			mpfr_inits2(getToolPrecision(), a, b, NULL);
			mpfr_set(a, domainStartScaled, MPFR_RNDN);
			mpfr_set(b, domainFinishScaled, MPFR_RNDN);
			
			if(1){
			    res=FPminimax(func, monom, formats, /*points*/NULL, a, b, FIXED, ABSOLUTESYM, constantPart, minimax);
			}else{
			    res=dirtierFPminimax(func, monom, formats, /*points*/NULL, a, b, FIXED, ABSOLUTESYM, constantPart, minimax);
			}
			setToolPrecision(prec);

			mpfr_clears(a,b, NULL);

			freeChain(monom,freeIntPtr);
			freeChain(formats,freeIntPtr);

			mpfr_t error;
			mpfr_init2(error, getToolPrecision());
			infNorm(error, func, res, domainStartScaled, domainFinishScaled, 71);

			if(mpfr_cmp_d(error, ldexp(1.0, -parent->m_rangeWF-3))<0){
			    // Error is good enough, can stop
			    success=true;
			    if(currDegree!=degree)
				mpfr_fprintf(stderr, "\n Exited with degree=%d and error=%Rg\n", currDegree, error);
			}

			mpfr_clear(error);

			if(success)
			    break;
		    }

		
		    free_memory(constantPart);
		}

		free_memory(func);
	
		return res;
	    }

	    sollya_node_t Segment::fpminimax(unsigned degree, int coeffWidths, sollya_node_t minimax)
	    {
		std::vector<int> coeffWidthsVec(degree+1, coeffWidths);
	
		return fpminimax(coeffWidthsVec, minimax);
	    }

	    NumberClass mpfr_class(mpfr_t x)
	    {
		assert(!mpfr_nan_p(x));
		if(mpfr_inf_p(x))
		    return mpfr_sgn(x)<0 ? NegInf : PosInf;
		if(mpfr_zero_p(x))
		    return mpfr_sgn(x)<0 ? NegZero : PosZero;
		return mpfr_sgn(x)<0 ? NegNormal : PosNormal;
	    }

	    // Maximum normal value in the domain
	    void Range::createDomainMax(mpfr_t val)
	    {
		// TODO : Is this actually correct?
		mpfr_set_d(val, 1.0, MPFR_RNDN);
		mpfr_nextbelow(val);
		mpfr_mul_2si(val, val, (1<<m_domainWE)/2-1, MPFR_RNDN);
	    }

	    // Minimum normal value in the domain
	    void Range::createDomainMin(mpfr_t val)
	    {
		// TODO : Is this actually correct?
		mpfr_set_d(val, 0.5, MPFR_RNDN);
		mpfr_mul_2si(val, val, -(1<<m_domainWE)/2+1, MPFR_RNDN);
	    }


	    Range::Range(const Function &f, int domainWF, int rangeWF, mpfr_t domainStart, mpfr_t domainFinish, int domainWE, int rangeWE, bool isDomainScaled)
		: m_function(f)
		, m_domainWF(domainWF)
		, m_rangeWF(rangeWF)
		, m_domainWE(domainWE)
		, m_rangeWE(rangeWE)
		, m_isDomainScaled(isDomainScaled)
	    {
		//if(mpfr_sgn(domainStart)<0)
		//	throw std::invalid_argument("Domain must be strictly positive.");
		if(mpfr_greater_p(domainStart, domainFinish))
		    throw std::invalid_argument("Domain must be ordered.");
	
		mpfr_init2(m_domainFractionEnd, domainWF+1);
		mpfr_init2(m_rangeFractionEnd, rangeWF+1);
	
		// TODO: What are these for?
		mpfr_set_d(m_domainFractionEnd, 0.5, MPFR_RNDN);
		mpfr_set_d(m_rangeFractionEnd, 0.5, MPFR_RNDN);
	
		mpfr_t tmpStart, tmpFinish;
		mpfr_inits2(domainWF+1, tmpStart, tmpFinish, NULL);

		NumberClass classStart=mpfr_class(domainStart);
		NumberClass classFinish=mpfr_class(domainFinish);

		if(classStart==NegInf){
		    // Create a segment for infinity...
		    m_segments.push_back(Segment(this, domainStart, domainStart));
		}

		if(classStart<=NegNormal && NegNormal<=classFinish){
		    if(classStart<NegNormal){
			createDomainMax(tmpStart);
			mpfr_neg(tmpStart, tmpStart, MPFR_RNDN);
		    }else{
			mpfr_set(tmpStart, domainStart, MPFR_RNDN);
		    }
		    if(NegNormal < classFinish){
			createDomainMin(tmpFinish);
			mpfr_neg(tmpFinish, tmpFinish, MPFR_RNDN);
		    }else{
			mpfr_set(tmpFinish, domainFinish, MPFR_RNDN);
		    }
		    m_segments.push_back(Segment(this, tmpStart, tmpFinish));
		}

		if(classStart<=NegZero && NegZero<=classFinish){
		    mpfr_set_d(tmpStart, 0.0, MPFR_RNDN);
		    mpfr_neg(tmpStart, tmpStart, MPFR_RNDN);	// Make negative zero
		    m_segments.push_back(Segment(this, tmpStart, tmpStart));
		}

		if(classStart<=PosZero && PosZero<=classFinish){
		    mpfr_set_d(tmpStart, 0.0, MPFR_RNDN);
		    m_segments.push_back(Segment(this, tmpStart, tmpStart));
		}

		if(classStart<=PosNormal && PosNormal<=classFinish){
		    if(classStart<PosNormal){
			createDomainMin(tmpStart);
		    }else{
			mpfr_set(tmpStart, domainStart, MPFR_RNDN);
		    }
		    if(PosNormal<classFinish){
			createDomainMax(tmpFinish);
		    }else{
			mpfr_set(tmpFinish, domainFinish, MPFR_RNDN);
		    }
		    m_segments.push_back(Segment(this, tmpStart, tmpFinish));
		}

		if(classFinish==PosInf){
		    m_segments.push_back(Segment(this, domainFinish, domainFinish));
		}

		mpfr_clears(tmpStart, tmpFinish, NULL);
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
		mpfr_init2(mid, m_domainWF+1);
		mpfr_add(mid, src->domainStart, src->domainFinish, MPFR_RNDD); // note that we round down here
		mpfr_mul_2si(mid, mid, -1, MPFR_RNDN);
	
		split_segment(src, mid);

		mpfr_clear(mid);
	    }



	    void Range::make_monotonic_or_range_flat()
	    {
		sollya_node_ptr_t deriv=make_shared_sollya(differentiate(m_function.getSollyaNode()));
	
		if(::flopoco::verbose>=DEBUG){
		    fprintf(stderr, "Deriv = ");
		    fprintTree(stderr, deriv.get());
		    fprintf(stderr, "\n");
		}
	
		sollya_range_t rr={ &(m_segments.begin()->domainStart), &(boost::prior(m_segments.end())->domainFinish) };
		sollya_chain_t zeros=fpFindZerosFunction(deriv.get(), rr, getToolPrecision());
	
		mpfr_t extremaLow, extremaHigh, extremaLowRange, extremaHighRange;
		mpfr_inits2(m_domainWF+1, extremaLow, extremaHigh, (mpfr_ptr)0);
		mpfr_inits2(m_rangeWF+1, extremaLowRange, extremaHighRange, (mpfr_ptr)0);
	
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
			eval_clamp(extremaLowRange, extremaLow);
			eval_clamp(extremaHighRange, extremaHigh);
			mpfr_fprintf(stderr, "    extremaLow=%Re, extremaHigh=%Re\n", extremaLow, extremaHigh);
			fprintf(stderr,"    lowInDomain=%d, highInDomain=%d\n", is_in_domain(extremaLow)?1:0, is_in_domain(extremaHigh)?1:0);
			
			// The function has an extrema at extremaLow < extrema < extremaHigh
			if((!is_in_domain(extremaLow)) || (!is_in_domain(extremaHigh))){
			    // We don't have to worry, as one or both of the extrema turns out
			    // not to be representable. If one is representable, then they other
			    // must be curving up or down quite happily
			    fprintf(stderr, "    one side of extrema not in domain.\n");
			}else{
			    segment_it_t seg=Range::find_segment(extremaLow);
			    mpfr_fprintf(stderr, "    [%Re,%Re]\n", seg->domainStart, seg->domainFinish);
				
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
		mpfr_init2(newFinish, m_domainWF+1);
	
		segment_it_t curr=m_segments.begin();
		while(curr!=m_segments.end()){
		    if(!curr->isDomainFlat){
			fprintf(stderr, "Flattening: ");
			curr->dump(stderr);
			fprintf(stderr, "\n");


			if(mpfr_sgn(curr->domainStart)==0){
			    throw std::string("Constraint violation - zero segment should be flat.");
			}else if(mpfr_sgn(curr->domainStart)<0){
			    // negative segment
			    if(mpfr_sgn(curr->domainFinish)>=0)
				throw std::string("Constraint violation - by this point segments should not span number classes.");

			    assert(mpfr_get_exp(curr->domainStart) > mpfr_get_exp(curr->domainFinish));
			    // We have [start,finish] with binade(start) > binade(finish) and start and finish both negative normal
			    // We need to find newFinish where  start<=newFinish and binade(start)==binade(newFinish)
			    int eStart=mpfr_get_exp(curr->domainStart);
			    // To get the endpoint, we can just do prev(2^(e+1))

			    mpfr_set_si_2exp(newFinish, -1, eStart-1, MPFR_RNDN);	// This is the smallest number in the binade

			    split_segment(curr, newFinish);
			}else{
			    // Positive segment
			    if(mpfr_sgn(curr->domainFinish)<=0)
				throw std::string("Constraint violation - by this point segments should not span number classes.");

			    // So currently we have [start,finish] with binade(start) < binade(finish) and start and finish both positive normal
			    assert(mpfr_get_exp(curr->domainStart) < mpfr_get_exp(curr->domainFinish));
			    // We need to find newFinish where  start<=newFinish and binade(start)==binade(newFinish)
			    int eStart=mpfr_get_exp(curr->domainStart);
			    // To get the endpoint, we can just do prev(2^(e+1))

			    mpfr_set_si_2exp(newFinish, 1, eStart, MPFR_RNDN);
			    mpfr_nextbelow(newFinish);	// This is the largest number in the binade

			    split_segment(curr, newFinish);
			}
			
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
		int lowS, highS, midS;
	
		mpfr_init_copy(low, start);
		mpfr_init_copy(high, finish);
		mpfr_init2(mid, m_domainWF+1);
		mpfr_init2(val, m_rangeWF+1);
	
		eval_clamp(val, low);
		lowE=mpfr_get_exp(val);
		lowS=mpfr_sgn(val);
		eval_clamp(val, high);
		highE=mpfr_get_exp(val);
		highS=mpfr_sgn(val);
	
		assert((lowE!=highE) || (lowS!=highS));
	
		while(true){
		    // mpfr_fprintf(stderr, "  low=%Re, high=%Re, lowE=%d lowSgn=%d, highE=%d highS=%d\n", low, high, lowE, lowS, highE, highS);
		
		    int finish=false;
		    mpfr_nextabove(low);
		    finish = mpfr_equal_p(low, high) != 0;
		    mpfr_nextbelow(low);
		
		    if(finish)
			break;
		
		    mpfr_add(mid, low, high, MPFR_RNDN);
		    mpfr_div_2ui (mid, mid, 1, MPFR_RNDN);
		
		    eval_clamp(val, mid);
		    midE=mpfr_get_exp(val);
		    midS=mpfr_sgn(val);
		
		    if(lowE==midE && lowS==midS){
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
		mpfr_init2(newFinish, m_domainWF+1);
	
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

		// Prune any out of range segments
		curr=m_segments.begin();
		while(curr!=m_segments.end()){
		    auto next=curr;
		    ++next;

		    assert(curr->isRangeFlat);
		    auto bRan=binade(curr->rangeStart);
		    if((bRan.first==PosNormal || bRan.first==NegNormal) &&  (bRan.second<((-1<<(m_rangeWE-1))+1))){
			throw std::string("Underflowing segment to -infinity - TODO, what should happen here?.\n");
			//m_segments.erase(curr);
		    }

		    curr=next;
		}
	
		mpfr_clear(newFinish);
	    }

	    void Range::eval(mpfr_t res, mpfr_t x)
	    {
		m_function.eval(res, x);
	    }

	    void Range::eval_clamp(mpfr_t res, mpfr_t x)
	    {
		m_function.eval(res, x);
		auto b=binade(res);
		if((b.first==PosNormal) || (b.first==NegNormal)){
		    if(b.second<((-1<<(m_rangeWE-1))+1)){
			mpfr_set_zero(res, b.first==PosNormal?+1:-1);
		    }else if(b.second>=(1<<(m_rangeWE-1))){
			mpfr_set_inf(res, b.first==PosNormal?+1:-1);
		    }
		}
	    }

	    void Range::init_eval(mpfr_t res, mpfr_t x)
	    {
		mpfr_init2(res, m_rangeWF+1);
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
		fprintf(dst, "Total segments = %u\n", (unsigned)m_segments.size());
	    }

	    // TODO : This is O(n)... fix it
	    bool Range::is_in_domain(mpfr_t x)
	    {
		bool res=false;
	
		mpfr_t tmp;
		mpfr_init2(tmp, m_domainWF+1);
		if(0==mpfr_set(tmp, x, MPFR_RNDN)){
		    segment_it_t curr=m_segments.begin();
		    while(curr!=m_segments.end()){
			if(mpfr_lessequal_p(curr->domainStart, tmp) && mpfr_lessequal_p(tmp,curr->domainFinish)){
			    res=true;
			    break;
			}				
			//mpfr_fprintf(stderr, "  not %Re in [%Re,%Re]\n", tmp, curr->domainStart, curr->domainFinish);
			++curr;
		    }
		}else{
		    mpfr_fprintf(stderr, "is_in_domain(%Re) - rounding results in %Re\n", x, tmp);
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

	    Range::segment_it_t Range::get_segment(unsigned index)
	    {
		if(index>=m_segments.size())
		    throw std::invalid_argument("get_segment - index out of range.");

		segment_it_t curr=m_segments.begin();
		while(index>0){
		    --index;
		    ++curr;
		}

		return curr;
	    }

	}; // float_approx
    }; // random
}; // flopoco

