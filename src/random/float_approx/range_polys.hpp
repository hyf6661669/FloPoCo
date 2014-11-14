#ifndef flopoco_random_float_approx_range_polys_hpp
#define flopoco_random_float_approx_range_polys_hpp

#include "range.hpp"

#include <boost/smart_ptr.hpp>
#include <boost/any.hpp>
#include <boost/format.hpp>
#include <boost/utility.hpp>

#include "random/utils/mpfr_vec.hpp"
#include "random/utils/comparable_float_type.hpp"

#include "FixFunctions/PolynomialEvaluator.hpp"
#include "random/poly/fixed_point_polynomial_evaluator.hpp"

#include "static_quantiser.hpp"

#include "mpreal.h"


namespace flopoco
{
    namespace random
    {

	FixedPointPolynomialEvaluator *CreateRoundingPolynomialEvaluator(
									 const fixed_format_t &inputFormat,
									 const std::vector<fixed_format_t> &coefficientFormats,
									 int outputLsb,
									 int guardBits,
									 const mpfr::mpreal &errorBudget,
									 Target *target
									 );
	


	namespace float_approx
	{

	    MPFRVec getPolyCoeffs(sollya_node_t poly, unsigned degree)
	    {
		MPFRVec res(degree+1, getToolPrecision());
		for(unsigned i=0;i<=degree;i++){
		    sollya_node_t c=getIthCoefficient(poly, i);
		    evaluateConstantExpression(res[i], c, getToolPrecision());
		    free_memory(c);
		}
		return res;
	    }

	    class RangePolys
	    {
	    public:	
		typedef Range::segment_it_t segment_it_t;
	    private:
	
		Range *m_range;
		unsigned m_degree;

		void infNorm(mpfr_t error, sollya_node_t func, sollya_node_t poly, mpfr_t a, mpfr_t b, unsigned points)
		{
		    sollya_node_t diff=makeSub(copyTree(func), copyTree(poly));
		
		    uncertifiedInfnorm(error, diff, a, b, points, getToolPrecision()); 
		    free_memory(diff);
		}

		void update_segment_minimax(segment_it_t curr)
		{
		    if(curr->has_property("minimax") && curr->has_property("minimax_error"))
			return;
		
		    sollya_node_ptr_t poly=make_shared_sollya(curr->minimax(m_degree));
		    sollya_node_t func=curr->get_scaled_flat_function();

		    MPFRVec error(1, getToolPrecision());
		
		    // TODO : THis should only be applied if the domain is not degenerate, to avoid warnings

	    
		    if(mpfr_equal_p(curr->domainStart,curr->domainFinish) || !mpfr_regular_p(curr->rangeStart)){
			if(mpfr_nan_p(curr->rangeStart)){
			    // THe output will be NaN
			    mpfr_set_d(error[0], 0.0, MPFR_RNDN);
			}else if(mpfr_inf_p(curr->rangeStart)){
			    // Output will be +-inf
			    mpfr_set_d(error[0], 0.0, MPFR_RNDN);
			}else{
			    mpfr_t correct, got;
			    mpfr_init2(correct, getToolPrecision());
			    mpfr_init2(got, getToolPrecision());
		  
			    evaluateFaithful(correct, func, curr->domainStartScaled, getToolPrecision());
			    evaluateFaithful(got, poly.get(), curr->domainStartScaled, getToolPrecision());
			    mpfr_sub(error[0], correct, got, MPFR_RNDN);
			    mpfr_abs(error[0],error[0], MPFR_RNDN);

			    mpfr_fprintf(stderr, "  minimax_error dom=[%Re,%Re], got=%Re, correct=%Re, err=%Re\n", curr->domainStart, curr->domainFinish, got, correct, error[0]);

			    if(mpfr_cmp_d(error[0],0.1)>0 || mpfr_nan_p(error[0])){
				throw std::string("Error is extremely large or nan. What to do?");
			    }
		
			    mpfr_clear(correct);
			    mpfr_clear(got);
			}
		    }else{

			infNorm(error[0], func, poly.get(), curr->domainStartScaled, curr->domainFinishScaled, 71);
		    }

		    if(mpfr_regular_p(curr->rangeStart)){
			if(mpfr_nan_p(error[0])){
			    fprintf(stderr, "\nfunc = ");
			    fprintTree(stderr, func);
			    fprintf(stderr, "\n\npoly = ");
			    fprintTree(stderr, poly.get());
			    mpfr_fprintf(stderr, "\n\n sDom=[%Re,%Re]\n", curr->domainStartScaled, curr->domainFinishScaled);
			}
		    }

		    assert(!mpfr_nan_p(error[0]));
		
		    curr->properties["minimax"]=poly;
		    curr->properties["minimax_error"]=error;

		    free_memory(func);
		}
	    public:
		RangePolys()
		    : m_range(NULL)
		{}

		RangePolys(Range *range, unsigned degree)
		    : m_range(range)
		    , m_degree(degree)
		{
		    segment_it_t it=m_range->m_segments.begin();
		    while(it!=m_range->m_segments.end()){
			update_segment_minimax(it);
			++it;
		    }
		}

		void dump(FILE *dst)
		{
		    m_range->dump(dst);
		}
		
		void split_to_error(double error)
		{
		    mpfr_t tol;
		    mpfr_init2(tol, getToolPrecision());
		    mpfr_set_d(tol, error, MPFR_RNDN);
		
		    unsigned i=0;
		    segment_it_t it=m_range->m_segments.begin();
		    while(it!=m_range->m_segments.end()){
			if(::flopoco::verbose>=INFO){
			    std::cerr<<"Checking poly "<<i<<" of "<<m_range->m_segments.size()<<"\n";
			}
			update_segment_minimax(it);

			MPFRVec serror=boost::any_cast<MPFRVec>(it->properties["minimax_error"]);

			if(::flopoco::verbose>=INFO){
			    mpfr_fprintf(stderr, "  error = %Re, tol = %Re\n", serror[0], tol);
			}
			if(mpfr_nan_p(it->rangeStart)){
			    ++it;
			    ++i;
			}else{
			    assert(!mpfr_nan_p(serror[0]));


			    if( mpfr_greaterequal_p(serror[0],tol)){
				if(::flopoco::verbose>=DEBUG){
				    mpfr_fprintf(stderr, "Splitting [%Rf,%Rf]",
						 it->domainStart, it->domainFinish
						 );
				    mpfr_fprintf(stderr, " error=%Rg\n", serror[0]);
				}
				m_range->split_segment(it);
				update_segment_minimax(it);
			    }			
			    assert(it->has_property("minimax_error"));
			
			    serror=boost::any_cast<MPFRVec>(it->properties["minimax_error"]);
			    if(mpfr_less_p(serror[0],tol)){
				++it;
				++i;
			    }
			}
		    }
		    mpfr_clear(tol);
		}
	
		// Find fixed-point polynomial coefficients such that the result is faithful (?) on rangeWF+guard bits
		segment_it_t calc_faithful_fixed_point(segment_it_t curr, int guard=0)
		{
		    std::string pnCoeffs=boost::str(boost::format("minimax_fixed%1%_coeffs")%guard);
		
		    if(curr->has_property(pnCoeffs))
			return curr;
		
		    update_segment_minimax(curr);
		
		    mpfr_t tol;
		    mpfr_init2(tol, getToolPrecision());
		    mpfr_set_d(tol, pow(2.0, -m_range->m_rangeWF-2), MPFR_RNDN);
		
		    MPFRVec serror=boost::any_cast<MPFRVec>(curr->properties["minimax_error"]);
		    if(mpfr_greater_p(serror[0],tol)){
			mpfr_fprintf(stderr, " error=%Rg\n", serror[0]);
			mpfr_clear(tol);

			// Ok, let's split here
			m_range->split_segment(curr);
			curr=calc_faithful_fixed_point(curr, guard);
			++curr;
			calc_faithful_fixed_point(curr, guard);
			return curr;
			//			throw std::runtime_error("calc_faithful_fixed_point - poly is not accurate enough to be faithful (?).");
		    }
		
		    // Ok, let's try to do this properly. Heuristic from ASAP 2010 is that where:
		    // k = leading zeros. Equivalent to index bits in a table-poly over [0,1). Here k=1, as we always have a leading zero due to range over [0,1)
		    // j = Coefficient degree
		    // p = target precision in bits (e.g. precision is 2^-p)
		    // then we need
		    //   lsb(c[i]) = -p + k*j
		    // In our case that simplifies to
		    //   lsb(c[i]) = -p + j
		    std::vector<int> lsbs(m_degree+1), fmts(m_degree+1);
		    for(unsigned j=0;j<=m_degree;j++){
			lsbs[j]= -(m_range->m_rangeWF+guard+1) + j -1;	// TODO : I seem to need one more LSB than the paper says, otherwise PolynomialEvaluator doesn't converge
			fmts[j]=-lsbs[j]; // format is negative of lsb
		    }
		    if(::flopoco::verbose>3)
			fprintf(stderr, " doing fpminimax...");
		    auto currMinimax= boost::any_cast<sollya_node_ptr_t>(curr->properties["minimax"]);
		    if(::flopoco::verbose>3){
			if(currMinimax)
			    fprintf(stderr, " (with init minimax)");
		    }
		    sollya_node_t fp=curr->fpminimax(fmts, currMinimax.get() );
		
		    MPFRVec error(1, getToolPrecision());

		    MPFRVec coeffs;
		
		    if(curr->isRangeConstant || !mpfr_regular_p(curr->rangeStart)){
			// We're ok, there is no polynomial
			coeffs=MPFRVec(m_degree+1, getToolPrecision());
			if(mpfr_regular_p(curr->rangeStart)){
			    mpfr_set(coeffs[0], curr->rangeStartFrac, MPFR_RNDN);
			}else{
			    mpfr_set(coeffs[0], curr->rangeStart, MPFR_RNDN);
			}
			for(int i=1;i<=(int)m_degree;i++){
			    mpfr_set_zero(coeffs[i],+1);
			}
			mpfr_set_zero(error[0],+1);
		    }else{

			fprintf(stderr, "infNorm...");
			sollya_node_t func=curr->get_scaled_flat_function();
			infNorm(error[0], func, fp, curr->domainStartScaled, curr->domainFinishScaled, 71);
			free_memory(func);

			fprintf(stderr, "\n");

			mpfr_abs(error[0],error[0],MPFR_RNDN);
			coeffs=getPolyCoeffs(fp, m_degree);
			assert((int)coeffs.size()==(int)m_degree+1);
			free_memory(fp);
		    }
		
		    mpfr_set_d(tol, pow(2.0, -m_range->m_rangeWF-1), MPFR_RNDN);
		    if(mpfr_greater_p(error[0],tol)){
			mpfr_fprintf(stderr, " didn't meet error: wanted=%Rg, got=%Rg\n", tol, error[0]);
			mpfr_clear(tol);

			// Ok, let's split here
			m_range->split_segment(curr);
			curr=calc_faithful_fixed_point(curr, guard);
			++curr;
			calc_faithful_fixed_point(curr, guard);
			return curr;
		    }else{
		
		std::string pnLsbs=boost::str(boost::format("minimax_fixed%1%_lsbs")%guard);
			std::string pnError=boost::str(boost::format("minimax_fixed%1%_error")%guard);
		
			curr->properties[pnCoeffs]=coeffs;
			curr->properties[pnLsbs]=lsbs;
			curr->properties[pnError]=error;
		
			mpfr_clear(tol);

			return curr;
		    }
		}
		
	
		void calc_faithful_fixed_point(int guard=0)
		{
		    // TODO : We split it so that the minimax is twice as accurate as needed. This means we use too
		    // many segments, and should be optimised. Yes, I suck at maths.
		    split_to_error(pow(2.0, -m_range->m_rangeWF-2));
		
		    unsigned i=0;
		    segment_it_t curr=m_range->m_segments.begin();
		    while(curr!=m_range->m_segments.end()){
			if(::flopoco::verbose>3)
			    mpfr_fprintf(stderr, "  calc_faithful_fixed_point([%Re,%Re]) %d of %d\n", curr->domainStart, curr->domainFinish, i, m_range->m_segments.size());


			curr=calc_faithful_fixed_point(curr, guard);
			++curr;
			++i;
		    }
		}
	
		MPFRVec m_concretePartition;
		std::vector<int> m_concreteExp;
		std::vector<int> m_concreteCoeffMsbs, m_concreteCoeffLsbs;
		MPFRVec m_concreteApproxError;	// worst polynomial error
	
		// Signed number with bits:  sign,  msb, msb-1, ... lsb+1, lsb
		// s, -1, -2, -3, -4 -> width=5=1+(-1- -4)=1+(msb-lsb)
		mpz_class ToTwosComplement(mpz_class x, int msb, int lsb)
		{
		    int w=2+msb-lsb;
		    //mpfr_fprintf(stderr, "   x=%Zd=%Zx, msb=%d, lsb=%d, width=%d\n", x.get_mpz_t(), x.get_mpz_t(), msb, lsb, w);
		
		    mpz_class ret;
		    if(x>=0){
			ret=x;
		    }else{
			for(int i=w-1;i>=0;i--){
			    ret=(ret<<1) + mpz_tstbit(x.get_mpz_t(), i);	// manually read out two's complements bits. Heh.
			}
		    }
		    //mpfr_fprintf(stderr, "  ret=%Zd = %Zx\n", ret.get_mpz_t(), ret.get_mpz_t());
		    mpz_class tmp;
		    mpz_ui_pow_ui(tmp.get_mpz_t(), 2, w);
		    assert(ret>=0);
		    assert(ret<tmp);
		    return ret;
		}
	
		void build_concrete(int guard=0)
		{
		    m_concreteApproxError=MPFRVec(1,getToolPrecision());
		    mpfr_set_d(m_concreteApproxError[0], 0.0, MPFR_RNDN);
		
		    mpfr_t tmp;
		    mpfr_init2(tmp, getToolPrecision());
		
		    std::string pnCoeffs=boost::str(boost::format("minimax_fixed%1%_coeffs")%guard);
		    std::string pnLsbs=boost::str(boost::format("minimax_fixed%1%_lsbs")%guard);
		    std::string pnError=boost::str(boost::format("minimax_fixed%1%_error")%guard);
		
		    m_concretePartition.resize(m_range->m_segments.size(), m_range->m_domainWF+1);
		    m_concreteExp.resize(0);
		
		    // These are chosen to coincide with setting for zero coefficients from PolynomialEvaluator.cpp
		    //	It is just a two-bit signed format of 0s.x00
		    std::vector<int> coeffMsbs(m_degree+1, 0);
		    std::vector<int> coeffLsbs(m_degree+1, -1);		
		
		    // First walk over all the coefficients, and find the min lsb and max msb of each coefficient
		    segment_it_t curr=m_range->m_segments.begin();
		    for(unsigned i=0;i<m_range->m_segments.size();i++){
			assert(curr->has_property(pnCoeffs));
			MPFRVec coeffs=boost::any_cast<MPFRVec>(curr->properties[pnCoeffs]);
			assert((int)coeffs.size()==(int)m_degree+1);
			std::vector<int> lsbs=boost::any_cast<std::vector<int> >(curr->properties[pnLsbs]);
			MPFRVec error=boost::any_cast<MPFRVec>(curr->properties[pnError]);
			
			if(::flopoco::verbose>=DEBUG){
			    mpfr_fprintf(stderr, "segment[%d] : [%Re,%Re] -> [%Re,%Re]\n", i, curr->domainStart, curr->domainFinish, curr->rangeStart, curr->rangeFinish);
			}

			if(i>0){
			    mpfr_t xx;
			    mpfr_init2(xx, m_range->m_domainWF+1);
			    mpfr_set(xx, curr->domainStart, MPFR_RNDN);
			    mpfr_nextbelow(xx);

			    auto prev=curr;
			    --prev;

			    if(!(mpfr_zero_p(prev->domainFinish) || mpfr_zero_p(curr->domainStart))){

				if(!mpfr_equal_p(prev->domainFinish,xx)){
				    mpfr_fprintf(stderr, "Non contiguous: prev(start[%d])=%Re, finish[%d]=%Re\n", i, xx, i-1, prev->domainFinish);
				    exit(1);
				}
			    }

			    mpfr_clear(xx);
			}

			for(unsigned j=0;j<=m_degree;j++){
			    coeffLsbs[j]=std::min(coeffLsbs[j], lsbs[j]);
			    int msb=std::max(lsbs[j], (int)mpfr_get_exp(coeffs[j])-1);// 0.5*2^0=0.5 -> mpfr_get_exp(0.5)==0,  0.5*2^1=1 -> mpfr_get_exp(1)==1
			    coeffMsbs[j]=std::max(coeffMsbs[j], msb);	
			    if(::flopoco::verbose>=DEBUG){
				mpfr_fprintf(stderr, "  x^%d : format[%d:%d] = %Rg\n", j, msb, lsbs[j], coeffs[j]);
			    }
			}
			if(::flopoco::verbose>=DEBUG){
			    mpfr_fprintf(stderr, "\n");
			}
			
			mpfr_set(m_concretePartition[i], curr->domainFinish, MPFR_RNDN);
			m_concreteExp.push_back(mpfr_get_exp(curr->rangeStart));
			
			if(mpfr_greater_p(error[0], m_concreteApproxError[0]))
			    m_concreteApproxError=error;
			
			++curr;
		    }
		
		    if(::flopoco::verbose>=INFO){
			mpfr_fprintf(stderr, "overall coeff formats :\n");
			for(unsigned j=0;j<=m_degree;j++){
			    mpfr_fprintf(stderr, "  format[%d:%d], width=%d\n", coeffMsbs[j], coeffLsbs[j], coeffMsbs[j]-coeffLsbs[j]+1);
			}
			mpfr_fprintf(stderr, "\n");
		    }
		
		    m_concreteCoeffMsbs=coeffMsbs;
		    m_concreteCoeffLsbs=coeffLsbs;
		}
	
		std::vector<mpfr::mpreal> get_polynomial(segment_it_t it, int guard)
		{
		    std::string pnCoeffs=boost::str(boost::format("minimax_fixed%1%_coeffs")%guard);

		    MPFRVec coeffs=boost::any_cast<MPFRVec>(it->properties[pnCoeffs]);
		    std::vector<mpfr::mpreal> res;
		    for(int i=0;i<coeffs.size();i++){
			res.push_back(mpfr::mpreal(0, getToolPrecision()));
			mpfr_set(get_mpfr_ptr(res.back()), coeffs[i], MPFR_RNDN);
		    }
		    return res;
		}

		MPFRVec calculate_polynomial_bounds(segment_it_t it, int guard, unsigned trials=1000)
		{
		    MPFRVec boundsHi(m_degree+1, getToolPrecision());
		    MPFRVec boundsLo(m_degree+1, getToolPrecision());

		    mpfr_t x, step, tmp, acc;
		    mpfr_inits2(getToolPrecision(), x,acc, tmp, step, NULL);
		    
		    mpfr_sub(step, it->domainFinishScaled, it->domainStartScaled, MPFR_RNDN);
		    mpfr_div_si(step, step, trials-1, MPFR_RNDN);

		    std::string pnCoeffs=boost::str(boost::format("minimax_fixed%1%_coeffs")%guard);
		    MPFRVec coeffs=boost::any_cast<MPFRVec>(it->properties[pnCoeffs]);

		    mpfr_set(x, it->domainStartScaled, MPFR_RNDN);
		    for(unsigned i=0;i<trials;i++){
			mpfr_set(acc, coeffs[m_degree], MPFR_RNDN);
			for(int j=m_degree-1;j>=0;j--){
			    mpfr_mul(acc, x, acc, MPFR_RNDN);
			    mpfr_add(acc, coeffs[j], acc, MPFR_RNDN);
			    
			    if(i==0){
				mpfr_set(boundsHi[j], acc, MPFR_RNDN);
				mpfr_set(boundsLo[j], acc, MPFR_RNDN);
			    }else{
				mpfr_max(boundsHi[j], boundsHi[j], acc, MPFR_RNDN);
				mpfr_min(boundsLo[j], boundsLo[j], acc, MPFR_RNDN);
			    }
			}
		    }

		    mpfr_clears(x,acc,step,tmp,NULL);

		    mpfr_fprintf(stderr, "  [%Re,%Re]->[%Re,%Re], [%Re,%Re]->[%Re,%Re], off=%Re, scale=%d, [%Re,%Re], rC=%s ", it->domainStart, it->domainFinish, it->rangeStart, it->rangeFinish, it->domainStartFrac, it->domainFinishFrac, it->rangeStartFrac, it->rangeFinishFrac, it->domainOffset, it->domainScale, it->domainStartScaled, it->domainFinishScaled,  it->isRangeConstant?"Y":"N");
		    fprintf(stderr, "\n");
		    for(int j=0;j<=m_degree;j++){
			mpfr_fprintf(stderr, " a%d:%Re;", j, coeffs[j]);
		    }

		    for(int j=0;j<=m_degree;j++){
			mpfr_abs(boundsHi[j],boundsHi[j],MPFR_RNDN);
			mpfr_abs(boundsLo[j],boundsLo[j],MPFR_RNDN);
			mpfr_max(boundsHi[j],boundsHi[j],boundsLo[j],MPFR_RNDN);
			fprintf(stderr, " %3d", mpfr_get_exp(boundsHi[j]));
			
		    }
		    fprintf(stderr, "\n");
		    return boundsHi;
		}

		MPFRVec calculate_polynomial_bounds(int guard, unsigned trials=1000)
		{
		    MPFRVec res(m_degree+1, getToolPrecision());
		    for(int i=0;i<=m_degree;i++){
			mpfr_set_inf(res[i],-1);
		    }

		    segment_it_t it=m_range->m_segments.begin();
		    while(it!=m_range->m_segments.end()){
			if(mpfr_regular_p(it->rangeStart)){
			    auto b=calculate_polynomial_bounds(it, guard, trials);
			    for(int j=0;j<=m_degree;j++){
				mpfr_max(res[j],res[j],b[j],MPFR_RNDN);
			    }
			}
			++it;
		    }

		    return res;
		}
		
		// Convert the polynomial coefficients into the actual data-values
		/* The format is  |prefix|an|...|a0|, where prefix is the floating point exponent+sign that goes on the front
		 */
		std::vector<mpz_class> build_ram_contents(int guard, int wRangeE)
		{
		    auto bounds=calculate_polynomial_bounds(guard);
		    for(unsigned i=0;i<=m_degree;i++){
			int bits=mpfr_get_exp(bounds[i])+1;
			mpfr_fprintf(stderr, " range[%d] = %Re, bits=%d\n", i, bounds[i], bits); 
		    }

		    std::string pnCoeffs=boost::str(boost::format("minimax_fixed%1%_coeffs")%guard);
		
		    mpfr_t tmp;
		    mpfr_init2(tmp, getToolPrecision());
		
		    int totalCoeffBits=0;
		    for(unsigned i=0;i<=m_degree;i++){
			totalCoeffBits+=(m_concreteCoeffMsbs[i]-m_concreteCoeffLsbs[i]+1)+1;
		    }

		    int totalScaleBits=0;
		    if(m_range->m_isDomainScaled){
			// The maximum shift is of m_domainWF-1
			// The maximum offset is m_domainWF
			totalScaleBits+=ceil(log(m_range->m_domainWF-1)/log(2.0));
			totalScaleBits+=m_range->m_domainWF;
		    }
		
		    if(::flopoco::verbose>=INFO){
			fprintf(stderr, "  Building table.\n");
		    }
		
		    std::vector<mpz_class> contents;
		
		    // Now we know the precise types of the coefficients - they are all signed with
		    // bits in the range coeffMsbs..coeffLsbs. Let's build the table...
		    segment_it_t curr=m_range->m_segments.begin();
		    while(curr!=m_range->m_segments.end()){
			assert(curr->domainScale==0);
			assert(mpfr_zero_p(curr->domainOffset));

			MPFRVec coeffs=boost::any_cast<MPFRVec>(curr->properties[pnCoeffs]);
			
			mpz_class acc, local;
			for(int i=m_degree;i>=0;i--){
			    mpfr_set(tmp, coeffs[i], MPFR_RNDN);
			    mpfr_mul_2si(tmp, tmp, -m_concreteCoeffLsbs[i], MPFR_RNDN);
			    mpfr_get_z(local.get_mpz_t(), tmp, MPFR_RNDN);
				
			    acc=(acc<<(2+m_concreteCoeffMsbs[i]-m_concreteCoeffLsbs[i])) + ToTwosComplement(local, m_concreteCoeffMsbs[i], m_concreteCoeffLsbs[i]);
			}
			
			// Now we tack the exponent on. Now have support for positive, zero, and negative numbers,
			// so it will be "010|(e-2^wE/2)|...
			// or            "001|0000000000|...
			// or            "000|0000000000|...
			// or            "011|(e-2^wE/2)|...
			
			bool isNeg=mpfr_sgn(curr->rangeStart)<0;
			if(mpfr_inf_p(curr->rangeStart)){
			    acc=(mpz_class(isNeg?(4+1):(4+0))<<(wRangeE+totalScaleBits+totalCoeffBits)) + /*exponent*/ 0 + acc;
			}else if(mpfr_nan_p(curr->rangeStart)){
			    acc=(mpz_class(6)<<(wRangeE+totalScaleBits+totalCoeffBits)) + /*exponent*/ 0 + acc;	    
			}else if(mpfr_zero_p(curr->rangeStart)){
			    acc=(mpz_class(isNeg ? 1 : 0)<<(wRangeE+totalScaleBits+totalCoeffBits)) + /*exponent*/ 0 + acc;
			}else{
			    mpz_class prefix=mpfr_get_exp(curr->rangeStart)-2+(1<<(wRangeE-1));
			    if((prefix<0) || ((mpz_class(1)<<wRangeE)<=prefix)){
				mpfr_fprintf(stderr, "  range=[%Re,%Re]\n", curr->rangeStart, curr->rangeFinish);
				std::cerr<<"  exponent="<<mpfr_get_exp(curr->rangeStart)<<", wRangeE="<<wRangeE<<"\n";
				throw std::string("Exponent out of range (increase wRangeE?).");
			    }
			    acc=(mpz_class(isNeg?3:2)<<(wRangeE+totalScaleBits+totalCoeffBits)) + (prefix<<(totalScaleBits+totalCoeffBits)) + acc;
			}
			
			contents.push_back(acc);
			
			if(::flopoco::verbose>=INFO){
			    mpfr_fprintf(stderr, "  table[%d]=0x%Zx\n", contents.size(), acc.get_mpz_t());
			}
			
			++curr;
		    }
		
		    if(::flopoco::verbose>=INFO){
			fprintf(stderr, "Done build concrete.\n");
		    }
		
		    mpfr_clear(tmp);
		
		    return contents;
		}
	
	
	
	
		PolynomialEvaluator *make_polynomial_evaluator(Target *target, map<string, double> inputDelays = map<string, double>())
		{
		    typedef PolynomialEvaluator::format_t format_t;
		
		    std::vector<format_t> coeffs(m_degree+1);
		    for(unsigned i=0;i<=m_degree;i++){
			coeffs[i].isSigned=true;
			coeffs[i].msb=m_concreteCoeffMsbs[i]+1;	// need +1 for sign bit
			coeffs[i].lsb=m_concreteCoeffLsbs[i];
		    }
		
		    format_t inputFormat;
		    inputFormat.isSigned=false;
		    inputFormat.msb=-1;
		    inputFormat.lsb=-m_range->m_domainWF;

		    return PolynomialEvaluator::Create(target, coeffs, inputFormat,
						       -m_range->m_rangeWF, //outputLsb
						       m_concreteApproxError[0], inputDelays);
		} 

		FixedPointPolynomialEvaluator *make_fixed_point_polynomial_evaluator(Target *target, map<string, double> inputDelays = map<string, double>())
		{
		    std::vector<fixed_format_t> coeffs(m_degree+1);
		    for(unsigned i=0;i<=m_degree;i++){
			coeffs[i].isSigned=true;
			coeffs[i].msb=m_concreteCoeffMsbs[i]+1;	// need +1 for sign bit
			coeffs[i].lsb=m_concreteCoeffLsbs[i];
		    }
		
		    fixed_format_t inputFormat;
		    inputFormat.isSigned=false;
		    inputFormat.msb=-1;
		    inputFormat.lsb=-m_range->m_domainWF;
		
		    mpfr::mpreal budget(pow(2.0, -m_range->m_rangeWF-1), getToolPrecision());
		    mpfr_sub(get_mpfr_ptr(budget), get_mpfr_ptr(budget), m_concreteApproxError[0], MPFR_RNDN);

		    if(1){
			return CreateRoundingPolynomialEvaluator(
								 inputFormat,
								 coeffs,
								 -m_range->m_rangeWF, //outputLsb
								 3,
								 budget,
								 target
								 );
		    }else{
			return CreateFixedPointPolynomialEvaluator(
								   inputFormat,
								   coeffs,
								   -m_range->m_rangeWF, //outputLsb
								   budget,
								   target
								   );
		    }
		} 
	
		mpz_class ToPositiveConstant(std::pair<int,int> format, mpfr_t x)
		{
		    if(mpfr_sgn(x)<=0)
			throw std::invalid_argument("ToPositiveConstant - non positive numbers not supported.");
		    if(!mpfr_regular_p(x))
			throw std::invalid_argument("ToPositiveConstant - only regular numbers are supported.");
		
		    std::string s=fp2bin(x, format.first, format.second);
		    mpz_class acc;
		    for(unsigned i=0;i<s.size();i++){
			acc = (acc<<1) + ((s[i]=='1')?1:0);
		    }
		    return acc;
		}
	
		StaticQuantiser *make_static_quantiser(Target *target, int wDomainE, map<string, double> inputDelays = map<string, double>())
		{
		    std::vector<mpz_class> boundaries;
		
		    ComparableFloatType type(wDomainE, m_range->m_domainWF);
		
		    segment_it_t curr=m_range->m_segments.begin();
		    while(curr!=m_range->m_segments.end()){
			boundaries.push_back(type.ToBits(curr->domainStart));
			++curr;
		    }
		    boundaries.push_back((mpz_class(1)<<type.Width())-1);	// Not used by quantiser
		
		    StaticQuantiser *sq=new StaticQuantiser(target, type.Width(), boundaries, inputDelays);
		
		    return sq;
		}
	
		unsigned find_segment(mpfr_t x)
		{
		    if(!m_range->is_in_domain(x))
			throw std::invalid_argument("find_segment - x not in domain.");

		    int i;
		    for(i=0;i<m_concretePartition.size();i++){
			if(mpfr_lessequal_p(x, m_concretePartition[i]))
			    break;
		    }
		    if(i==m_concretePartition.size())
			throw std::invalid_argument("eval_concrete - out of range.");
		    return i;
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
	
		/*
		  void exhaust_concrete()
		  {
		  mpfr_t x, got, want, actual, err, tol;
		
		  mpfr_init2(got, m_range->m_rangeWF);
		  mpfr_init2(want, m_range->m_rangeWF);
		  mpfr_init2(actual, getToolPrecision());
		  mpfr_init2(err, getToolPrecision());
		  mpfr_init2(tol, getToolPrecision());
		  mpfr_init_copy(x, m_range->m_segments.begin()->domainStart);
		
		  mpfr_set_si_2exp(tol, 1, -m_range->m_rangeWF, MPFR_RNDN);
		
		  int start=time(NULL);
		  int delta_print=2;
		  int next_print=start+delta_print;
		
		  while(mpfr_lessequal_p(x, boost::prior(m_range->m_segments.end())->domainFinish)){
		  unsigned seg=eval_concrete(got, x);
		  m_range->m_function.eval(actual, x);
		  mpfr_set(want, actual, MPFR_RNDN);
			
		  mpfr_sub(err, got, actual, MPFR_RNDN);
		  mpfr_div(err, err, actual, MPFR_RNDN);
		  mpfr_abs(err, err, MPFR_RNDN);
			
		  int u=ulps(got, want);
			
		  if(u>1){
		  mpfr_fprintf(stdout, "%u, %Rg, %Rg, %Rg, %Rg, %Rg, %s, %d\n", seg, x, got, want, actual, err, mpfr_greater_p(err, tol) ? "FAIL" : "pass", u);
		  }
			
		  mpfr_nextabove(x);
			
		  if(time(NULL) > next_print){
		  mpfr_fprintf(stderr, "  curr=%Rg, range=[%Rg,%Rg]\n", x, m_range->m_segments.begin()->domainStart, boost::prior(m_range->m_segments.end())->domainFinish);
		  delta_print=std::min(60, delta_print*2);
		  next_print=time(NULL)+delta_print;
		  }
		  }
		
		  mpfr_clears(x, got, want, (mpfr_ptr)0);
		  }
		*/
	    };

	}; // float_approx
    };	// random
}; // flopco

#endif
