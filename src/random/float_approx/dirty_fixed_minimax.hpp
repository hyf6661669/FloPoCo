#ifndef flopoco_random_float_approx_dirty_fixed_minimax_hpp
#define flopoco_random_float_approx_dirty_fixed_minimax_hpp

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
	
	/*! This is a really dumb version of fpminima that just iterates through the
		coefficients from lowest to highest degree, and rounds them. 
		I realise the full version is much better, but it just fails to converge
		too often.
	*/
	sollya_node_t dirtyFPminimax(sollya_node_t f,
		sollya_chain_t monomials,
		sollya_chain_t formats,
		sollya_chain_t points,
		mpfr_t a, mpfr_t b,
		int fp, int absrel,
		sollya_node_t consPart,
		sollya_node_t minimax
	){
		if(absrel!=ABSOLUTESYM)
			throw std::string("dirtyFPMinimax - only absolute error is supported.");
		if(fp!=FIXED)
			throw std::string("dirtyFPMinimax - only fixed-point is supported.");
		if(consPart){
			if(!isConstant(consPart)){
				throw std::string("dirtyFPMinimax - consPart must be truly constant.");
			}
			mpfr_t tmp;
			mpfr_init2(tmp, getToolPrecision());
			evaluateConstantExpression(tmp, consPart, getToolPrecision());
			int isZero=mpfr_zero_p(tmp);
			mpfr_clear(tmp);
			if(!isZero)
				throw std::string("dirtyFPMinimax - constPart must be zero if it exists.");
		}
		if(monomials==0)
			throw std::string("dirtyFPMinimax - monomials must be non-NULL.");
		if(formats==0)
			throw std::string("dirtyFPMinimax - formats must be non-NULL.");
		if(lengthChain(monomials)!=lengthChain(formats))
			throw std::string("dirtyFPMinimax - formats and monomials must be the same length.");
		
		// Parameter checking done, we are assured of boring formats. Now
		// unpack the lsbs of the coefficients
		std::vector<int> lsbs;
		sollya_chain_t mt=monomials, ft=formats;
		while(mt){
			int mm=*(int*)first(mt);
			int ff=-*(int*)first(ft); // Passed in as te number of fractional bits
			if(mm!=(int)lsbs.size())
				throw std::string("dirtyFPMinimax - Monomials must be a contiguous ascending list from 0 up.");
			lsbs.push_back(ff);
			mt=tail(mt);
			ft=tail(ft);

			fprintf(stderr, "  lsbs[%lu]=%d\n", lsbs.size()-1, ff);
		}		
		
		// We'll slowly build them up in here
		mpfr_t *coeffs=(mpfr_t*)malloc(sizeof(mpfr_t)*lsbs.size());
		for(unsigned i=0;i<lsbs.size();i++){
			mpfr_init2(coeffs[i], getToolPrecision());
		}
		
		mpfr_t quality;
		mpfr_init2(quality, getToolPrecision());
		mpfr_set_d(quality, ldexp(1.0, lsbs[0]-1), MPFR_RNDN);
		sollya_node_t weight=makeConstantDouble(1.0);
		
		// Now just iterate over the coefficients, rounding one each time
		for(unsigned i=0;i<lsbs.size();i++){
			sollya_node_t currFunc=0;
			sollya_chain_t currMonom=makeIntPtrChainFromTo(i, lsbs.size()-1);
			if(i==0){
				currFunc=copyTree(f);
			}else{
				sollya_node_t currPoly=makePolynomial(coeffs, i-1);
				currFunc=makeSub(copyTree(f), currPoly);
			}

			fprintf(stderr, "currFunc = ");
			fprintTree(stderr, currFunc);
			fprintf(stderr, "\n");

			
			sollya_node_t solution=remez(currFunc, weight, currMonom, a, b, &quality, getToolPrecision());
			sollya_node_t x=getIthCoefficient(solution, i);
			evaluateConstantExpression(coeffs[i], x, getToolPrecision());

			
			free_memory(x);
			free_memory(solution);
			free_memory(currFunc);
			freeChain(currMonom,freeIntPtr);

			mpfr_fprintf(stderr, " coeff[%d] = %Re (pre-rounding)\n", i, coeffs[i]);
			
			// Now do the rounding
			mpfr_mul_2si(coeffs[i], coeffs[i], -lsbs[i], MPFR_RNDN);
			mpfr_round(coeffs[i], coeffs[i]);
			mpfr_mul_2si(coeffs[i], coeffs[i], lsbs[i], MPFR_RNDN);

			mpfr_fprintf(stderr, " coeff[%d] = %Re (post-rounding)\n", i, coeffs[i]);
		}
		
		sollya_node_t res=makePolynomial(coeffs, lsbs.size()-1);
		
		mpfr_clear(quality);
		free_memory(weight);
		for(unsigned i=0;i<lsbs.size();i++){
			mpfr_clear(coeffs[i]);
		}
		free(coeffs);

		fprintf(stderr, " Original = ");
		fprintTree(stderr, minimax);
		fprintf(stderr, "\n\n dirty = ");
		fprintTree(stderr, res);
		fprintf(stderr, "\n");

		return res;
	}

    /* So the dirty version didn't work (I probably need to
       step the weight function down while moving the decided
       part to the other side). So... I'm just rounding the
       f***ing coefficients. Do not care any more; we will
       probably waste some bits, it may mean a bigger guard
       band further up (or even failures for higher degrees),
       but life is too short.
    */
	sollya_node_t dirtierFPminimax(sollya_node_t f,
		sollya_chain_t monomials,
		sollya_chain_t formats,
		sollya_chain_t points,
		mpfr_t a, mpfr_t b,
		int fp, int absrel,
		sollya_node_t consPart,
		sollya_node_t minimax
	){
		if(absrel!=ABSOLUTESYM)
			throw std::string("dirtyFPMinimax - only absolute error is supported.");
		if(fp!=FIXED)
			throw std::string("dirtyFPMinimax - only fixed-point is supported.");
		if(consPart){
			if(!isConstant(consPart)){
				throw std::string("dirtyFPMinimax - consPart must be truly constant.");
			}
			mpfr_t tmp;
			mpfr_init2(tmp, getToolPrecision());
			evaluateConstantExpression(tmp, consPart, getToolPrecision());
			int isZero=mpfr_zero_p(tmp);
			mpfr_clear(tmp);
			if(!isZero)
				throw std::string("dirtyFPMinimax - constPart must be zero if it exists.");
		}
		if(monomials==0)
			throw std::string("dirtyFPMinimax - monomials must be non-NULL.");
		if(formats==0)
			throw std::string("dirtyFPMinimax - formats must be non-NULL.");
		if(lengthChain(monomials)!=lengthChain(formats))
			throw std::string("dirtyFPMinimax - formats and monomials must be the same length.");
		
		// Parameter checking done, we are assured of boring formats. Now
		// unpack the lsbs of the coefficients
		std::vector<int> lsbs;
		sollya_chain_t mt=monomials, ft=formats;
		while(mt){
			int mm=*(int*)first(mt);
			int ff=-*(int*)first(ft); // Passed in as te number of fractional bits
			if(mm!=(int)lsbs.size())
				throw std::string("dirtyFPMinimax - Monomials must be a contiguous ascending list from 0 up.");
			lsbs.push_back(ff);
			mt=tail(mt);
			ft=tail(ft);

			//fprintf(stderr, "  lsbs[%d]=%d\n", lsbs.size()-1, ff);
		}		

		mpfr_t quality;
		mpfr_init2(quality, getToolPrecision());
		mpfr_set_d(quality, ldexp(1.0, lsbs[0]-1), MPFR_RNDN);
		sollya_node_t weight=makeConstantDouble(1.0);

		sollya_node_t solution=minimax;
		if(!minimax)
		    solution=remez(f, weight, monomials, a, b, &quality, getToolPrecision());
	
		mpfr_t *coeffs=(mpfr_t*)malloc(sizeof(mpfr_t)*lsbs.size());
		for(unsigned i=0;i<lsbs.size();i++){
			mpfr_init2(coeffs[i], getToolPrecision());

			sollya_node_t x=getIthCoefficient(solution, i);
			evaluateConstantExpression(coeffs[i], x, getToolPrecision());
			
			free_memory(x);		

			mpfr_mul_2si(coeffs[i], coeffs[i], -lsbs[i], MPFR_RNDN);
			mpfr_round(coeffs[i], coeffs[i]);
			mpfr_mul_2si(coeffs[i], coeffs[i], lsbs[i], MPFR_RNDN);

			// mpfr_fprintf(stderr, " coeff[%d] = %Re (post-rounding)\n", i, coeffs[i]);
		}
		
		sollya_node_t res=makePolynomial(coeffs, lsbs.size()-1);
		
		mpfr_clear(quality);
		free_memory(weight);
		for(unsigned i=0;i<lsbs.size();i++){
			mpfr_clear(coeffs[i]);
		}
		free(coeffs);

		if(!minimax)
		    free_memory(solution);


		return res;
	}

}; // float_approx
}; // random
}; // flopoco

#endif
