#ifndef flopoco_random_utils_poly_min_max_hpp
#define flopoco_random_utils_poly_min_max_hpp

#include "sollya.h"
#include "mpreal/mpreal.h"
#include <vector>
#include <assert.h>
#include <stdlib.h>

namespace flopoco{
namespace random{

std::pair<mpfr::mpreal,mpfr::mpreal> PolyMinMax(const std::vector<mpfr::mpreal> &coeffs, mpfr::mpreal minX, mpfr::mpreal maxX)
{
	mpfr_prec_t prec=std::max(minX.get_prec(), maxX.get_prec());
	for(unsigned  i=0;i<coeffs.size();i++){
		prec=std::max(prec, coeffs[i].get_prec());
	}
	
	int degree=coeffs.size()-1;
	
	if(coeffs.size()>1 && coeffs.back()==0){
		std::vector<mpfr::mpreal> tmp(coeffs.begin(), coeffs.end()-1);
		return PolyMinMax(tmp, minX, maxX);
	}
	
	assert((coeffs.size()==1) || (coeffs.back()!=0));
	
	if(coeffs.size()==0){
		throw std::string("PolyMinMax - Can't have zero degree polynomial.");
		
	}else if(coeffs.size()==1){
		return std::make_pair(coeffs[0],coeffs[0]);
		
	}else if(coeffs.size()==2){
		mpfr::mpreal start=coeffs[0]+coeffs[1]*minX, finish=coeffs[0]+coeffs[1]*maxX;
		return std::make_pair(std::min(start, finish), std::max(start, finish));
		
	}else if(coeffs.size()==3){
		mpfr::mpreal start=coeffs[0]+minX*(coeffs[1]+minX*coeffs[2]);
		mpfr::mpreal finish=coeffs[0]+maxX*(coeffs[1]+maxX*coeffs[2]);
		
		// We have p(x) = c0+x*c1+x^2*c2. If there is a maxima or minima in the centre then it will
		// be where the derivative equals zero:
		//  0 = p'(x) = c1 + 2*x*c2
		//	 - c1 = 2*x*c2
		//   x=-c1 / (2*c2)
		mpfr::mpreal low=std::min(start, finish), high=std::max(start, finish);
		
		mpfr::mpreal rx=-coeffs[1] / (2*coeffs[2]);
		if((rx>minX) && (rx<maxX)){
			mpfr::mpreal y=coeffs[0]+rx*(coeffs[1]+rx*coeffs[2]);
			low=std::min(low, y);
			high=std::max(high, y);
		}
		return std::make_pair(low, high);
		
	}else if(coeffs.size()==4){
		mpfr::mpreal start=coeffs[0]+minX*(coeffs[1]+minX*(coeffs[2]+minX*coeffs[3]));
		mpfr::mpreal finish=coeffs[0]+maxX*(coeffs[1]+maxX*(coeffs[2]+maxX*coeffs[3]));
		
		mpfr::mpreal low=std::min(start, finish);
		mpfr::mpreal high=std::max(start, finish);
		
		// We need to consider solutions to 0=a1+a2*2*x+3*a3*x^2
		mpfr::mpreal a=3*coeffs[3], b=2*coeffs[2], c=coeffs[1];
		mpfr::mpreal disc=b*b - 4*a*c;
		if(disc>=0){			
			mpfr::mpreal tmp=(-b+sqrt(disc))/(2*a);
			if((minX < tmp) && (tmp < maxX)){
				mpfr::mpreal y = coeffs[0]+tmp*(coeffs[1]+tmp*(coeffs[2]+tmp*coeffs[3]));
				low=std::min(low, y);
				high=std::max(high, y);
			}
			
			tmp=(-b-sqrt(disc))/(2*a);
			if((minX < tmp) && (tmp < maxX)){
				mpfr::mpreal y = coeffs[0]+tmp*(coeffs[1]+tmp*(coeffs[2]+tmp*coeffs[3]));
				low=std::min(low, y);
				high=std::max(high, y);
			}
		}
		return std::make_pair(low, high);
	
	}else{
		mpfr::mpreal start=coeffs.back(), finish=coeffs.back();
		sollya_node_t poly=makeConstant(mpfr::mpreal(coeffs[degree]).mpfr_ptr());
		for(int i=degree-1;i>=0;i--){
			start=coeffs[i] + minX*start;
			finish=coeffs[i] + maxX*finish;
			poly=makeAdd(makeConstant(mpfr::mpreal(coeffs[i]).mpfr_ptr()), makeMul(makeVariable(), poly));
		}
		sollya_node_t diff=differentiate(poly);
		
		mpfr::mpreal low=std::min(start, finish), high=std::max(start, finish);
		
		sollya_range_t range;
		// TODO : Is this right?
		mpfr_t rawMinX, rawMaxX;
		mpfr_inits2(prec, rawMinX, rawMaxX, (mpfr_ptr)0);
		mpfr_set(rawMinX, miget_mpfr_ptr(nX), MPFR_RNDN);
		mpfr_set(rawMaxX, maget_mpfr_ptr(xX), MPFR_RNDN);
		range.a=&rawMinX;
		range.b=&rawMaxX;
		sollya_chain_t zeros=fpFindZerosFunction(diff, range, prec);
		sollya_chain_t curr=zeros;
		//std::cerr<<"diff=\n";
		//printTree(diff);
		//std::cerr<<"\n";
		while(curr!=0){
			mpfr::mpreal x(*(const mpfr_t*)curr->value);
			//std::cerr<<"root="<<x<<"\n";
			if((minX < x) && (x<maxX)){
				mpfr::mpreal acc=coeffs[degree];
				for(int i=degree-1;i>=0;i--){
					acc=coeffs[i] + acc*x;
				}
				low=std::min(low, acc);
				high=std::max(high, acc);
			}
			
			curr=curr->next;
		}
		
		struct ff{
			static void freeMpfrPtr(void *p)
			{
				mpfr_clear(*(mpfr_t*)p);
				free(p);
			}
		};
		freeChain(zeros, ff::freeMpfrPtr);
		
		mpfr_clears(rawMinX, rawMaxX, (mpfr_ptr)0);
		free_memory(diff);
		free_memory(poly);
		
		return std::make_pair(low, high);
	}
	
}

}; // random
}; // flopoco

#endif
