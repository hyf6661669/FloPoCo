

class DumbPolynomial
{
private:
	int m_degree;	
	int m_lsbX;	// precision of input x is 2^m_lsbX
	int m_msbX;	// weight of lsb for x
	bool m_isSignedX;	// If true, then m_msbX is the sign bit

	int m_lsbY;	// precision of output y
	mpfr::mpreal m_epsPoly;	// Maximum absolute error of the raw polynomial (i.e. gauranteed norm out of sollya)
	
	std::vector<int> m_lsbCoeffs;	// lsb weight of each coefficient

	std::vector<int> m_lsbStages;	// lsb after each multiply-accumulate and round.
	std::vector<int> m_msbStages;	// msb after each multiply-accumulate and round (if signed, then this will be the sign bit)
	std::vector<bool> m_isSignedStages;	// Whether the output of a given stage is signed or unsigned
	std::vector<mpfr::mpreal> m_epsStages;	// Maximum cumulative rounding error present at each stage

	/*
		We are trying to determine some y_r which is a rounded version of y_t with lsbRes, and an error of
		less than epsULP=2^-lsbRes compared to the true y
		
		The input is some known polynomial, with absolute error guaranteed to be less than epsP < epsULP/2
		y_t=f(x)
		y_p=a0 + x*a1 + x^2*a2 + ... x^d*ad
		| y_p-y_t | < epsP

		We are now trying to find some rounding scheme that results in an error epsR < epsULP/2, so that
		when we finally round, causing an additional error of epsULP/2, then the total error is <epsULP.
		For the final level of the polynomial we have to round to lsbULP, but there could be internal
		rounding errors as well (though we could always keep all intermediate bits, which would mean
		no error in evaluation of the polynomial).

		y_r = a0 + x * acc1 + epsP + epsULP/2
		|y_r - y_t | < epsP+epsULP/2 < epsULP
			
		acc1 = round(a1 + x * acc2, prec1)
			= a1 + x*acc2 + eps1, where eps1=2^-prec1
		acc2 = round(a2 + x * acc3, prec2)
			= a2 + x*acc3 + eps2, where eps2=2^-prec2
		acc3 = a3

		acc0 = a0 + x*(a1 + x*acc2 + eps1)
			= a0 + x*a1 + x^2 *acc2 + x*eps1
			= a0 + x*a1 + x^2 *(a2 + x*acc3 + eps2) + x*eps1
			= a0 + x*a1 + x^2 * a2 + x^3 *acc3 + x^2 * eps2 + x*eps1
			= a0 + x*a1 + x^2 * a2 + x^3 *a3 + x^2 * eps2 + x*eps1
			
		So we need:
		|y_r - y_t| < epsd*|x^d| + ... + eps2*|x^2| + eps1*|x| + epsP + epsULP/2 < epsULP
		which requires:
		epsd*|x^d| + ... + eps2*|x^2| + eps1*|x| < epsULP-epsULP/2 -epsP = epsULP/2 - epsP
	*/

	mpfr::mpreal MaxAbsError(const mpfr::mpreal &maxX, const mpfr::mpreal &epsP, const std::vector<int> &prec)
	{
		mpfr::mpreal curr=epsP
		for(int i=0;i<prec.size();i++){
			if(prec[i]!=INT_MAX){
				curr += pow(maxX, i) * pow(2.0, -prec-1);
			}
		}
		return curr;
	}
	
	void ChoosePrec(
		int lsbsTarget,
		const mpfr::mpreal &maxX, const mpfr::mpreal &epsP, 
		const std::vector<int> &lsbX,
		const std::vector<int> &lsbsCoeffs
	){
		double epsULP=pow(2.0, -lsbsTarget);
		
		mpfr::mpreal epsBudget=epsULP/2 - epsP;
		if(epsBudget < 0)
			throw std::string("ChosePrec - Not enough precision in polynomial to meet target.");
		
		// Set up guaranteed-good rounding - will always work, as there is no rounding at all (!)
		std::vector<int> lsbs(degree);
		int lsbsPrec=lsbsCoeffs[degree];
		for(int i=degree-1;i>=0;i--){
			int lsbsCurr=std::max(lsbsPrev + lsbX, lsbsCoeffs[i]);
			lsbs[i]=lsbsCurr;
		}
		
		m_epsStages[degree]=0;
		
		// We're just going to walk along from high degree to low degree
		// accumulating epsBudget/(degree) or less per stage.
		for(int i=degree-1;i>1;i--){
			assert(epsBudget >= 0);	// if invariant is true, then we can always get a solution
			
			mpfr::mpreal epsStageMax=epsBudget / i;
			mpfr::mpreal epsCurr=0;	// Always zero for this stage using initial rounding
			mpfr::mpreal powX=pow(maxX, i);
			while(true){
				int lsbsAlt=lsbs[i]+1;
				mpfr::mpreal epsAlt=ldexp(powX, -lsbsAlt-1);
				if(epsAlt < epsStageMax){
					epsCurr=epsAlt;
					lsbs[i]=lsbsNew;
				}else{
					break;
				}
			}
			m_epsStages[i]=epsCurr + m_epsStages[i+1];
			epsBudget -= epsCurr;
		}
		assert(epsBudget >= 0);
		
		lsbs[0]=m_lsbX;	// This is the final epsULP/2 rounding
		m_epsStages[0]=epsP + epsULP/2 + m_epsStages[1];
		
		return prec;
	}
	
	void AddPolynomialToRanges(const std::vector<mpfr::mpreal> &coeffs)
	{
		assert(coeffs.size()==m_degree+1);
		
		sollya_node_t poly = makeConstant(get_mpfr_ptr(coeffs[m_degree]));
		for(unsigned i=degree-1;i>=0;i--){			
			poly = makeAdd(makeConstant(get_mpfr_ptr(coeffs[m_degree])), makeMul(makeVariable(), chain));
			
			
			sollya_node_t chainUp = makeAdd(copyTree(chain), makeConstant(get_mpfr_ptr(m_epsStages[i])));
			
			free_memory(chainUp);
			
			sollya_node_t chainDown = makeNeg(makeSub(negateNode(copyTree(chain)), makeConstant(m_epsStages)));
			m_minValStages min=-supnorm(chainDown);
		}
	}

	mpfr::mpreal Evaluate(mpfr::mpreal x)
	{
		mpfr::mpreal acc=m_coeffsRaw[m_d];
		for(int i=m_d-1;i>=0;i--){
			mpfr::mpreal mm=acc*x;
			acc = round_to_lsb(acc + mm, m_lsbStages[i]);
		}
		return acc;
	}
};
