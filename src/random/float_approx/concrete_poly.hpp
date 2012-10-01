#ifndef flopoco_random_float_approx_concrete_poly_hpp
#define flopoco_random_float_approx_concrete_poly_hpp

#include <stdio.h>
#include <mpfr.h>
#include <vector>
#include <limits.h>
#include "../../sollya.h"

namespace flopoco
{
namespace random
{
namespace float_approx
{

class MPFRVec
{
private:
	int len, allocLen;
	mpfr_t *values;
public:
	MPFRVec()
		: len(0)
		, allocLen(0)
		, values(NULL)
	{}

	MPFRVec(int size, int prec)
		: len(size)
		, allocLen(size)
		, values((mpfr_t*)malloc(sizeof(mpfr_t)*len))
	{
		for(int i=0;i<size;i++){
			mpfr_init2(values[i], prec);
		}
	}
	
	MPFRVec(const MPFRVec &o)
		: len(o.len)
		, allocLen(o.len)
		, values((mpfr_t*)malloc(sizeof(mpfr_t)*len))
	{
		for(int i=0;i<len;i++){
			mpfr_init_copy(values[i], o.values[i]);
		}
	}
	
	MPFRVec &operator=(const MPFRVec &o)
	{
		if(this!=&o){
			clear();
			len=o.len;
			allocLen=o.len;
			values=(mpfr_t*)realloc(values, sizeof(mpfr_t)*len);
			for(int i=0;i<len;i++){
				mpfr_init_copy(values[i], o.values[i]);
			}
		}
		return *this;
	}
	
	void resize(int newLen, int prec)
	{
		if(newLen<len){
			while(newLen<len){
				mpfr_clear(values[len-1]);
				--len;
			}
		}
		if(newLen>len){
			if(newLen > allocLen){
				values=(mpfr_t*)realloc(values, sizeof(mpfr_t)*newLen);
				allocLen=newLen;
			}
			while(newLen>len){
				mpfr_init2(values[len], prec);
				++len;
			}
		}
	}
	
	void clear()
	{
		resize(0, 0);
	}
	
	void push_back(mpfr_t x)
	{
		if(allocLen <= len){
			allocLen=std::max(8, allocLen*2);
			values=(mpfr_t*)realloc(values, sizeof(mpfr_t)*allocLen);
		}
		assert(len < allocLen);
		mpfr_init_copy(values[len], x);
		++len;
	}
	
	const mpfr_t *begin() const
	{ return values; }
	
	const mpfr_t *end() const
	{ return values+len; }
	
	~MPFRVec()
	{
		clear();
		if(values){
			free(values);
			values=0;
			allocLen=0;
		}
	}
	
	int size() const
	{ return len; }
	
	mpfr_t &operator[](int i)
	{
		assert(i>=0);
		assert(i<len);
		return values[i];
	}
	
	const mpfr_t &operator[](int i) const
	{
		assert(i>=0);
		assert(i<len);
		return values[i];
	}
};

// Find maxima/minima of polynomial over range to working precision
// The returned set will always include the largest extrema, but may contain
// rubbish as well.
MPFRVec find_polynomial_extremal_points(MPFRVec coeffs, MPFRVec range)
{
	assert(coeffs.size()>0);
	assert(range.size()==2);
	
	unsigned degree=coeffs.size()-1;
	MPFRVec res=range;	// include endpoints by default
	
	if(degree<=1){
		// for constant and linear just the endpoints are fine
	}else{
		// Get coefficients for derivative, [0,degree)
		MPFRVec dc(degree, getToolPrecision());
		for(unsigned i=0;i<degree;i++){
			mpfr_mul_d(dc[i], coeffs[i+1], i+1, MPFR_RNDN);
		}
		
		// Create polynomial for derivative, and find zeros (extrema)
		sollya_node_t deriv=makePolynomial(&dc[0], degree-1);
		sollya_range_t rr={&range[0], &range[1]};
		sollya_chain_t zeros=fpFindZerosFunction(deriv, rr, getToolPrecision());
		free_memory(deriv);
		
		sollya_chain_t tmp=zeros;
		while (tmp != NULL) {
			res.push_back( *(mpfr_t*)(tmp->value) );
			tmp=tmp->next;
		}
		freeChain(zeros, freeMpfrPtr);
	}
	return res;
}

MPFRVec find_polynomial_and_subpolynomial_extremal_points(MPFRVec coeffs, MPFRVec range)
{
	unsigned degree=coeffs.size()-1;
	MPFRVec res;
	
	if(degree<=1){
		res=range;
	}else{
		// iterate over sub-polynomials
		for(unsigned i=2;i<=degree;i++){
			MPFRVec part(i+1, getToolPrecision());
			for(unsigned j=0;j<=i;j++){
				/* degree=4, i=2,  0->2, 1->3, 2->4
					degree=4, i=3, 0->1, 2->3, ... 3->4
					degree=4, i=4, 0->0    4->4 */
				mpfr_set(part[j], coeffs[j+(degree-i)], MPFR_RNDN);
			}
			
			MPFRVec tmp=find_polynomial_extremal_points(part, range);
			assert(tmp.size()>=2);
			
			for(int j=0;j<tmp.size();j++){
				res.push_back(tmp[j]);
			}
		}
	}
	
	assert(res.size()>0);		
	
	return res;
}

// Completely specified fixed-point horner-form polynomial
class ConcretePoly
{
	unsigned degree;
	MPFRVec coeffs;
	MPFRVec range;	// 2 element vector with [start,finish] range
	
	int inputFracBits, inputMsb;
	
	std::vector<int> coeffFracBits;		// [0..degree] fractional-precision of each fixed-point constant
	std::vector<int> coeffMsbs;		// [0..degree] msb of each fixed-point signed constant
	
	std::vector<int> calcFracBits;		// [0..degree) number of fractional bits kept after each addition (i.e. end of stage)
	std::vector<int> calcMsbs;	// empty, or [0..degree) list of msbs+1 of output of calculation stages
	

	void update_intermediates(mpfr_t x, mpfr_t rx, mpfr_t tmp)
	{
		mpfr_set(tmp, coeffs[degree], MPFR_RNDN);
		mpfr_set(rx, x, MPFR_RNDN);
		
		for(int i=degree-1;i>=0;i--){
			mpfr_mul(tmp, tmp, rx, MPFR_RNDN);
			
			mpfr_add(tmp, tmp, coeffs[i], MPFR_RNDN);
			mpfr_fix(tmp, calcFracBits[i]);
			
			calcMsbs[i]=std::max(calcMsbs[i], (int)mpfr_get_exp(tmp));
		}
	}
	
	/* Calculate the maximum MSB+1 of the result of each multiply-accumulate stage
		Returns a vector [0..degree) which is the position of the highest set one that
		can occur. e.g. if the maximum achieved is [0.5,1), then the msb will be 0
	*/
	void calc_msbs()
	{	
		calcMsbs=calcFracBits;	// By default they have no MSBs
		for(int i=0;i<calcMsbs.size();i++){
			calcMsbs[i] *= -1;
		}
		
		MPFRVec extrema=find_polynomial_and_subpolynomial_extremal_points(coeffs, range); // get possibly redundant set of points
		assert(extrema.size()>0);
		
		mpfr_t x, rx, tmp;
		mpfr_inits2(getToolPrecision(), x, rx, tmp, (mpfr_ptr)0);
		
		fprintf(stderr, "Extrema:\n");
		for(int i=0;i<extrema.size();i++){
			mpfr_fprintf(stderr, "[%Rf,  %Rf,  %Rf]\n", range[0], extrema[i], range[1]);
			
			mpfr_sub_d(x, extrema[i], 10*pow(2.0, -inputFracBits), MPFR_RNDD);
			mpfr_fix(x, inputFracBits, MPFR_RNDD);
			
			for(int j=0;j<20;j++){
				if(mpfr_greaterequal_p(range[0], x) && mpfr_lessequal_p(x, range[1])){
					update_intermediates(x, rx, tmp);
				}
				mpfr_nextabove(x);
			}			
		}
		
		mpfr_clears(x, rx, tmp, (mpfr_ptr)0);
	}
public:	
	ConcretePoly()
		: degree(INT_MAX)
	{}
	
	ConcretePoly(MPFRVec _coeffs, MPFRVec _range, int _inputFracBits, std::vector<int> _coeffFracBits, int _calcFracBits, int _resFracBits)
		: degree(_coeffs.size()-1)
		, coeffs(_coeffs)
		, range(_range)
		, inputFracBits(_inputFracBits)
		, coeffFracBits(_coeffFracBits)
		, calcFracBits(degree, _calcFracBits)
	{
		fprintf(stderr, "ConcretePoly(degree=%d, calcFracBits=%d, outFracBits=%d)\n", _coeffs.size()-1, _calcFracBits, _resFracBits);
		
		inputMsb=std::max(mpfr_get_exp(_range[0]) , mpfr_get_exp(_range[1]));
		
		assert(coeffFracBits.size()==degree+1);
		
		calcFracBits[0]=_resFracBits;
		
		mpfr_t tmp;
		mpfr_init2(tmp, getToolPrecision());
		
		for(unsigned i=0;i<=degree;i++){			
			mpfr_set(coeffs[i], _coeffs[i], MPFR_RNDN);
			mpfr_fix(coeffs[i], coeffFracBits[i]);
			
			if(!mpfr_equal_p(coeffs[i], _coeffs[i]))
				throw std::invalid_argument("ConcretePoly - coefficient does not fit in specified format.");
			
			coeffMsbs.push_back(mpfr_get_exp(coeffs[i]));
		}
		
		calc_msbs();
		
		mpfr_clear(tmp);
	}
	
	void eval(mpfr_t res, mpfr_t x)
	{
		mpfr_t tmp, rx;
		mpfr_init2(tmp, getToolPrecision());	// assumed to be big enough not to truncate anything...
		mpfr_init2(rx, getToolPrecision());
		
		mpfr_set(tmp, coeffs[degree], MPFR_RNDN);	// should automatically be right bit-width
		mpfr_set(rx, x, MPFR_RNDN);
		
		for(int i=degree-1;i>=0;i--){
			
			mpfr_mul(tmp, tmp, rx, MPFR_RNDN);
			
			mpfr_add(tmp, tmp, coeffs[i], MPFR_RNDN);
			mpfr_fix(tmp, calcFracBits[i]);
		}
		
		mpfr_set(res, tmp, MPFR_RNDN);
		
		mpfr_clears(tmp, rx, (mpfr_ptr)0);
	}
	
	void dump(FILE *dst, const char *prefix="") const
	{
		fprintf(dst, "%sPoly[degree=%u,inputWidth=%u,inputFracBits=%u]\n", prefix, degree, inputMsb+inputFracBits, inputFracBits);
		mpfr_fprintf(dst, "%s  bits=%u, coeff=%Rf\n", prefix, coeffMsbs[degree]+coeffFracBits[degree]+1, coeffs[degree]);
		for(int i=degree-1;i>=0;i--){
			mpfr_fprintf(dst, "%s  stage %d : mul_in_bits=[%ux%u], stage_out_bits=%u, coeff=%Rf\n", prefix, i,
				((inputMsb+inputFracBits)+1),((coeffMsbs[i]+coeffFracBits[i])+1), 
				calcMsbs[i]+calcFracBits[i]+1, coeffs[i]
			);			
		}
	}
};


/*
struct format_t
{
	bool isSigned;
	int msb;
	int lsb;
	
	int width() const
	{ return (msb-lsb)+(isSigned?1:0); }
};

format_t to_signed(format_t x)
{
	if(x.isSigned)
		return x;
	else
		format_t(true, x.msb+1, lsb);
}

format_t mul(format_t a, format_t b)
{
	if(a.isSigned||b.isSigned){
		a=to_signed(a);
		b=to_signed(b);
	}
	return format_t(a.isSigned, 
}

// Input is always unsigned, coefficients are always signed
class ConcretePolyOperator
{
	std::string sign_extend(int dstW, int dstL, std::string src, int srcW, int srcL)
	{
		assert(dstL <= srcL);
		assert(dstL+dstW >= srcL+srcW);
		
		std::stringstream acc;
		acc<<"(";
		if(dstL+dstW > srcL+srcW){
			for(int i=0;i<(dstW+dstL)-(srcW+srcL);i++){
				acc<<name<<"("<<srcW-1<<")&";
			}
		}
		acc<<name;
		if(dstL < srcL){
			acc << " & \"";
			for(int i=0;i<(srcL-dstL);i++){
				acc << "0";
			}
			acc<<"\"";
		}		
		acc<<")";
		return acc.str();
	}
	
	std::string round_signed(int dstL, std::string src, int srcW, int srcL)
	{
		assert(dstL >= srcL);
		
		std::stringstream acc;
		acc<<"(";
		if(srcL==dstL){
			acc<<name;
		}else{
			acc<<"\"";
			for(int i=srcW-1;i>=srcL;i--){
				acc<<"0";
			}
			acc<<"0";
			#error "Here"
			
			acc<<"
		}
		acc<<")";
		return acc.str();
	}
	
	fixed_signal_t render_mac_stage(
		unsigned i,
		std::string prefix,
		std::string prevName, int prevWidth, int prevLsb,	// prev is always signed
		std::string inputName, int inputWidth, int inputLsb, // input is always unsigned
		std::string coeffName, int coeffWidth, int coeffLsb, // coeff is always signed
		std::string outputName, int outputWidth, int outputLsb	// output is always signed
	){
		//stage_mul_i <= prev * input;
		//stage_macc_i <= stage_mul_i + coeff;
		//stage_res_i <= round(stage_macc_i);
		
		int mul_width=prevWidth+1+inputWidth;
		int mul_lsb=prevLsb+inputLsb;
		vhdl << declare(prefix+"_mul", mul_width) << "<=" << prevName << " * (\'0'&"<<inputName<<");\n";
		
		int add_lsb=std::min(mul_lsb,coeffLsb);
		int add_width=std::max(mul_width+mul_lsb,coeffWidth+coeffLsb)-add_lsb;
		
		vhdl << declare(prefix+"_mul_ext", add_width) << "<=" << sign_extend(add_width, add_lsb, prefix+"_mul", mul_width, mul_lsb) << ";\n";
		vhdl << declare(prefix+"_coeff_ext", add_width) << "<=" << sign_extend(add_width, add_lsb, coeffName, coeffWidth, coeffLsb) << ";\n";
		vhdl << declare(prefix+"_add"), add_width) << "<=" << (prefix+"_mul_ext") << "+" << (prefix+"_coeff_ext") << ";\n";
		
		vhdl << declare(outputName, outputWidth) << "<=" << round_signed(outputWidth, outputLsb, (prefix+"_add"), add_width, add_lsb) <<";\n";
	}
	
	ConcretePolyOperator(
		Target* target,
		unsigned nPoly,
		ConcretePoly *aPoly
	)
		: Operator(target)
	{
		setName("ConcretePoly", "");
		setCopyrightString("Imperial College 2012");
	}
};
*/

}; // float_approx
}; // random
}; // flopoco

#endif
