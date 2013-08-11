#ifndef FPNormalCDF_HPP
#define FPNormalCDF_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"
#include "FPNumber.hpp"

namespace flopoco{
	class FPNormalCDF : public Operator
	{
	public:
		/**
		 * The FPSqrt constructor
		 * @param[in]		target		the target device
		 * @param[in]		wInE		the width of the exponent for the input f-p number X
		 * @param[in]		wInF		the width of the fraction for the input f-p number X
		* @param[in]			wOutE	the width of the exponent for the output f-p number R
		 * @param[in]		wOutF	the width of the fraction for the output f-p number R
		* @param[in]			minInputValue	for values in [-inf,minInputValue], the implementation is allowed to return zero. By default is -8. Must be less than zero.
		 */
		FPNormalCDF(Target* target, int wInE, int wInF, int wOutE, int wOutF, double underflowPoint=0);

		~FPNormalCDF();

		void emulate(TestCase * tc);
	
		void buildStandardTestCases(TestCaseList* tcl);
		TestCase* buildRandomTestCase(int i);

	private:
		int wInE, wInF, wOutE, wOutF;
		double minInputValue;
	
		mpfr_t m_minInputX;
	
		int wSqrXF;	// Width of signal sqrX
		int wBxE, wBxF;	// Width of the signal Bx
		int fixXLSB, fixXMSB;	// Fixed-point type going into function eval
		int FxLSB, FxMSB;		// Type coming out of function eval
		int wFxE, wFxF;	// Floating point type of function eval
	
		bool m_debugOutputs;
	
		FunctionEvaluator *opFx;
	
		Function *funcFx;	// Function requested from FunctionEvaluator
		std::vector<mpz_class> emulate_Fx(mpz_class x);	// Fixed-point evaluation of Fx
		void emulate_Fx(mpfr_t r, mpfr_t x);	// Converts to and from fixed-point, and returns one of the legal outputs
	
		void NCD(mpfr_t r, mpfr_t x, mpfr_rnd_t mode=MPFR_RNDN);
		void NCDInv(mpfr_t r, mpfr_t x);
		void F(mpfr_t r, mpfr_t x, mpfr_rnd_t mode=MPFR_RNDN);
		void B(mpfr_t r, mpfr_t x, mpfr_rnd_t mode=MPFR_RNDN);
	
		void CalcMinX();
		int CalcFxOutputLSB(int relativeBitsNeeded);
		int CalcFxInputMSB();
	};
}
#endif //FPNormalCDF_HPP
