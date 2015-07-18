#ifndef FixFunctionByPiecewisePoly_HPP
#define FixFunctionByPiecewisePoly_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "../Operator.hpp"
#include "../Table.hpp"
#include "FixFunction.hpp"
#include "PiecewisePolyApprox.hpp"

namespace flopoco{


	/** The FixFunctionByPiecewisePoly class */
	class FixFunctionByPiecewisePoly : public Operator
	{
	public:

		/** a subclass to generate a table of coefficients */
		class CoeffTable: public Table {
		public:
			CoeffTable(Target* target, int wIn, int wOut, PiecewisePolyApprox* polyApprox, bool addFinalRoundBit, int finalRoundBitPos);
			mpz_class function(int x);
			PiecewisePolyApprox* polyApprox; // don't understand why C++ won't let me use that of FixFunctionByPiecewisePoly
			bool addFinalRoundBit;
			int finalRoundBitPos;
		};


		/**
		 * The FixFunctionByPiecewisePoly constructor
			 @param[string] func    a string representing the function, input range should be [0,1]
			 @param[int]    lsbIn   input LSB weight
			 @param[int]    msbOut  output MSB weight, used to determine wOut
			 @param[int]    lsbOut  output LSB weight
			 @param[int]    degree  degree of polynomial approximation used. Controls the trade-off between tables and multipliers.
			 @param[bool]   finalRounding: if false, the operator outputs its guard bits as well, saving the half-ulp rounding error. 
			                 This makes sense in situations that further process the result with further guard bits.
			 @param[bool]   plainStupidVHDL: if true, generate * and +; if false, use BitHeap-based FixMultAdd
			 
			 One could argue that MSB weight is redundant, as it can be deduced from an analysis of the function. 
			 This would require quite a lot of work for non-trivial functions (isolating roots of the derivative etc).
			 So this is currently left to the user.
		 */
		FixFunctionByPiecewisePoly(Target* target, std::string func, int lsbIn, int msbOut, int lsbOut, int degree, bool finalRounding = true,  std::map<std::string, double> inputDelays = emptyDelayMap);

		/**
		 * FixFunctionByPiecewisePoly destructor
		 */
		~FixFunctionByPiecewisePoly();
		
		void emulate(TestCase * tc);

		void buildStandardTestCases(TestCaseList* tcl);

	private:
		int degree;
		PiecewisePolyApprox *polyApprox;
		FixFunction *f; 
		bool finalRounding;
	};

}

#endif
