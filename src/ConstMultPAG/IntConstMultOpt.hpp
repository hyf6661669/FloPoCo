#ifndef INTCONSTMULTOPT_HPP
#define INTCONSTMULTOPT_HPP

#ifdef HAVE_PAGLIB

#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include <cstdlib>

#include "../Operator.hpp"
#include "IntAddSubCmp/IntAdder.hpp"
#include "../ConstMultPAG/ConstMultPAG.hpp"

/**
	Integer constant multiplication using minimum number of adders due to

	Gustafsson, O., Dempster, A., Johansson, K., Macleod, M., & Wanhammar, L. (2006).
	Simplified Design of Constant Coefficient Multipliers. Circuits, Systems, and Signal Processing, 25(2), 225–251.


	All constants up to 19 bit will be realized optimal using precomputed tables provided by the SPIRAL project (http://spiral.ece.cmu.edu/mcm/).

*/


namespace flopoco{

    class IntConstMultOpt : public ConstMultPAG
	{
	public:
		/** The standard constructor, inputs the number to implement */ 
        IntConstMultOpt(Target* target, int wIn, int c, bool syncInOut=true);

//		void emulate(TestCase* tc);
//        void buildStandardTestCases(TestCaseList* tcl);

        static OperatorPtr parseArguments( Target *target, vector<string> &args );
        static void registerFactory();

	private:
        int coeff;  /**< The constant */
        int wIn;

        void generateAOp(int a, int b, int c, int eA, int eB, int signA, int signB, int preFactor=1);
		void buildAdderGraph(int c, int preFactor=1);
		void generateVHDLOptSCM(int c);

        stringstream adderGraph;
	};
}

#endif //HAVE_PAGLIB
#endif
