#ifndef IntConstMultOptTernary_HPP
#define IntConstMultOptTernary_HPP

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
    Integer constant multiplication using minimum number of adders using ternary adders

*/


namespace flopoco{

    class IntConstMultOptTernary : public ConstMultPAG
	{
	public:
		/** The standard constructor, inputs the number to implement */ 
        IntConstMultOptTernary(Target* target, int wIn, int coeff, bool syncInOut=true);

//		void emulate(TestCase* tc);
//		void buildStandardTestCases(TestCaseList* tcl);

        static OperatorPtr parseArguments( Target *target, vector<string> &args );
        static void registerFactory();

	private:
		int coeff;  /**< The constant */

	};
}

#endif //HAVE_PAGLIB
#endif
