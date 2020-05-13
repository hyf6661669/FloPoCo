#ifndef __FPLOG_HPP
#define __FPLOG_HPP
#include <vector>
#include <sstream>

#include "Operator.hpp"

namespace flopoco{

#define MAXRRSTAGES 2000 // 4000 bits of accuracy should be enough for anybody

	class FPLog : public Operator
	{
	public:
		//		Overloading the virtual functions of Operator
		static void emulate(TestCase * tc, int wE, int wF);
		static void buildStandardTestCases(OperatorPtr op, TestCaseList* tcl, int wE, int wF);
		static TestList unitTest(int index);
		/**Overloading the function of Operator with a function that tests only positive FP numbers (full range)*/
		static TestCase* buildRandomTestCase(OperatorPtr op, int i, int wE, int wF);
		// User-interface stuff
		/** Factory method */
		static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target , vector<string> &args);
		static void registerFactory();
	protected:
		int wE, wF;
	};
	
}
#endif
