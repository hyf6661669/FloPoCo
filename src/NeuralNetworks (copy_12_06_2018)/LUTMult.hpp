/* Each Operator declared within the flopoco framework has 
   to inherit the class Operator and overload some functions listed below*/
#include "Operator.hpp"

/* This file contains a lot of useful functions to manipulate vhdl */
#include "utils.hpp"
#include <vector>

/*  All flopoco operators and utility functions are declared within
    the flopoco namespace.
    You have to use flopoco:: or using namespace flopoco in order to access these
    functions.
*/

namespace flopoco {

    class LUTMultCell : public Operator
    {
    public:
        LUTMultCell(Target* target,unsigned int inputBitWidth_, unsigned int LutBitWidth_ = 4, bool implementFullAddder_ = true);
        ~LUTMultCell() {};
        unsigned int inputBitWidth;
        unsigned int LutBitWidth;
        bool implementFullAddder;
    };


	// new operator class declaration
//    class LUTMult : public Operator
//    {
//	public:
//		/* operatorInfo is a user defined parameter (not a part of Operator class) for
//           stocking information about the operator. The user is able to defined any number of parameter in this class, as soon as it does not affect Operator parameters undeliberatly*/
//        unsigned int LutBitWidth;
//        unsigned int inputWordsize;
//	public:
//		// definition of some function for the operator

//		// constructor, defined there with two parameters (default value 0 for each)
//        LUTMult(Target* target,unsigned int inputWordsize_, unsigned int LutBitWidth_ = 5);

//		// destructor
//        ~LUTMult() {};


//		// Below all the functions needed to test the operator
//		/* the emulate function is used to simulate in software the operator
//		   in order to compare this result with those outputed by the vhdl opertator */
//		void emulate(TestCase * tc);

//		/* function used to create Standard testCase defined by the developper */
//		void buildStandardTestCases(TestCaseList* tcl);


//		/* function used to bias the (uniform by default) random test generator
//		   See FPExp.cpp for an example */
//		// TestCase* buildRandomTestCase(int i);

//		/** Factory method that parses arguments and calls the constructor */
//		static OperatorPtr parseArguments(Target *target , vector<string> &args);
		
//		/** Factory register method */
//		static void registerFactory();


//	};


}//namespace
