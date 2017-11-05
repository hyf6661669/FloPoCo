/* Each Operator declared within the flopoco framework has 
   to inherit the class Operator and overload some functions listed below*/
#include "Operator.hpp"

/* This file contains a lot of useful functions to manipulate vhdl */
#include "utils.hpp"

/*  All flopoco operators and utility functions are declared within
    the flopoco namespace.
    You have to use flopoco:: or using namespace flopoco in order to access these
    functions.
*/


namespace flopoco {

	// new operator class declaration
    class SOP_KCM : public Operator {

        int input_bit_width;
        int Const_bit_width;
        int No_of_Products;
        int output_bit_width;
        int g;// nomber of gard bits

        bool signed_calculation; // if false there are just positive constnats and inputs
        bool faithful_rounding; // if false there are just positive constnats and inputs

	public:

        SOP_KCM(Target* target,int inputWordSize,int constantWordSize, int parameterNo, int _defaultConstant);
        ~SOP_KCM() {};

        //void emulate(TestCase * tc);
        //void buildStandardTestCases(TestCaseList* tcl);

		/** Factory method that parses arguments and calls the constructor */
		static OperatorPtr parseArguments(Target *target , vector<string> &args);
		/** Factory register method */ 
		static void registerFactory();

    private:
        int LUT_bit_width;
        string generateInitStringFor(int weight, unsigned int LUTNo);


	};


}//namespace
