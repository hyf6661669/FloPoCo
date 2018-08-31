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
    class Multiplier_B_B_E : public Operator
    {

    public:

        Multiplier_B_B_E(Target* target, int inputWordSize);
        ~Multiplier_B_B_E(){}

        //void emulate(TestCase * tc);
        //void buildStandardTestCases(TestCaseList* tcl);

        /** Factory method that parses arguments and calls the constructor */
        static OperatorPtr parseArguments(Target *target , vector<string> &args);
        /** Factory register method */
        static void registerFactory();

    private:
        int input_bit_width;

    };


}//namespace
