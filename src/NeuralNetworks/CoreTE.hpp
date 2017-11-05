#include "Operator.hpp"

/* This file contains a lot of useful functions to manipulate vhdl */
#include "utils.hpp"


namespace flopoco {

    // new operator class declaration
    class CoreTE : public Operator
    {
    public:
        int input_bit_width;
        int const_bit_width;
        int output_bit_width;
        int no_of_products;


    public:
        // definition of some function for the operator

        // constructor, defined there with two parameters (default value 0 for each)
        CoreTE(Target* target, int input_bit_width_, int Const_bit_width_, int No_of_Products_);

        // destructor
        ~CoreTE() {};

        /** Factory method that parses arguments and calls the constructor */
        static OperatorPtr parseArguments(Target *target , vector<string> &args);

        /** Factory register method */
        static void registerFactory();
    };

};
