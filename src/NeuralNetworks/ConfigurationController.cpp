// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "ConfigurationController.hpp"

using namespace std;
namespace flopoco
{




    ConfigurationController::ConfigurationController(Target* target, int input_bit_width_, int Const_bit_width_, int No_of_Products_) : Operator(target), input_bit_width(input_bit_width_), Const_bit_width(Const_bit_width_),No_of_Products(No_of_Products_)
    {
        srcFileName="ConfigurationController";
		// definition of the name of the operator
		ostringstream name;
        name << "ConfigurationController_" << Const_bit_width_ << "_" << No_of_Products_;
        setName(name.str());
        // Copyright
        setCopyrightString("UNIVERSITY of Kassel 2017");


    }

    void ConfigurationController::emulate(TestCase * tc)
    {
	}

    void ConfigurationController::buildStandardTestCases(TestCaseList * tcl)
    {
		// please fill me with regression tests or corner case tests!
	}

    OperatorPtr ConfigurationController::parseArguments(Target *target, vector<string> &args)
    {
        int param0, param1, param2;
        UserInterface::parseInt(args, "input_bit_width", &param0); // param0 has a default value, this method will recover it if it doesnt't find it in args,
        UserInterface::parseInt(args, "const_bit_width", &param1);
        UserInterface::parseInt(args, "no_of_products", &param2);
        if(param1 == -1)
        {
            param1 = param0;
        }

        return new ConfigurationController(target, param0, param1, param2);
	}
	
    void ConfigurationController::registerFactory()
    {
        UserInterface::add("ConfigurationController", // name
                                             "My first ConfigurationController.", // description, string
											 "NeuralNetworks", // category, from the list defined in UserInterface.cpp
											 "", //seeAlso
											 // Now comes the parameter description string.
											 // Respect its syntax because it will be used to generate the parser and the docs
											 // Syntax is: a semicolon-separated list of parameterDescription;
											 // where parameterDescription is parameterName (parameterType)[=defaultValue]: parameterDescriptionString 
                                             "input_bit_width(int)=16: input word size; \
                                             const_bit_width(int)=-1: coefficient word size, per default the same as the input bit width;\
                                             no_of_products(int)=1: the nomber of products to accumulate",
											 // More documentation for the HTML pages. If you want to link to your blog, it is here.
											 "Feel free to experiment with its code, it will not break anything in FloPoCo. <br> Also see the developper manual in the doc/ directory of FloPoCo.",
                                             ConfigurationController::parseArguments
											 ) ;
	}

}//namespace
