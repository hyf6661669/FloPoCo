// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "ConfigurationGenerator.hpp"

using namespace std;
namespace flopoco {




	ConfigurationGenerator::ConfigurationGenerator(Target* target, int coeffWordSize_) : Operator(target), coeffWordSize(coeffWordSize_) {
		/* constructor of the ConfigurationGenerator
		   Target is the targeted FPGA : Stratix, Virtex ... (see Target.hpp for more informations)
		   coeffWordSize and param1 are some parameters declared by this Operator developpers, 
		   any number can be declared, you will have to modify 
		   -> this function,  
		   -> the prototype of this function (available in ConfigurationGenerator.hpp)
		   ->  main.cpp to uncomment the RegisterFactory line
		*/
		/* In this constructor we are going to generate an operator that 
			 - takes as input three bit vectors X,Y,Z of lenght coeffWordSize, 
			 - treats them as unsigned integers, 
			 - sums them 
			 - and finally outputs the concatenation of the the most significant bit of the sum and its last param1 bits.
			 Don't ask why

			 All the vhdl code needed by this operator has to be generated in this constructor */

		// definition of the source file name, used for info and error reporting using REPORT 
		srcFileName="ConfigurationGenerator";
		this->useNumericStd();
		// definition of the name of the operator
		ostringstream name;
		name << "ConfigurationGenerator_" << coeffWordSize;
		setName(name.str());
		// Copyright 
		setCopyrightString("UNIVERSITY of Kassel 2017");

		/* SET UP THE IO SIGNALS
		   Each IO signal is declared by addInput(name,n) or addOutput(name,n) 
		   where name is a string that stands for the name of the variable and 
		   n is an integer (int)   that stands for the length of the corresponding 
		   input/output */

		// declaring inputs
		addInput ("ce_i");
		addInput ("coeff_i" , coeffWordSize);
		addInput ("ctrl_int" , 3);

		// declaring output
		addOutput("cdo_lsb_o" , coeffWordSize+5);
		addOutput("cdo_msb_o" , coeffWordSize+5);

        if((targetID == "Virtex5") || (targetID == "Virtex6") || (targetID == "Virtex7") || (targetID == "Spartan6"))
        {
					//use coeff2lut_prim

        }
        else
        {
					//use coeff2lut_no_prim
        }
    }

	
	void ConfigurationGenerator::emulate(TestCase * tc) {
	}


	void ConfigurationGenerator::buildStandardTestCases(TestCaseList * tcl) {
	}





	OperatorPtr ConfigurationGenerator::parseArguments(Target *target, vector<string> &args) {
		int coeffWordSize;
		UserInterface::parseInt(args, "coeffWordSize", &coeffWordSize); // coeffWordSize has a default value, this method will recover it if it doesnt't find it in args, 
		return new ConfigurationGenerator(target, coeffWordSize);
	}
	
	void ConfigurationGenerator::registerFactory(){
		UserInterface::add("ConfigurationGenerator", // name
											 "An heavily commented example operator to start with FloPoCo.", // description, string
											 "NeuralNetworks", // category, from the list defined in UserInterface.cpp
											 "", //seeAlso
											 // Now comes the parameter description string.
											 // Respect its syntax because it will be used to generate the parser and the docs
											 // Syntax is: a semicolon-separated list of parameterDescription;
											 // where parameterDescription is parameterName (parameterType)[=defaultValue]: parameterDescriptionString 
											 "coeffWordSize(int): Coefficient word size;",
											 // More documentation for the HTML pages. If you want to link to your blog, it is here.
											 "",
											 ConfigurationGenerator::parseArguments
											 ) ;
	}

}//namespace
