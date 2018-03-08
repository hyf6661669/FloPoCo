//
// Created by Viktor Schmidt.
//

#include "MonotoneFunctionLUT.hpp"
// general c++ library for manipulating streams
#include <iostream>
#include <sstream>
#include <bitset>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

using namespace std;
namespace flopoco {
    MonotoneFunctionLUT::MonotoneFunctionLUT(Target* target, int inputWidth_, int outputWidth_) : inputWidth(inputWidth_), outputWidth(outputWidth_), Table(target, inputWidth_, outputWidth_) {
        /* constructor of the U`serDefinedOperator
           Target is the targeted FPGA : Stratix, Virtex ... (see Target.hpp for more informations)
           param0 and param1 are some parameters declared by this Operator developpers,
           any number can be declared, you will have to modify
           -> this function,
           -> the prototype of this function (available in UserDefinedOperator.hpp)
           ->  main.cpp to uncomment the RegisterFactory line
        */
        /* In this constructor we are going to generate an operator that
             - takes as input three bit vectors X,Y,Z of lenght param0,
             - treats them as unsigned integers,
             - sums them
             - and finally outputs the concatenation of the the most significant bit of the sum and its last param1 bits.
             Don't ask why

             All the vhdl code needed by this operator has to be generated in this constructor */

        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="MonotoneFunctionLUT";

        // definition of the name of the operator
        ostringstream name;
        name << "MonotoneLUT" << inputWidth << "_" << outputWidth;
        setName(name.str());
        // Copyright
        setCopyrightString("Viktor Schmidt 2018");

        /* SET UP THE IO SIGNALS
           Each IO signal is declared by addInput(name,n) or addOutput(name,n)
           where name is a string that stands for the name of the variable and
           n is an integer (int)   that stands for the length of the corresponding
           input/output */

        // declaring inputs
//        addInput("i" , inputWidth);
//        addFullComment(" addFullComment for a large comment ");
//        addComment("addComment for small left-aligned comment");

        // declaring output
//        addOutput("o" , outputWidth);

        /* Some piece of information can be delivered to the flopoco user if  the -verbose option is set
           [eg: flopoco -verbose=0 UserDefinedOperator 10 5 ]
           , by using the REPORT function.
           There is three level of details
           -> INFO for basic information ( -verbose=1 )
           -> DETAILED for complete information, includes the INFO level ( -verbose=2 )
           -> DEBUG for complete and debug information, to be used for getting information
           during the debug step of operator development, includes the INFO and DETAILED levels ( -verbose=3 )
        */
        // basic message
        REPORT(INFO,"Declaration of MonotoneLUT \n");

        // more detailed message
        REPORT(DETAILED, "this operator has received two parameters " << inputWidth << " and " << outputWidth);

        // debug message for developper
        REPORT(DEBUG,"debug of MonotoneLUT");


        /* vhdl is the stream which receives all the vhdl code, some special functions are
           available to smooth variable declaration and use ...
           -> when first using a variable (Eg: T), declare("T",64) will create a vhdl
           definition of the variable : signal T and includes it it the header of the architecture definition of the operator

           Each code transmited to vhdl will be parsed and the variables previously declared in a previous cycle will be delayed automatically by a pipelined register.
        */
        sollya_lib_parse_string("[0;1]");
        fS= sollya_lib_parse_string("sin(3.14/4*x);");

        //declare("lut_in", inputWidth);
//        vhdl << tab << "case i is" << endl;

//        mpz_class lut_in(0), lut_out;


//        for(int i = 0; i < pow(2, inputWidth); ++i) {
//            REPORT(DEBUG,"eval#" << i);
//            lut_in = mpz_class(i);
//            mpfr_set_d(mpX, i, GMP_RNDN);
//            eval(lut_out, lut_in);
//            mpfr_get_str(in_str, &in_exp, 2, inputWidth, mpX, GMP_RNDN);
//            mpfr_get_str(out_str, &out_exp, 2, outputWidth, mpR, GMP_RNDN);
//            vhdl << tab << tab << "when " << "\"" << lut_in.get_str(2) << "\" => o <= \"" << lut_out.get_str(2) << "\";" << endl;

//        }

//        vhdl << tab << "end case;" << endl;
    };


    void MonotoneFunctionLUT::emulate(TestCase * tc) {
        /* This function will be used when the TestBench command is used in the command line
           we have to provide a complete and correct emulation of the operator, in order to compare correct output generated by this function with the test input generated by the vhdl code */
        /* first we are going to format the entries */
        mpz_class sx = tc->getInputValue("X");

        /* then we are going to manipulate our bit vectors in order to get the correct output*/
        mpz_class sr;
        eval(sr, sx);

        /* at the end, we indicate to the TestCase object what is the expected
           output corresponding to the inputs */
        tc->addExpectedOutput("Y",sr);
    }


    void MonotoneFunctionLUT::buildStandardTestCases(TestCaseList * tcl) {
        // please fill me with regression tests or corner case tests!
    }


    void MonotoneFunctionLUT::eval(mpz_class& r, mpz_class x) const
    {
        mpfr_t mpX, mpR;
        mpfr_init2(mpX, inputWidth+2);
        mpfr_init2(mpR, outputWidth*3);
        sollya_lib_set_prec(sollya_lib_constant_from_int(inputWidth*2));

        mpfr_set_z(mpX, x.get_mpz_t(), GMP_RNDN);
        mpfr_div_2si(mpX, mpX, inputWidth, GMP_RNDN);

        sollya_lib_evaluate_function_at_point(mpR, fS, mpX, NULL);

        mpfr_mul_2si(mpR, mpR, outputWidth, GMP_RNDN);
        mpfr_get_z(r.get_mpz_t(), mpR, GMP_RNDN);

        REPORT(FULL,"function() input is:"<<mpfr_get_d(mpX, GMP_RNDN));
        REPORT(FULL,"function() output is:"<<mpfr_get_d(mpR, GMP_RNDN));

        mpfr_clear(mpX);
        mpfr_clear(mpR);
    }


    OperatorPtr MonotoneFunctionLUT::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args) {
        int param0, param1;
        UserInterface::parseInt(args, "inputWidth", &param0); // param0 has a default value, this method will recover it if it doesnt't find it in args,
        UserInterface::parseInt(args, "outputWidth", &param1);
        return new MonotoneFunctionLUT(target, param0, param1);
    }

    void MonotoneFunctionLUT::registerFactory(){
        UserInterface::add("MonotoneFunctionLUT", // name
                           "Generates a LUT.", // description, string
                           "Miscellaneous", // category, from the list defined in UserInterface.cpp
                           "", //seeAlso
                // Now comes the parameter description string.
                // Respect its syntax because it will be used to generate the parser and the docs
                // Syntax is: a semicolon-separated list of parameterDescription;
                // where parameterDescription is parameterName (parameterType)[=defaultValue]: parameterDescriptionString
                           "inputWidth(int)=16: Input bit count; \
                        outputWidth(int)=8: Output bit count",
                // More documentation for the HTML pages. If you want to link to your blog, it is here.
                           "Feel free to experiment with its code, it will not break anything in FloPoCo. <br> Also see the developer manual in the doc/ directory of FloPoCo.",
                           MonotoneFunctionLUT::parseArguments
        ) ;
    }

    mpz_class MonotoneFunctionLUT::function(int x) {
        mpz_class lut_in(x), lut_out;

        eval(lut_out, lut_in);

        return lut_out;
    }

}//namespace