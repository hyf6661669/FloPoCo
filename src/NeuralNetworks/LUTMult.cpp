// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "LUTMult.hpp"
#include "../PrimitiveComponents/Xilinx/Xilinx_CFGLUT5.hpp"

using namespace std;
namespace flopoco {

    LUTMultCell::LUTMultCell (Target* target, unsigned int inputBitWidth_, unsigned int LutBitWidth_, bool implementFullAddder_) : Operator(target)
    {
        inputBitWidth = inputBitWidth_;
        LutBitWidth = LutBitWidth_;
        implementFullAddder = implementFullAddder_;

        srcFileName="LUTMult";

        ostringstream name;
        name << "LUTMult" << inputBitWidth << "_" << LutBitWidth  << "_" << implementFullAddder;
        setName(name.str());
        setCopyrightString("ACME and Co 2010");

        addInput("ce"); //clk enable
        addInput("cdi");//configuration data in
        addInput ("x" ,inputBitWidth);
        addOutput("result");

        if(implementFullAddder)
        {
            addInput ("FA");
            addInput ("Cin");
            addOutput("Cout");
        }

        //addFullComment(" addFullComment for a large comment ");
        //addComment("addComment for small left-aligned comment");

        Xilinx_CFGLUT5 *myCFGLUT = new Xilinx_CFGLUT5(target);
            addToGlobalOpList(myCFGLUT);
            string instanceName = "CGFLUT5_inst";


            inPortMap(myCFGLUT, "ce","ce");
            inPortMap(myCFGLUT, "i0","x(0)");
            inPortMap(myCFGLUT, "i1","x(1)");
            inPortMap(myCFGLUT, "i2","x(2)");
            inPortMap(myCFGLUT, "i3","x(3)");
            inPortMap(myCFGLUT, "i4","x(4)");

            outPortMap(myCFGLUT, "o5","x()");
            outPortMap(myCFGLUT, "o6","x()");
            outPortMap(myCFGLUT, "cout","x()");

            addInput("ce"); //clk enable
            addInput("cdi");//configuration data in

//            //selector inputs
//            addInput("i0");
//            addInput("i1");
//            addInput("i2");
//            addInput("i3");
//            addInput("i4");

//            // declaring output
//            addOutput( "o5"); //4 LUT output
//            addOutput( "o6"); //5 LUT output
//            addOutput( "cdo"); //configuration data out

            //string instanceName = "LUTMult_inst";
            vhdl << instance(myCFGLUT,instanceName);
    }



//    LUTMult::LUTMult(Target* target, int inputWordsize_, int LutBitWidth_=5) : Operator(target), inputWordsize(inputWordsize_), LutBitWidth(LutBitWidth_)
//{
//        srcFileName="LUTMult";

//		// definition of the name of the operator
//		ostringstream name;
//        name << "LUTMult" << inputWordsize_ << "_" << LutBitWidth_;
//		setName(name.str());
//		// Copyright
//		setCopyrightString("ACME and Co 2010");

//        vhdl << tab << declare("ce",1,false) << "<= '0';" << endl;
//        vhdl << tab << declare("cdi",1,false) << " <= T(0);" << endl;
//        vhdl << tab << declare("x",5) << ";" << std::endl;  // << "<= (others=> '0');" << endl;
//        //vhdl << tab << declare("LUToutput") << endl;


//        int LUT_per_stage = 8;
//        int No_of_stages = 3;
//        int LUT_bit_width = 5;
//        int input_bit_width = 5;


//        //for(int stage = 0; stage < No_of_stages; ++stage)
//        //{
//        //    vhdl << tab << declare(join("LUT_stage_",to_string(stage)),LUT_per_stage) << ";" << std::endl;
//        //}

//        for(int stage = 0; stage < No_of_stages; ++stage)
//        {
//            for(int LUT_No = 0; LUT_No< LUT_per_stage; ++LUT_No)
//            {
//                Xilinx_CFGLUT5 *myCFGLUT = new Xilinx_CFGLUT5(target);
//                addToGlobalOpList(myCFGLUT);
//                //inPortMap(myCFGLUT, "clk","clk");
//                inPortMap(myCFGLUT, "ce","ce");
//                inPortMap(myCFGLUT, "cdi","cdi");

//                for(int i=0; i< LUT_bit_width; ++i)
//                {
//                     inPortMap(myCFGLUT, join("i", to_string(i)),join("x(",to_string(i+stage*LUT_bit_width),")"));
//                }

//                string outputSignalName = "LUT_stage_" + to_string(stage) + "(" +to_string(LUT_No) +")";
//                outPortMap(myCFGLUT, "o6",outputSignalName,false);

//                string instanceName = "CFGLUT_inst_" + to_string(stage) + "_LUT_No_" +to_string(LUT_No);
//                vhdl << instance(myCFGLUT,instanceName);

//                //Full adder:
//                vhdl << tab << declare(join("instanceName",Fullad),1,false) << "<= '0';" << endl;
//            }
//        }
        
//	};

//    void LUTMult::emulate(TestCase * tc) {
//	}

//    void LUTMult::buildStandardTestCases(TestCaseList * tcl) {
//		// please fill me with regression tests or corner case tests!
//	}

//    OperatorPtr LUTMult::parseArguments(Target *target, vector<string> &args) {
//        //int param0, param1;
//        //UserInterface::parseInt(args, "param0", &param0); // param0 has a default value, this method will recover it if it doesnt't find it in args,
//        //UserInterface::parseInt(args, "param1", &param1);
//        //return new LUTMult(target, param0, param1);
//        return new LUTMult(target, 8);
//	}
	
//    void LUTMult::registerFactory(){
//        UserInterface::add("LUTMult", // name
//                                             "My first LUTMult.", // description, string
//											 "NeuralNetworks", // category, from the list defined in UserInterface.cpp
//											 "", //seeAlso
//											 // Now comes the parameter description string.
//											 // Respect its syntax because it will be used to generate the parser and the docs
//											 // Syntax is: a semicolon-separated list of parameterDescription;
//											 // where parameterDescription is parameterName (parameterType)[=defaultValue]: parameterDescriptionString
//											 "param0(int)=16: A first parameter, here used as the input size; \
//                        param1(int): A second parameter, here used as the output size",
//											 // More documentation for the HTML pages. If you want to link to your blog, it is here.
//											 "Feel free to experiment with its code, it will not break anything in FloPoCo. <br> Also see the developper manual in the doc/ directory of FloPoCo.",
//                                             LUTMult::parseArguments
//											 ) ;
//	}

}//namespace
