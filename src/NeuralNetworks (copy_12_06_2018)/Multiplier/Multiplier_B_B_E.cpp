// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "../NeuralNetworks/Multiplier/Multiplier_B_B_E.hpp"
#include "../NeuralNetworks/Multiplier/CompressorTypeB.hpp"
#include "../NeuralNetworks/Multiplier/CompressorTypeE.hpp"
#include <math.h>

using namespace std;
namespace flopoco {




    Multiplier_B_B_E::Multiplier_B_B_E(Target* target, int inputWordSize) : Operator(target)
    {
        input_bit_width = inputWordSize;
        srcFileName= join("Multiplier_B_B_E_", input_bit_width);

        ostringstream name;
        name << "Multiplier_B_B_E" + to_string(input_bit_width);
        setName(name.str());
        // Copyright
        setCopyrightString("UNIVERSITY of Kassel 2018");

        addInput("X", inputWordSize);
        addInput("conf", 6);
        addOutput("Y",inputWordSize+1);

        declare("sel_1",2);
        declare("sel_2",2);
        declare("sel_3",2);

        vhdl << tab << "sel_1 <= " << "conf(1 downto 0)" << std::endl;
        vhdl << tab << "sel_2 <= " << "conf(3 downto 2)" << std::endl;
        vhdl << tab << "sel_3 <= " << "conf(5 downto 4)" << std::endl;

        CompressorTypeB *myCompB1 = new CompressorTypeB(target, inputWordSize,3);
        addToGlobalOpList(myCompB1);

        inPortMapCst(myCompB1, "A1","X");
        inPortMapCst(myCompB1, "A2","X");
        inPortMapCst(myCompB1, "A3","X");
        inPortMapCst(myCompB1, "B0","X");
        inPortMapCst(myCompB1, "S" ,"sel_1");
        outPortMap(myCompB1, "Y", "Y1");
        vhdl << instance(myCompB1,"operation_B_1") << std::endl;
        nextCycle();

        CompressorTypeB *myCompB2 = new CompressorTypeB(target, inputWordSize+1,6);
        addToGlobalOpList(myCompB2);

        inPortMapCst(myCompB2, "A1","Y1");
        inPortMapCst(myCompB2, "A2","Y1");
        inPortMapCst(myCompB2, "A3","Y1");
        inPortMapCst(myCompB2, "B0","Y1");
        inPortMapCst(myCompB2, "S" ,"sel_2");
        outPortMap(myCompB2, "Y", "Y2");
        vhdl << instance(myCompB2,"operation_B_2") << std::endl;
        nextCycle();


        CompressorTypeE *myCompE3 = new CompressorTypeE(target, inputWordSize+2,6);
        addToGlobalOpList(myCompE3);

        inPortMapCst(myCompE3, "A1","Y2");
        inPortMapCst(myCompE3, "A2","Y2");
        inPortMapCst(myCompE3, "A3","Y2");
        inPortMapCst(myCompE3, "B0","Y2");
        inPortMapCst(myCompE3, "S" ,"sel3");
        outPortMap(myCompE3, "Y", "Y",false);
        vhdl << instance(myCompE3,"operation_E_3") << std::endl;
        nextCycle();

    }


    OperatorPtr Multiplier_B_B_E::parseArguments(Target *target, vector<string> &args)
    {
        int param0;
        UserInterface::parseInt(args, "inputWordSize", &param0);
        return new Multiplier_B_B_E(target, param0);
    }

    void Multiplier_B_B_E::registerFactory()
    {
        UserInterface::add("Multiplier_B_B_E", // name
                                     "My first Multiplier_B_B_E.", // description, string
                                     "NeuralNetworks", // category, from the list defined in UserInterface.cpp
                                     "", //seeAlso
                                     // Now comes the parameter description string.
                                     // Respect its syntax because it will be used to generate the parser and the docs
                                     // Syntax is: a semicolon-separated list of parameterDescription;
                                     // where parameterDescription is parameterName (parameterType)[=defaultValue]: parameterDescriptionString
                                     "inputWordSize(int)=8: Defined the input word size;",
                                     // More documentation for the HTML pages. If you want to link to your blog, it is here.
                                     "Feel free to experiment with its code, it will not break anything in FloPoCo. <br> Also see the developper manual in the doc/ directory of FloPoCo.",
                                     Multiplier_B_B_E::parseArguments
                                             ) ;
    }



}//namespace
