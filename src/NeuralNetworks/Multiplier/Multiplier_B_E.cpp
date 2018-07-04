// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "../NeuralNetworks/Multiplier/Multiplier_B_E.hpp"
#include "../NeuralNetworks/Multiplier/CompressorTypeB.hpp"
#include "../NeuralNetworks/Multiplier/CompressorTypeE.hpp"
#include <math.h>

using namespace std;
namespace flopoco {




    Multiplier_B_E::Multiplier_B_E(Target* target, int inputWordSize) : Operator(target)
    {
        useNumericStd();
        input_bit_width = inputWordSize;
        int input_1B_witdh=inputWordSize+3;
        int output_1B_witdh=input_1B_witdh+1;

        int input_2E_witdh=output_1B_witdh+3;
        int output_2E_witdh=input_2E_witdh+1;

        srcFileName= join("Multiplier_B_E_", input_bit_width);

        ostringstream name;
        name << "Multiplier_B_E" + to_string(input_bit_width);
        setName(name.str());
        // Copyright
        setCopyrightString("UNIVERSITY of Kassel 2018");

        addInput("X", inputWordSize);
        addInput("conf", 4);
        addOutput("Y",output_2E_witdh);

        declare("sel_1",2);
        declare("sel_2",2);

        vhdl << tab << "sel_1 <= " << "conf(1 downto 0);" << std::endl;
        vhdl << tab << "sel_2 <= " << "conf(3 downto 2);" << std::endl;
        vhdl << std::endl;

        vhdl << tab << declare("B1_1",input_1B_witdh) << " <= std_logic_vector(shift_left(resize(signed(X), B1_1'length),2));" << std::endl;

        vhdl << std::endl;

        CompressorTypeB *myCompB1 = new CompressorTypeB(target, input_1B_witdh,6);
        addToGlobalOpList(myCompB1);

        vhdl << tab << declare("A1_1",input_1B_witdh) << " <= std_logic_vector(shift_left(resize(signed(X), A1_1'length),0));"<< std::endl;
        vhdl << tab << declare("A2_1",input_1B_witdh) << " <= std_logic_vector(shift_left(resize(signed(X), A2_1'length),1));"<< std::endl;
        vhdl << tab << declare("A3_1",input_1B_witdh) << " <= std_logic_vector(shift_left(resize(signed(X), A3_1'length),3));"<< std::endl;

        inPortMapCst(myCompB1, "A1","A1_1");
        inPortMapCst(myCompB1, "A2","A2_1");
        inPortMapCst(myCompB1, "A3","A3_1");
        inPortMapCst(myCompB1, "B1","B1_1");
        inPortMapCst(myCompB1, "S" ,"sel_1");
        outPortMap(myCompB1, "Y", "Y1");
        vhdl << instance(myCompB1,"operation_B_1") << std::endl;
        nextCycle();

        vhdl << tab << declare("A1_2",output_2E_witdh) << " <= std_logic_vector(shift_left(resize(signed(Y1), A1_2'length),0));"<< std::endl;
        vhdl << tab << declare("A2_2",output_2E_witdh) << " <= std_logic_vector(shift_left(resize(signed(Y1), A2_2'length),3));"<< std::endl;
        vhdl << tab << declare("B1_2",output_2E_witdh) << " <= std_logic_vector(shift_left(resize(signed(X), B1_2'length),2));"<< std::endl;

        CompressorTypeE *myCompE2 = new CompressorTypeE(target, output_2E_witdh,6);
        addToGlobalOpList(myCompE2);

        inPortMapCst(myCompE2, "A1","A1_2");
        inPortMapCst(myCompE2, "A2","A2_2");
        inPortMapCst(myCompE2, "B1","B1_2");
        inPortMapCst(myCompE2, "S" ,"sel_2");
        outPortMap(myCompE2, "Y", "Y",false);
        vhdl << instance(myCompE2,"operation_E_2") << std::endl;
        nextCycle();

    }


    OperatorPtr Multiplier_B_E::parseArguments(Target *target, vector<string> &args)
    {
        int param0;
        UserInterface::parseInt(args, "inputWordSize", &param0);
        return new Multiplier_B_E(target, param0);
    }

    void Multiplier_B_E::registerFactory()
    {
        UserInterface::add("Multiplier_B_E", // name
                                     "My first Multiplier_B_E.", // description, string
                                     "NeuralNetworks", // category, from the list defined in UserInterface.cpp
                                     "", //seeAlso
                                     // Now comes the parameter description string.
                                     // Respect its syntax because it will be used to generate the parser and the docs
                                     // Syntax is: a semicolon-separated list of parameterDescription;
                                     // where parameterDescription is parameterName (parameterType)[=defaultValue]: parameterDescriptionString
                                     "inputWordSize(int)=8: Defined the input word size;",
                                     // More documentation for the HTML pages. If you want to link to your blog, it is here.
                                     "Feel free to experiment with its code, it will not break anything in FloPoCo. <br> Also see the developper manual in the doc/ directory of FloPoCo.",
                                     Multiplier_B_E::parseArguments
                                             ) ;
    }



}//namespace
