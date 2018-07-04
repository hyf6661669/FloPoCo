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
        useNumericStd();
        input_bit_width = inputWordSize;
        int input_1B_witdh=inputWordSize+2;
        int output_1B_witdh=input_1B_witdh+1;

        int input_2B_witdh=output_1B_witdh+2;
        int output_2B_witdh=input_2B_witdh+1;

        int input_3E_witdh=output_2B_witdh+2;
        int output_3E_witdh=input_3E_witdh+1;

        srcFileName= join("Multiplier_B_B_E_", input_bit_width);

        ostringstream name;
        name << "Multiplier_B_B_E" + to_string(input_bit_width);
        setName(name.str());
        // Copyright
        setCopyrightString("UNIVERSITY of Kassel 2018");

        addInput("X", inputWordSize);
        addInput("conf", 6);
        addOutput("Y",output_3E_witdh);

        declare("sel_1",2);
        declare("sel_2",2);
        declare("sel_3",2);

        vhdl << tab << "sel_1 <= " << "conf(1 downto 0);" << std::endl;
        vhdl << tab << "sel_2 <= " << "conf(3 downto 2);" << std::endl;
        vhdl << tab << "sel_3 <= " << "conf(5 downto 4);" << std::endl;
        vhdl << std::endl;

        vhdl << tab << declare("X0",input_1B_witdh) << " <= std_logic_vector(resize(signed(X), X0'length));" << std::endl;
        vhdl << tab << declare("X1",input_2B_witdh) << " <= std_logic_vector(resize(signed(X), X1'length));"<< std::endl;
        vhdl << tab << declare("X2",input_3E_witdh) << " <= std_logic_vector(resize(signed(X), X2'length));"<< std::endl;

        vhdl << std::endl;

        CompressorTypeB *myCompB1 = new CompressorTypeB(target, input_1B_witdh,3);
        addToGlobalOpList(myCompB1);

        vhdl << tab << declare("X_shifted_1",input_1B_witdh) << " <= std_logic_vector(shift_left(resize(signed(X), X_shifted_1'length),0));"<< std::endl;
        vhdl << tab << declare("X_shifted_2",input_1B_witdh) << " <= std_logic_vector(shift_left(resize(signed(X), X_shifted_2'length),1));"<< std::endl;
        vhdl << tab << declare("X_shifted_3",input_1B_witdh) << " <= std_logic_vector(shift_left(resize(signed(X), X_shifted_3'length),2));"<< std::endl;

        inPortMapCst(myCompB1, "A1","X_shifted_1");
        inPortMapCst(myCompB1, "A2","X_shifted_2");
        inPortMapCst(myCompB1, "A3","X_shifted_3");
        inPortMapCst(myCompB1, "B1","X0");
        inPortMapCst(myCompB1, "S" ,"sel_1");
        outPortMap(myCompB1, "Y", "Y1");
        vhdl << instance(myCompB1,"operation_B_1") << std::endl;
        nextCycle();

        CompressorTypeB *myCompB2 = new CompressorTypeB(target, input_2B_witdh,6);
        addToGlobalOpList(myCompB2);

        vhdl << tab << declare("Y1_shifted_1",input_2B_witdh) << " <= std_logic_vector(shift_left(resize(signed(Y1), Y1_shifted_1'length),0));"<< std::endl;
        vhdl << tab << declare("Y1_shifted_2",input_2B_witdh) << " <= std_logic_vector(shift_left(resize(signed(Y1), Y1_shifted_2'length),1));"<< std::endl;
        vhdl << tab << declare("Y1_shifted_3",input_2B_witdh) << " <= std_logic_vector(shift_left(resize(signed(Y1), Y1_shifted_3'length),2));"<< std::endl;

        inPortMapCst(myCompB2, "A1","Y1_shifted_1");
        inPortMapCst(myCompB2, "A2","Y1_shifted_2");
        inPortMapCst(myCompB2, "A3","Y1_shifted_3");
        inPortMapCst(myCompB2, "B1","X1");
        inPortMapCst(myCompB2, "S" ,"sel_2");
        outPortMap(myCompB2, "Y", "Y2");
        vhdl << instance(myCompB2,"operation_B_2") << std::endl;
        nextCycle();

        vhdl << tab << declare("Y2_shifted_1",input_3E_witdh) << " <= std_logic_vector(shift_left(resize(signed(Y2), Y2_shifted_1'length),0));"<< std::endl;
        vhdl << tab << declare("Y2_shifted_2",input_3E_witdh) << " <= std_logic_vector(shift_left(resize(signed(Y2), Y2_shifted_1'length),1));"<< std::endl;

        CompressorTypeE *myCompE3 = new CompressorTypeE(target, input_3E_witdh,6);
        addToGlobalOpList(myCompE3);

        inPortMapCst(myCompE3, "A1","Y2_shifted_1");
        inPortMapCst(myCompE3, "A2","Y2_shifted_2");
        inPortMapCst(myCompE3, "B1","X2");
        inPortMapCst(myCompE3, "S" ,"sel_3");
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
