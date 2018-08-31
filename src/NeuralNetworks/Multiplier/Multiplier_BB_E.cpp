// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "../NeuralNetworks/Multiplier/Multiplier_BB_E.hpp"
#include "../NeuralNetworks/Multiplier/CompressorTypeB.hpp"
#include "../NeuralNetworks/Multiplier/CompressorTypeE.hpp"
#include <math.h>

using namespace std;
namespace flopoco {




    Multiplier_BB_E::Multiplier_BB_E(Target* target, int inputWordSize) : Operator(target)
    {
        useNumericStd();
        input_bit_width = inputWordSize;

        int io_1B1_witdh = input_bit_width+1+3;
        int io_1B2_witdh = input_bit_width+1+3;

        int io_2E1_witdh = io_1B1_witdh+1+3;

        srcFileName= join("Multiplier_BB_E", input_bit_width);

        ostringstream name;
        name << "Multiplier_BB_E" + to_string(input_bit_width);
        setName(name.str());
        // Copyright
        setCopyrightString("UNIVERSITY of Kassel 2018");

        addInput("X", inputWordSize);
        addInput("conf",6);
        addOutput("Y",io_2E1_witdh);

        declare("sel_1",2);
        declare("sel_2",2);
        declare("sel_3",2);


        vhdl << tab << "sel_1 <= " << "conf(1 downto 0);" << std::endl;
        vhdl << tab << "sel_2 <= " << "conf(3 downto 2);" << std::endl;
        vhdl << tab << "sel_3 <= " << "conf(5 downto 4);" << std::endl;
        vhdl << std::endl;

        vhdl << tab << declare("X0",io_1B1_witdh) << " <= std_logic_vector(resize(signed(X), X0'length));" << std::endl;
        vhdl << tab << declare("X1",io_1B2_witdh) << " <= std_logic_vector(resize(signed(X), X1'length));"<< std::endl;
        vhdl << std::endl;


        CompressorTypeB *myCompB1 = new CompressorTypeB(target, io_1B1_witdh,4);
        addToGlobalOpList(myCompB1);

        vhdl << tab << declare("X1_shifted_1",io_1B1_witdh) << " <= std_logic_vector(shift_left(resize(signed(X), X1_shifted_1'length),0));"<< std::endl;
        vhdl << tab << declare("X1_shifted_2",io_1B1_witdh) << " <= std_logic_vector(shift_left(resize(signed(X), X1_shifted_2'length),2));"<< std::endl;
        vhdl << tab << declare("X1_shifted_3",io_1B1_witdh) << " <= std_logic_vector(shift_left(resize(signed(X), X1_shifted_3'length),3));"<< std::endl;

        vhdl << tab << declare("X1_shifted_4",io_1B1_witdh) << " <= std_logic_vector(shift_left(resize(signed(X), X1_shifted_4'length),3));"<< std::endl;

        inPortMapCst(myCompB1, "A1","X1_shifted_1");
        inPortMapCst(myCompB1, "A2","X1_shifted_2");
        inPortMapCst(myCompB1, "A3","X1_shifted_3");
        inPortMapCst(myCompB1, "B1","X1_shifted_4");
        inPortMapCst(myCompB1, "S" ,"sel_1");
        outPortMap(myCompB1, "Y", "Y1");
        vhdl << instance(myCompB1,"operation_B_1") << std::endl;

        CompressorTypeB *myCompB2 = new CompressorTypeB(target, io_1B2_witdh,3);
        addToGlobalOpList(myCompB2);

        vhdl << tab << declare("X2_shifted_1",io_1B2_witdh) << " <= std_logic_vector(shift_left(resize(signed(X), X2_shifted_1'length),0));"<< std::endl;
        vhdl << tab << declare("X2_shifted_2",io_1B2_witdh) << " <= std_logic_vector(shift_left(resize(signed(X), X2_shifted_2'length),1));"<< std::endl;
        vhdl << tab << declare("X2_shifted_3",io_1B2_witdh) << " <= std_logic_vector(shift_left(resize(signed(X), X2_shifted_3'length),3));"<< std::endl;

        vhdl << tab << declare("X2_shifted_4",io_1B2_witdh) << " <= std_logic_vector(shift_left(resize(signed(X), X2_shifted_4'length),0));"<< std::endl;

        inPortMapCst(myCompB2, "A1","X2_shifted_1");
        inPortMapCst(myCompB2, "A2","X2_shifted_2");
        inPortMapCst(myCompB2, "A3","X2_shifted_3");
        inPortMapCst(myCompB2, "B1","X2_shifted_4");
        inPortMapCst(myCompB2, "S" ,"sel_2");
        outPortMap(myCompB2, "Y", "Y2");
        vhdl << instance(myCompB2,"operation_B_2") << std::endl;

        nextCycle();

        vhdl << tab << declare("Y1_shifted_1", io_2E1_witdh) << " <= std_logic_vector(shift_left(resize(signed(Y1), Y1_shifted_1'length),0));"<< std::endl;
        vhdl << tab << declare("Y1_shifted_2", io_2E1_witdh) << " <= std_logic_vector(shift_left(resize(signed(Y1), Y1_shifted_2'length),3));"<< std::endl;

        vhdl << tab << declare("Y2_shifted_3", io_2E1_witdh) << " <= std_logic_vector(shift_left(resize(signed(Y2), Y2_shifted_3'length),0));"<< std::endl;

        CompressorTypeE *myCompE3 = new CompressorTypeE(target, io_2E1_witdh,6);
        addToGlobalOpList(myCompE3);

        inPortMapCst(myCompE3, "A1","Y1_shifted_1");
        inPortMapCst(myCompE3, "A2","Y1_shifted_2");
        inPortMapCst(myCompE3, "B1","Y2_shifted_3");
        inPortMapCst(myCompE3, "S" ,"sel_3");
        outPortMap(myCompE3, "Y", "Y",false);
        vhdl << instance(myCompE3,"operation_E_3") << std::endl;
        nextCycle();

    }


    OperatorPtr Multiplier_BB_E::parseArguments(Target *target, vector<string> &args)
    {
        int param0;
        UserInterface::parseInt(args, "inputWordSize", &param0);
        return new Multiplier_BB_E(target, param0);
    }

    void Multiplier_BB_E::registerFactory()
    {
        UserInterface::add("Multiplier_BB_E", // name
                                     "My first Multiplier_BB_E.", // description, string
                                     "NeuralNetworks", // category, from the list defined in UserInterface.cpp
                                     "", //seeAlso
                                     // Now comes the parameter description string.
                                     // Respect its syntax because it will be used to generate the parser and the docs
                                     // Syntax is: a semicolon-separated list of parameterDescription;
                                     // where parameterDescription is parameterName (parameterType)[=defaultValue]: parameterDescriptionString
                                     "inputWordSize(int)=8: Defined the input word size;",
                                     // More documentation for the HTML pages. If you want to link to your blog, it is here.
                                     "Feel free to experiment with its code, it will not break anything in FloPoCo. <br> Also see the developper manual in the doc/ directory of FloPoCo.",
                                     Multiplier_BB_E::parseArguments
                                             ) ;
    }



}//namespace
