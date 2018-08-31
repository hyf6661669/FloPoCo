// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "OutputLayer.hpp"

using namespace std;
namespace flopoco {




    OutputLayer::OutputLayer(Target* target, int howMany_, int wordSize_) :
        Layer(target), howMany(howMany_), wordSize(wordSize_) {

		// definition of the source file name, used for info and error reporting using REPORT 
		srcFileName="OutputLayer";

        // author
        setCopyrightString("Nicolai Fiege, 2017");

		// definition of the name of the operator
		ostringstream name;
        name << "OutputLayer";
        setName(name.str());

        addInput("newStep",1);
        for(int i=0; i<howMany; i++)
        {
            addInput("X_"+to_string(i),wordSize);
            addOutput("R_"+to_string(i),wordSize);
            addInput("validData_i_"+to_string(i),1);

            //declare("X_temp_"+to_string(i),wordSize);

            this->vhdl << tab << "Feature" << i << " : process(clk)" << endl
                       << tab << "begin" << endl
                       << tab << tab << "if(rising_edge(clk)) then" << endl
                       << tab << tab << tab << "if(rst='1' or newStep=\"1\") then" << endl
                       << tab << tab << tab << tab << "R_" << i << " <= (others => '0');" << endl
                       << tab << tab << tab << "else" << endl
                       << tab << tab << tab << tab << "if(validData_i_" << i << "=\"1\") then" << endl
                       << tab << tab << tab << tab << tab << "R_" << i << " <= X_" << i << ";" << endl
                       << tab << tab << tab << tab << "end if;" << endl
                       << tab << tab << tab << "end if;" << endl
                       << tab << tab << "end if;" << endl
                       << tab << "end process;" << endl;
        }
    }

    string OutputLayer::getOutputSignalName(int feature)
    {
        return "R_"+to_string(feature);
    }

    string OutputLayer::getInputSignalName(int feature)
    {
        return "X_"+to_string(feature);
    }

}//namespace flopoco
