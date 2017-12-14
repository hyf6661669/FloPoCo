// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

//for reading the .txt file
#include <fstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "Register.hpp"

using namespace std;
namespace flopoco {




    Register::Register(Target* target, unsigned int wordSize_, unsigned int delays_) :
        Operator(target), wordSize(wordSize_), delays(delays_) {

        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="Register";

        // author
        setCopyrightString("Nicolai Fiege, 2017");

        // definition of the name of the operator
        ostringstream name;
        name << "Register_wordSize_" << wordSize << "_delays_" << delays_;
        setName(name.str());

        //add in/out
        addInput("X", wordSize);
        addOutput("R", wordSize);

        for(unsigned int d=0; d<=delays; d++)
        {
            declare("Delay"+to_string(d),wordSize);
        }
        vhdl << "Delays0 <= X;" << endl;

        // DELAY THAT SHIT FOR "DELAYS_" CYCLES!
        if(delays>0)
        {
            vhdl << "process(clk)" << endl
                 << "begin" << endl
                 << tab << "if(rising_edge(clk)) then" << endl
                 << tab << tab << "if(rst='1') then" << endl;
            for(unsigned int d=0; d<=delays; d++)
            {
                vhdl << tab << tab << tab << "Delays" << d << " <= (others => '0');" << endl;
            }
            vhdl << tab << tab << "else" << endl;
            for(unsigned int d=0; d<delays; d++)
            {
                vhdl << tab << tab << tab << "Delays" << d+1 << " <= " << "Delays" << d << ";" << endl;
            }
            vhdl << tab << tab << "end if;" << endl
                 << tab << "end if;" << endl
                 << "end process;" << endl;
        }

		
        vhdl << "R <= Delays" << delays << ";";
    }

}//namespace flopoco
