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
#include "ModuloCounter.hpp"

using namespace std;
namespace flopoco {




    ModuloCounter::ModuloCounter(Target* target, unsigned int mod_) :
        Operator(target), mod(mod_) {

        useNumericStd();

        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="ModuloCounter";

        // author
        setCopyrightString("Nicolai Fiege, 2017");

        // definition of the name of the operator
        ostringstream name;
        name << "ModuloCounter_" << mod;
        setName(name.str());
		
        unsigned int wordSize=ceil(log2(mod));

        //add in/out
        addInput("enable", 1);
        addOutput("counter", wordSize);
        addInput("manualReset", 1);

        vhdl << "counter <= " << declare("counter_temp", wordSize) << ";" << endl;
		
        vhdl << "process(clk,rst,manualReset)" << endl
             << "begin" << endl
             << tab << "if(rst='1' or manualReset=\"1\") then" << endl
             << tab << tab << "counter_temp <= (others => '0');" << endl
             << tab << "else" << endl
             << tab << tab << "if(rising_edge(clk)) then" << endl
             << tab << tab << tab << "if(enable=\"1\") then" << endl
             << tab << tab << tab << tab << "if(unsigned(counter_temp)<" << mod-1 << ") then" << endl
             << tab << tab << tab << tab << tab << "counter_temp <= std_logic_vector(unsigned(counter_temp)+1);" << endl
             << tab << tab << tab << tab << "else" << endl
             << tab << tab << tab << tab << tab << "counter_temp <= (others => '0');" << endl
             << tab << tab << tab << tab << "end if;" << endl
             << tab << tab << tab << "end if;" << endl
             << tab << tab << "end if;" << endl
             << tab << "end if;" << endl
             << "end process;" << endl;
		
    }

}//namespace flopoco
