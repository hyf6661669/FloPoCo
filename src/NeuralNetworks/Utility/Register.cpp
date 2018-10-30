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




    Register::Register(Target* target, unsigned int wordSize_, unsigned int delays_, bool enabled, bool externalReset, bool syncReset) :
        Operator(target), wordSize(wordSize_), delays(delays_) {

        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="Register";

        // author
        setCopyrightString("Nicolai Fiege, 2017");

        // definition of the name of the operator
        ostringstream name;
        name << "Register_wordSize_" << wordSize << "_delays_" << delays_ << "_isEnabled_" << enabled << "_externalReset_" << externalReset << (syncReset==true?"_sync":"_async");
        setName(name.str());

        //add in/out
        addInput("X", wordSize);
        addOutput("R", wordSize);
        if(enabled==true)
        {
            addInput("E", 1);
        }
        if(externalReset==true)
        {
            addInput("ExternalReset", 1);
        }

        for(unsigned int d=0; d<=delays; d++)
        {
            declare("Delay"+to_string(d),wordSize);
        }
        vhdl << "Delay0 <= X;" << endl;


        if(delays>0)
        {
            vhdl << "process(clk" << (syncReset==true?"":",rst") << ((externalReset==true && syncReset==false)?",ExternalReset":"") << ")" << endl;
            vhdl << "begin" << endl;
            if(syncReset==true)
            {
                vhdl << tab << "if(rising_edge(clk)) then" << endl;
                vhdl << tab << tab << "if(rst='1'" << ((externalReset==true)?(" or ExternalReset=\"1\""):("")) << ") then" << endl;
                for(unsigned int d=1; d<=delays; d++)
                {
                    vhdl << tab << tab << tab << "Delay" << d << " <= (others => '0');" << endl;
                }
                if(enabled==false)
                {
                    vhdl << tab << tab << "else" << endl;
                }
                else
                {
                    vhdl << tab << tab << "elsif(E=\"1\") then" << endl;
                }
                for (unsigned int d = 0; d < delays; d++)
                {
                    vhdl << tab << tab << tab << "Delay" << d + 1 << " <= " << "Delay" << d << ";" << endl;
                }
                vhdl << tab << tab << "end if;" << endl;
                vhdl << tab << "end if;" << endl;
            }
            else
            {
                vhdl << tab << "if(rst='1'" << ((externalReset==true)?(" or ExternalReset=\"1\""):("")) << ") then" << endl;
                for(unsigned int d=1; d<=delays; d++)
                {
                    vhdl << tab << tab << "Delay" << d << " <= (others => '0');" << endl;
                }
                vhdl << tab << "elsif(rising_edge(clk)) then" << endl;
                if(enabled==false)
                {
                    for (unsigned int d = 0; d < delays; d++)
                    {
                        vhdl << tab << tab << "Delay" << d + 1 << " <= " << "Delay" << d << ";" << endl;
                    }
                }
                else
                {
                    vhdl << tab << tab << "if(E=\"1\") then" << endl;
                    for (unsigned int d = 0; d < delays; d++)
                    {
                        vhdl << tab << tab << tab << "Delay" << d + 1 << " <= " << "Delay" << d << ";" << endl;
                    }
                    vhdl << tab << tab << "end if;" << endl;
                }
                vhdl << tab << "end if;" << endl;
            }

            vhdl << "end process;" << endl;
        }

		
        vhdl << "R <= Delay" << delays << ";" << endl;
    }

}//namespace flopoco
