// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "GlobalController.hpp"

using namespace std;
namespace flopoco {




    GlobalController::GlobalController(Target* target, unsigned int numberOfInputs_) :
        Operator(target), numberOfInputs(numberOfInputs_) {

        this->useNumericStd();

        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="GlobalController";

        // author
        setCopyrightString("Nicolai Fiege, 2017");

        // definition of the name of the operator
        ostringstream name;
        name << "GlobalController";
        setName(name.str());

        //add in/out
        addOutput("newStep",1);

        for(unsigned int i=0; i<numberOfInputs; i++)
        {
            addInput("finished_"+to_string(i),1);
        }

        // signal to get rising edge of all finished-signals
        vhdl << declare("finished_all",1) << " <= ";
        for(unsigned int i=0; i<numberOfInputs; i++)
        {
            if(i>0)
            {
                vhdl << " and ";
            }
            vhdl << "finished_" + to_string(i);
        }
        vhdl << ";" << endl;

        // process for newStep and to detect rising edge of 'finished_all'
        vhdl << "process(clk)" << endl
             << "begin" << endl
             << tab << "if(rising_edge(clk)) then" << endl
             << tab << tab << "if(rst='1') then" << endl
             << tab << tab << tab << declare("finished_all_last",1) << " <= \"0\";" << endl
             << tab << tab << tab << "newStep <= \"0\";" << endl
             << tab << tab << "else" << endl
             << tab << tab << tab << "finished_all_last <= finished_all;" << endl
             << tab << tab << tab << "if(finished_all_last=\"0\" and finished_all=\"1\") then" << endl
             << tab << tab << tab << tab << "newStep <= \"1\";" << endl
             << tab << tab << tab << "else" << endl
             << tab << tab << tab << tab << "newStep <= \"0\";" << endl
             << tab << tab << tab << "end if;" << endl
             << tab << tab << "end if;" << endl
             << tab << "end if;" << endl
             << "end process;" << endl;
    }

}//namespace flopoco
