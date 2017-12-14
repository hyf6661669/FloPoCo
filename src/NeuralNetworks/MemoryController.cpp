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
#include "MemoryController.hpp"

using namespace std;
namespace flopoco {




    MemoryController::MemoryController(Target* target, unsigned int dataWidth_, unsigned int addressWidth_) :
        Operator(target), dataWidth(dataWidth_), addressWidth(addressWidth_) {

        useNumericStd();

        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="MemoryController";

        // author
        setCopyrightString("Nicolai Fiege, 2017");

        // definition of the name of the operator
        ostringstream name;
        name << "MemoryController_dataWidth_" << dataWidth << "_addressWidth_" << addressWidth;
        setName(name.str());

        // add in/out
        // communication between layers
        addInput("newStep",1);
        addInput("validData_i",1);
        addInput("getNewData",1);
        addInput("data_i",dataWidth);
        addOutput("validData_o",1);
        addOutput("data_o",dataWidth);

        // communication with bram
        addInput("bramDataRead",dataWidth);
        addOutput("bramDataWrite",dataWidth);
        addOutput("bramAddressRead",addressWidth+1);
        addOutput("bramAddressWrite",addressWidth+1);
        addOutput("bramWriteEnable",1);


        // vhdl
        vhdl << "bramAddressRead <= " << declare("addressRead",addressWidth+1) << ";" << endl;
        vhdl << "bramAddressWrite <= " << declare("addressWrite",addressWidth+1) << ";" << endl;
        declare("read0",1);

        vhdl << "addressRead <= " << "read0 & " << declare("addressReadCounter",addressWidth) << ";" << endl;
        vhdl << "addressWrite <= " << "(not read0) & " << declare("addressWriteCounter",addressWidth) << ";" << endl;

        vhdl << "process(clk,rst)" << endl
             << "begin" << endl
             << tab << "if(rising_edge(clk)) then" << endl
             << tab << tab << "if(rst='1') then" << endl
             << tab << tab << tab << "validData_o <= \"0\";" << endl
             << tab << tab << tab << "data_o <= (others => '0');" << endl
             << tab << tab << tab << "bramDataWrite <= (others => '0');" << endl
             << tab << tab << tab << "addressReadCounter <= (others => '0');" << endl
             << tab << tab << tab << "addressWriteCounter <= (others => '0');" << endl
             << tab << tab << tab << "read0 <= \"0\";" << endl
             << tab << tab << "else" << endl
             << tab << tab << tab << "if(validData_i=\"1\") then" << endl
             << tab << tab << tab << tab << "bramDataWrite <= data_i;" << endl
             << tab << tab << tab << tab << "bramWriteEnable <= \"1\";" << endl
             << tab << tab << tab << tab << "addressWriteCounter <= std_logic_vector(unsigned(addressWriteCounter)+1);" << endl
             << tab << tab << tab << "else" << endl
             << tab << tab << tab << tab << "bramWriteEnable <= \"1\";" << endl
             << tab << tab << tab << "end if;" << endl
             << tab << tab << tab << "if(getNewData=\"1\") then" << endl
             << tab << tab << tab << tab << "data_o <= bramDataRead;" << endl
             << tab << tab << tab << tab << "validData_o <= \"1\";" << endl
             << tab << tab << tab << tab << "addressReadCounter <= std_logic_vector(unsigned(addressReadCounter)+1);" << endl
             << tab << tab << tab << "else" << endl
             << tab << tab << tab << tab << "validData_o <= \"0\";" << endl
             << tab << tab << tab << "end if;" << endl
             << tab << tab << tab << "if(newStep=\"1\") then" << endl
             << tab << tab << tab << tab << "addressReadCounter <= (others => '0');" << endl
             << tab << tab << tab << tab << "addressWriteCounter <= (others => '0');" << endl
             << tab << tab << tab << tab << "read0 <= not read0;" << endl
             << tab << tab << tab << "end if;" << endl
             << tab << tab << "end if;" << endl
             << tab << "end if;" << endl
             << "end process;" << endl;
    }

}//namespace flopoco
