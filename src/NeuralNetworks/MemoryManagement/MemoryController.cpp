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




    MemoryController::MemoryController(Target* target, unsigned int dataWidth_, unsigned int addressWidth_, bool outputOverwriteValue, bool needReadCounterReset, bool hasSpecificAddressOnResetPort, bool hasSpecificAddressOnGetNewDataResetPort) :
        Operator(target), dataWidth(dataWidth_), addressWidth(addressWidth_) {

        useNumericStd();

        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="MemoryController";

        // author
        setCopyrightString("Nicolai Fiege, 2017");

        // name of the operator
        ostringstream name;
        name << "MemoryController_dataWidth_" << dataWidth << "_addressWidth_" << addressWidth << ((outputOverwriteValue==true)?("_outputDataToOverwrite"):("")) << ((needReadCounterReset==true)?("_true"):("_false")) << ((hasSpecificAddressOnResetPort==true)?("_true"):("_false")) << ((hasSpecificAddressOnGetNewDataResetPort==true)?("_true"):("_false"));
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

        vhdl << "bramDataWrite <= data_i;" << endl;
        vhdl << "bramWriteEnable <= validData_i;" << endl;

        vhdl << "data_o <= bramDataRead;" << endl;

        if(outputOverwriteValue==true)
        {
                addOutput("data2_o",dataWidth);
                addInput("bramData2Read",dataWidth);
                addOutput("bramAddress2Read",addressWidth+1);
                addInput("getNewData2",1);
                addInput("getNewData2_reset",1);
                vhdl << "data2_o <= bramData2Read;" << endl;
                vhdl << "bramAddress2Read <= " << declare("address2Read",addressWidth+1) << ";" << endl;
                vhdl << "address2Read <= " << "(not read0) & " << declare("addressOverwriteCounter",addressWidth) << ";" << endl;
        }

        if(needReadCounterReset==true)
        {
            addInput("getNewData_reset",1);
        }

        if(hasSpecificAddressOnGetNewDataResetPort==true)
        {
            addInput("getNewData_reset_address",this->addressWidth);
        }

        if(hasSpecificAddressOnResetPort==true)
        {
            addInput("reset_address",this->addressWidth);
        }

        vhdl << "process(clk)" << endl;
        vhdl << "begin" << endl;
        vhdl << tab << "if(rising_edge(clk)) then" << endl;
        vhdl << tab << tab << "if(rst='1') then" << endl;
        vhdl << tab << tab << tab << "validData_o <= \"0\";" << endl;
        vhdl << tab << tab << tab << "addressReadCounter <= (others => '0');" << endl;
        vhdl << tab << tab << tab << "addressWriteCounter <= (others => '0');" << endl;
        if(outputOverwriteValue==true)
        {
            vhdl << tab << tab << tab << "addressOverwriteCounter <= (others => '0');" << endl;
        }
        vhdl << tab << tab << tab << "read0 <= \"0\";" << endl;
        vhdl << tab << tab << "else" << endl;
        vhdl << tab << tab << tab << "if(newStep=\"1\") then" << endl; // read channel
        vhdl << tab << tab << tab << tab << "addressReadCounter <= (others => '0');" << endl;
        vhdl << tab << tab << tab << tab << "validData_o <= \"0\";" << endl;
        if(needReadCounterReset==true)
        {
            vhdl << tab << tab << tab << "elsif(getNewData_reset=\"1\") then" << endl; // read channel
            if(hasSpecificAddressOnGetNewDataResetPort==true)
            {
                vhdl << tab << tab << tab << tab << "addressReadCounter <= getNewData_reset_address;" << endl;
            }
            else
            {
                vhdl << tab << tab << tab << tab << "addressReadCounter <= (others => '0');" << endl;
            }
            vhdl << tab << tab << tab << tab << "validData_o <= \"0\";" << endl;
        }
        vhdl << tab << tab << tab << "elsif(getNewData=\"1\") then" << endl;
        vhdl << tab << tab << tab << tab << "addressReadCounter <= std_logic_vector(unsigned(addressReadCounter)+1);" << endl;
        vhdl << tab << tab << tab << tab << "validData_o <= \"1\";" << endl;
        vhdl << tab << tab << tab << "else" << endl;
        vhdl << tab << tab << tab << tab << "validData_o <= \"0\";" << endl;
        vhdl << tab << tab << tab << "end if;" << endl;
        vhdl << tab << tab << tab << "if(newStep=\"1\") then" << endl;
        vhdl << tab << tab << tab << tab << "read0 <= not read0;" << endl;
        vhdl << tab << tab << tab << "end if;" << endl;
        if(hasSpecificAddressOnResetPort==true)
        {
            vhdl << tab << tab << tab << "if(newStep=\"1\") then" << endl; // write channel
            vhdl << tab << tab << tab << tab << "addressWriteCounter <= (others => '0');" << endl;
            vhdl << tab << tab << tab << tab << "addressOverwriteCounter <= (others => '0');" << endl;
            vhdl << tab << tab << tab << "elsif(getNewData2_reset=\"1\") then" << endl;
            vhdl << tab << tab << tab << tab << "addressWriteCounter <= reset_address;" << endl;
            vhdl << tab << tab << tab << tab << "addressOverwriteCounter <= reset_address;" << endl;
        }
        else
        {
            vhdl << tab << tab << tab << "if(newStep=\"1\"" << ((outputOverwriteValue==true)?(" or getNewData2_reset=\"1\""):("")) << ") then" << endl; // write channel
            vhdl << tab << tab << tab << tab << "addressWriteCounter <= (others => '0');" << endl;
            if(outputOverwriteValue==true)
            {
                vhdl << tab << tab << tab << tab << "addressOverwriteCounter <= (others => '0');" << endl;
            }
        }
        vhdl << tab << tab << tab << "else" << endl;
        vhdl << tab << tab << tab << tab << "if(validData_i=\"1\") then" << endl;
        vhdl << tab << tab << tab << tab << tab << "addressWriteCounter <= std_logic_vector(unsigned(addressWriteCounter)+1);" << endl;
        vhdl << tab << tab << tab << tab << "end if;" << endl;
        if(outputOverwriteValue==true)
        {
            vhdl << tab << tab << tab << tab << "if(getNewData2=\"1\") then" << endl;
            vhdl << tab << tab << tab << tab << tab << "addressOverwriteCounter <= std_logic_vector(unsigned(addressOverwriteCounter)+1);" << endl;
            vhdl << tab << tab << tab << tab << "end if;" << endl;
        }
        vhdl << tab << tab << tab << "end if;" << endl;
        vhdl << tab << tab << "end if;" << endl; //reset
        vhdl << tab << "end if;" << endl; //clk
        vhdl << "end process;" << endl;
    }

}//namespace flopoco
