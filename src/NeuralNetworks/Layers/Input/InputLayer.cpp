// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "InputLayer.hpp"
#include "NeuralNetworks/MemoryManagement/BlockRam.hpp"

using namespace std;
namespace flopoco {




    InputLayer::InputLayer(Target* target, int numberOfFeatures_, int wordSize_, int bramAddressWidth_) :
        Layer(target), numberOfFeatures(numberOfFeatures_), wordSize(wordSize_), bramAddressWidth(bramAddressWidth_) {

        useNumericStd();

		// definition of the source file name, used for info and error reporting using REPORT 
		srcFileName="InputLayer";

        // author
        setCopyrightString("Nicolai Fiege, 2017");

		// definition of the name of the operator
		ostringstream name;
        name << "InputLayer";
        setName(name.str());
		
        addInput("newStep",1);
        stringstream typeStream;
        typeStream << "(waiting, writing);" << endl;

        for(int i=0;i<numberOfFeatures;i++)
        {
            // write-side
            addInput("X_"+to_string(i),wordSize);
            addInput("validData_i_"+to_string(i),1);
            addInput("newDataSet_"+to_string(i),1);

            // read-side
            addOutput("R_"+to_string(i),wordSize);
            addOutput("validData_o_"+to_string(i),1);
            addInput("getNewData_"+to_string(i),1);

            // flag signal to decide which memory area to read and which to write
            declare("writeFlag_"+to_string(i),1);

            // small FSM to wait for a valid input frame
            typeStream << "signal state" << i << ": state_t;" << endl;

            // one BRAM
            BlockRam* bram = new BlockRam(target,wordSize,bramAddressWidth);
            addSubComponent(bram);
            inPortMap(bram,"Address_W_in",declare("WriteAddress_"+to_string(i),bramAddressWidth+1));
            inPortMap(bram,"Address_R_in",declare("ReadAddress_"+to_string(i),bramAddressWidth+1));
            inPortMap(bram,"Data_in","X_"+to_string(i));
            inPortMap(bram,"WriteEnable","validData_i_"+to_string(i));
            outPortMap(bram,"Data_out",declare("bramRead_"+to_string(i)));
            vhdl << instance(bram,"BlockRam_"+to_string(i));
        }
        addType("state_t", typeStream.str());

        // FSM to wait for a valid input frame
        vhdl << "process(clk)" << endl
             << "begin" << endl
             << tab << "if(rising_edge(clk)) then" << endl
             << tab << tab << "if(rst='1' or newStep=\"1\") then" << endl;
        for(int i=0;i<numberOfFeatures;i++)
        {
            vhdl << tab << tab << tab << "state" << i << " <= waiting;" << endl;
        }
        vhdl << tab << tab << "else" << endl;
        for(int i=0;i<numberOfFeatures;i++)
        {
            vhdl << tab << tab << tab << "if(state" << i << "=waiting and newDataSet_" << i << "=\"1\") then" << endl
                 << tab << tab << tab << tab << "state" << i << " <= writing;" << endl
                 << tab << tab << tab << "end if;" << endl;
        }
        vhdl << tab << tab << "end if;" << endl
             << tab << "end if;" << endl
             << "end process;" << endl;

        // Flag Signal to decide which Addresses to read and which to write
        for(unsigned int i=0;i<numberOfFeatures;i++)
        {
            declare("addressReadCounter_"+to_string(i),bramAddressWidth);
            declare("addressWriteCounter_"+to_string(i),bramAddressWidth);
        }

        vhdl << "process(clk)" << endl
             << "begin" << endl
             << tab << "if(rising_edge(clk)) then" << endl
             << tab << tab << "if(rst='1') then" << endl;
        for(unsigned int i=0;i<numberOfFeatures;i++)
        {
            vhdl << tab << tab << tab << "addressReadCounter_" << i << " <= (others => '0');" << endl
                 << tab << tab << tab << "addressWriteCounter_" << i << " <= (others => '0');" << endl
                 << tab << tab << tab << "validData_o_" << i << " <= \"0\";" << endl;
        }
        vhdl << tab << tab << "else" << endl
             << tab << tab << tab << "if(newStep=\"1\") then" << endl;
        for(unsigned int i=0;i<numberOfFeatures;i++)
        {
            vhdl << tab << tab << tab << tab << "addressReadCounter_" << i << " <= (others => '0');" << endl
                 << tab << tab << tab << tab << "validData_o_" << i << " <= \"0\";" << endl;
        }
        vhdl << tab << tab << tab << "else" << endl;
        for(unsigned int i=0;i<numberOfFeatures;i++)
        {
            vhdl << tab << tab << tab << tab << "if(getNewData_" << i << "=\"1\") then" << endl
                 << tab << tab << tab << tab << tab << "addressReadCounter_" << i << " <= std_logic_vector(unsigned(addressReadCounter_" << i << ")+1);" << endl
                 << tab << tab << tab << tab << tab << "validData_o_" << i << " <= \"1\";" << endl
                 << tab << tab << tab << tab << "end if;" << endl;
        }
        vhdl << tab << tab << tab << "end if;" << endl;
        for(unsigned int i=0;i<numberOfFeatures;i++)
        {
            vhdl << tab << tab << tab << "if(newDataSet_" << i << "=\"1\") then" << endl
                 << tab << tab << tab << tab << "addressWriteCounter_" << i << " <= (others => '0');" << endl
                 << tab << tab << tab << "else" << endl
                 << tab << tab << tab << tab << "if(validData_i_" << i << "=\"1\" and state=writing) then" << endl
                 << tab << tab << tab << tab << tab << "addressWriteCounter_" << i << " <= std_logic_vector(unsigned(addressWriteCounter_" << i << ")+1);" << endl
                 << tab << tab << tab << tab << "end if;" << endl
                 << tab << tab << tab << "end if;" << endl;
        }
        vhdl << tab << tab << "end if;" << endl
             << tab << "end if;" << endl
             << "end process;" << endl;

        for(unsigned int i=0;i<numberOfFeatures;i++)
        {
            vhdl << "WriteAddress_" << i << " <= writeFlag & addressWriteCounter_" << i << ";" << endl;
            vhdl << "ReadAddress_" << i << " <= (not writeFlag) & addressReadCounter_" << i << ";" << endl;
        }

        // Input/Output

    }

    string InputLayer::getOutputSignalName(int feature)
    {
        return "R_"+to_string(feature);
    }

    string InputLayer::getInputSignalName(int feature)
    {
        return "X_"+to_string(feature);
    }

}//namespace flopoco
