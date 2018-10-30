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




    InputLayer::InputLayer(Target* target, NeuralNetwork* parent, int numberOfFeatures_, int wordSize_, int bramAddressWidth_, bool needGetNewDataReset) :
        Layer(target,parent), numberOfFeatures(numberOfFeatures_), wordSize(wordSize_), bramAddressWidth(bramAddressWidth_) {

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
        typeStream << "(waiting, reading);" << endl;

        for(int i=0;i<numberOfFeatures;i++)
        {
            // write-side
            addInput("X_"+to_string(i),wordSize);
            addInput(Layer::getValidDataName(i)+"_i",1);
            addInput("newDataSet_"+to_string(i),1);

            // read-side
            addOutput("R_"+to_string(i),wordSize);
            addOutput(Layer::getValidDataName(i)+"_o",1);
            addInput(Layer::getGetNewDataName(i),1);

            if(needGetNewDataReset==true)
            {
                addInput(Layer::getGetNewDataName(i)+"_reset",1);
            }

            // flag signal to decide which memory area to read and which to write
            declare("writeFlag_"+to_string(i),1);

            // small FSM to wait for a valid input frame
            typeStream << "signal state_" << i << ": state_t";
            if(i<numberOfFeatures-1)
            {
                typeStream << ";" << endl;
            }

            // one BRAM
            BlockRam* bram = new BlockRam(target,wordSize,bramAddressWidth+1);
            addSubComponent(bram);
            inPortMap(bram,"Address_W_in",declare("WriteAddress_"+to_string(i),bramAddressWidth+1));
            inPortMap(bram,"Address_R_in",declare("ReadAddress_"+to_string(i),bramAddressWidth+1));
            inPortMap(bram,"Data_in","X_"+to_string(i));
            inPortMap(bram,"WriteEnable",Layer::getValidDataName(i)+"_i");
            outPortMap(bram,"Data_out","R_"+to_string(i),false);
            //outPortMap(bram,"Data_out","R_temp_"+to_string(i),true);
            vhdl << instance(bram,"BlockRam_"+to_string(i));
        }
        addType("state_t", typeStream.str());

        // FSM to wait for a valid input frame
        vhdl << "process(clk)" << endl
             << "begin" << endl
             << tab << "if(rising_edge(clk)) then" << endl
             << tab << tab << "if(rst='1') then" << endl;
        for(int i=0;i<numberOfFeatures;i++)
        {
            vhdl << tab << tab << tab << "state_" << i << " <= waiting;" << endl;
            vhdl << tab << tab << tab << "writeFlag_" << i << " <= \"0\";" << endl;
        }
        vhdl << tab << tab << "elsif(newStep=\"1\") then" << endl;
        for(int i=0;i<numberOfFeatures;i++)
        {
            vhdl << tab << tab << tab << "state_" << i << " <= waiting;" << endl;
        }
        vhdl << tab << tab << "else" << endl;
        for(int i=0;i<numberOfFeatures;i++)
        {
            vhdl << tab << tab << tab << "if(state_" << i << "=waiting and newDataSet_" << i << "=\"1\") then" << endl
                 << tab << tab << tab << tab << "state_" << i << " <= reading;" << endl
                 << tab << tab << tab << tab << "writeFlag_" << i << " <= not writeFlag_" << i << ";" << endl
                 << tab << tab << tab << "end if;" << endl;
        }
        vhdl << tab << tab << "end if;" << endl
             << tab << "end if;" << endl
             << "end process;" << endl;


        // Flag Signal to decide which Addresses to read and which to write
        for(int i=0;i<numberOfFeatures;i++)
        {
            declare("addressReadCounter_"+to_string(i),bramAddressWidth);
            declare("addressWriteCounter_"+to_string(i),bramAddressWidth);
        }

        vhdl << "process(clk)" << endl
             << "begin" << endl
             << tab << "if(rising_edge(clk)) then" << endl
             << tab << tab << "if(rst='1') then" << endl;
        for(int i=0;i<numberOfFeatures;i++)
        {
            vhdl << tab << tab << tab << "addressReadCounter_" << i << " <= (others => '0');" << endl
                 << tab << tab << tab << "addressWriteCounter_" << i << " <= (others => '0');" << endl
                 //<< tab << tab << tab << "R_" << i << " <= (others => '0');" << endl
                 << tab << tab << tab << Layer::getValidDataName(i)+"_o" << " <= \"0\";" << endl;
        }
        vhdl << tab << tab << "else" << endl
             << tab << tab << tab << "if(newStep=\"1\") then" << endl;
        for(int i=0;i<numberOfFeatures;i++)
        {
            vhdl << tab << tab << tab << tab << "addressReadCounter_" << i << " <= (others => '0');" << endl
                 //<< tab << tab << tab << tab << "R_" << i << " <= (others => '0');" << endl
                 << tab << tab << tab << tab << Layer::getValidDataName(i)+"_o" << " <= \"0\";" << endl;
        }
        vhdl << tab << tab << tab << "else" << endl;
        for(int i=0;i<numberOfFeatures;i++)
        {
            vhdl << tab << tab << tab << tab << "if(" << Layer::getGetNewDataName(i) << "=\"1\" and state_" << i << "=reading) then" << endl;
            vhdl << tab << tab << tab << tab << tab << "addressReadCounter_" << i << " <= std_logic_vector(unsigned(addressReadCounter_" << i << ")+1);" << endl;
            vhdl << tab << tab << tab << tab << tab << Layer::getValidDataName(i)+"_o" << " <= \"1\";" << endl;
                 //<< tab << tab << tab << tab << tab << "R_" << i << " <= R_temp_" << i << ";" << endl
            if(needGetNewDataReset==true)
            {
                vhdl << tab << tab << tab << tab << "elsif(" << Layer::getGetNewDataName(i) << "_reset=\"1\") then" << endl;
                vhdl << tab << tab << tab << tab << tab << "addressReadCounter_" << i << " <= (others => '0');" << endl;
                vhdl << tab << tab << tab << tab << tab << Layer::getValidDataName(i)+"_o" << " <= \"0\";" << endl;
            }
            vhdl << tab << tab << tab << tab << "else" << endl;
            vhdl << tab << tab << tab << tab << tab << Layer::getValidDataName(i)+"_o" << " <= \"0\";" << endl;
            vhdl << tab << tab << tab << tab << "end if;" << endl;
        }
        vhdl << tab << tab << tab << "end if;" << endl;
        for(int i=0;i<numberOfFeatures;i++)
        {
            vhdl << tab << tab << tab << "if(newDataSet_" << i << "=\"1\") then" << endl
                 << tab << tab << tab << tab << "addressWriteCounter_" << i << " <= (others => '0');" << endl
                 << tab << tab << tab << "else" << endl
                 << tab << tab << tab << tab << "if(" << Layer::getValidDataName(i)+"_i" << "=\"1\") then" << endl
                 << tab << tab << tab << tab << tab << "addressWriteCounter_" << i << " <= std_logic_vector(unsigned(addressWriteCounter_" << i << ")+1);" << endl
                 << tab << tab << tab << tab << "end if;" << endl
                 << tab << tab << tab << "end if;" << endl;
        }
        vhdl << tab << tab << "end if;" << endl
             << tab << "end if;" << endl
             << "end process;" << endl;

        for(int i=0;i<numberOfFeatures;i++)
        {
            vhdl << "WriteAddress_" << i << " <= writeFlag_" << i << " & addressWriteCounter_" << i << ";" << endl;
            vhdl << "ReadAddress_" << i << " <= (not writeFlag_" << i << ") & addressReadCounter_" << i << ";" << endl;
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
