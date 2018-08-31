// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "WindowShiftRegister.hpp"

using namespace std;
namespace flopoco {




    WindowShiftRegister::WindowShiftRegister(Target* target, unsigned int wordSize_, unsigned int windowSize_, unsigned int rowLength_) :
        Operator(target), wordSize(wordSize_), windowSize(windowSize_), rowLength(rowLength_) {


        //author
        setCopyrightString("Nicolai Fiege, 2017");

        //definition of the source file name, used for info and error reporting using REPORT
        srcFileName="WindowShiftRegister";

        //definition of the operator name
        ostringstream name;
        name << "WindowShiftRegister_wordSize" << wordSize << "_windowSize_" << windowSize << "_rowLength_" << rowLength_;
        setName(name.str());

        //add inputs and outputs
        addInput("X",wordSize);
        addInput("enable",1);
        addInput("newStep",1);
        unsigned int numberOfOutputs=windowSize*windowSize;
        for(unsigned int i=0; i<numberOfOutputs; i++)
        {
            string outNameTemp = "R" + to_string(i);
            outputNames.push_back(outNameTemp);
            addOutput(outNameTemp,wordSize);
        }

        //calculate buffer length and declare the array
        numberOfFlipFlops = (windowSize-1)*(rowLength+1);
        stringstream typeStream;
        typeStream << "array(0 to " << numberOfFlipFlops << ") of std_logic_vector(" << wordSize-1 << " downto 0);" << endl
                   << "signal reg : reg_t";
        addType("reg_t", typeStream.str());

        //process for the array (assign all signals)
        this->vhdl << tab << "process(clk)" << endl
                   << tab << "begin" << endl
                   << tab << tab << "if(rising_edge(clk)) then" << endl
                   << tab << tab << tab << "if(rst='1') then" << endl
                   << tab << tab << tab << tab << "reg <= (others => (others => '0'));" << endl
                   << tab << tab << tab << "else" << endl
                   << tab << tab << tab << tab << "if(newStep=\"1\") then" << endl
                   << tab << tab << tab << tab << tab << "reg <= (others => (others => '0'));" << endl
                   << tab << tab << tab << tab << "else" << endl
                   << tab << tab << tab << tab << tab << "if(enable=\"1\") then" << endl
                   << tab << tab << tab << tab << tab << tab << "reg(0) <= X;" << endl;
        for(unsigned int i=0; i<numberOfFlipFlops; i++)
        {
            this->vhdl << tab << tab << tab << tab << tab << tab << "reg(" << i+1 << ") <= reg(" << i << ");" << endl;
        }
        this->vhdl << tab << tab << tab <<tab << tab << "end if;" << endl
                   << tab << tab << tab << tab << "end if;" << endl
                   << tab << tab << tab << "end if;" << endl
                   << tab << tab << "end if;" << endl
                   << tab << "end process;" << endl;

        //assign outputs
        for(unsigned int y=0; y<windowSize; y++)
        {
            for(unsigned int x=0; x<windowSize; x++)
            {
                this->vhdl << tab << "R" << y*windowSize+x << " <= reg(" << y*rowLength+x << ");" << endl;
            }
        }
    }

}//namespace flopoco
