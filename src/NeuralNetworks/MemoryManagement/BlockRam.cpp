// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "BlockRam.hpp"

using namespace std;
namespace flopoco {
	
    BlockRam::BlockRam(Target* target, unsigned int dataWidth_, unsigned int addressWidth_) :
        Operator(target), dataWidth(dataWidth_), addressWidth(addressWidth_) {

        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="BlockRam";

        // use numeric_std
        useNumericStd();

        // author
        setCopyrightString("Nicolai Fiege, 2017");

        // definition of the name of the operator
        ostringstream name;
        name << "BlockRam_dataWidth_" << dataWidth << "_addressWidth_" << addressWidth;
        setName(name.str());
		
        //compute depth of the RAM
        depth = mpz_class(1) << addressWidth;

        addInput("Address_W_in", addressWidth);
        addInput("Address_R_in", addressWidth);
        addInput("Data_in", dataWidth);
        addInput("WriteEnable",1);
        addOutput("Data_out", dataWidth);

        //type declaration for the array
        stringstream typeStream;
        typeStream << "array(0 to " << depth-mpz_class(1) << ") of std_logic_vector(" << dataWidth-1 << " downto 0);" << endl
                           << "signal ram : ram_t";
        addType("ram_t", typeStream.str());


         vhdl
         << tab << "process(clk)" << endl
         << tab << "begin" << endl
         << tab << tab << "if (rising_edge(clk)) then" << endl
         << tab << tab << tab << "if (WriteEnable=\"1\") then" << endl
         << tab << tab << tab << tab << "ram(to_integer(unsigned(Address_W_in))) <= Data_in;" << endl
         << tab << tab << tab << "end if;" << endl
         << tab << tab << tab << "Data_out <= ram(to_integer(unsigned(Address_R_in)));" << endl
         << tab << tab << "end if;" << endl
         << tab << "end process;" << endl;
		
    }

}//namespace flopoco
