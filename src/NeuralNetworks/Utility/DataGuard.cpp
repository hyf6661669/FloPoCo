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
#include "DataGuard.hpp"

using namespace std;
namespace flopoco {




    DataGuard::DataGuard(Target* target, int wordSize_, int wordSizeGuard_, unsigned int numberToCompare_) :
        Operator(target), wordSize(wordSize_), wordSizeGuard(wordSizeGuard_), numberToCompare(numberToCompare_) {

        useNumericStd();

        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="DataGuard";

        // author
        setCopyrightString("Nicolai Fiege, 2018");

        // definition of the name of the operator
        ostringstream name;
        name << "DataGuard_" << this->wordSize << "_" << this->wordSizeGuard << "_" << this->numberToCompare;
        setName(name.str());
		

        // add in/out
        addInput("Data_in", this->wordSize);
        addInput("Guard_in", this->wordSizeGuard);
        addOutput("Data_out", this->wordSize);
		
		// logic
		this->vhdl << "with to_integer(unsigned(Guard_in)) select Data_out <= " << endl
                   << tab << "Data_in when " << this->numberToCompare << "," << endl
                   << tab << "(others => '0') when others;" << endl;
        
		
    }

}//namespace flopoco
