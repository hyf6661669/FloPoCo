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
#include "ReLU.hpp"

using namespace std;
namespace flopoco {




    ReLU::ReLU(Target* target, unsigned int wordSize_) :
        Operator(target), wordSize(wordSize_) {

		// definition of the source file name, used for info and error reporting using REPORT 
		srcFileName="ReLU";

        // author
        setCopyrightString("Nicolai Fiege, 2017");

		// definition of the name of the operator
		ostringstream name;
        name << "ReLU_wordSize_" << wordSize;
        setName(name.str());

        //add in/out
        addInput("X", wordSize);
        addOutput("R", wordSize);

        vhdl << declare("X_top", wordSize) << " <= (others => X(" << wordSize-1 << "));" << endl;
        vhdl << "R <= X and (not X_top);" << endl;
    }

}//namespace flopoco
