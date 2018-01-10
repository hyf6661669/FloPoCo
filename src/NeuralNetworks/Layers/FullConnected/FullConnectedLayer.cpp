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
#include "FullConnectedLayer.hpp"

using namespace std;
namespace flopoco {




    FullConnectedLayer::FullConnectedLayer(Target* target, unsigned int numberOfInputs_, unsigned int numberOfNeurons_, string id_) :
        Operator(target), numberOfInputs(numberOfInputs_), numberOfNeurons(numberOfNeurons_), id(id_) {

		// definition of the source file name, used for info and error reporting using REPORT 
		srcFileName="FullConnectedLayer";

        // author
        setCopyrightString("Nicolai Fiege, 2017");

		// definition of the name of the operator
		ostringstream name;
        name << "FullConnectedLayer_" << id;
        setName(name.str());
		
        addInput("a",8);
        addOutput("b",8);

        vhdl << "b <= a;" << endl;
    }

}//namespace flopoco
