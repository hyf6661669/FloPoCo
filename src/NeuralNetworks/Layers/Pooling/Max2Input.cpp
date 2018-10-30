// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "Max2Input.hpp"

#include "IntAddSubCmp/IntComparator.hpp"
#include "PrimitiveComponents/GenericMux.hpp"

using namespace std;
namespace flopoco {


    Max2Input::Max2Input(Target* target, unsigned int wordSize_, bool sign_, bool alsoReturnIndex, unsigned int indexWordSize) :
        Operator(target), wordSize(wordSize_), sign(sign_) {

        this->useNumericStd();

        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="Max2Input";

        // author
        setCopyrightString("Nicolai Fiege, 2017");

        // definition of the name of the operator
        ostringstream name;
        name << "Max2Input_wordSize_" << wordSize << "_sign_" << (sign==true?"true":"false") << "_alsoReturnIndex_" << alsoReturnIndex << "_wIndex_" << indexWordSize;
        setName(name.str());

        // add in/out
        addInput("X0", wordSize);
		addInput("X1", wordSize);
        addOutput("R", wordSize);

        if(alsoReturnIndex==true)
        {
            addInput("I0",indexWordSize);
            addInput("I1",indexWordSize);
            addOutput("IMax",indexWordSize);
        }
        nextCycle(); // input registers

        // vhdl code
        if(sign==true)
        {
            this->vhdl << "R <= X0 when (signed(X0) > signed(X1)) else X1;" << endl;
            if(alsoReturnIndex==true)
            {
                this->vhdl << "IMax <= I0 when (signed(X0) > signed(X1)) else I1;" << endl;
            }
        }
        else
        {
            this->vhdl << "R <= X0 when (unsigned(X0) > unsigned(X1)) else X1;" << endl;
            if(alsoReturnIndex==true)
            {
                this->vhdl << "IMax <= I0 when (unsigned(X0) > unsigned(X1)) else I1;" << endl;
            }
        }

    }

}//namespace flopoco
