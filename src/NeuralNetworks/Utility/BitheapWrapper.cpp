//
// Created by nfiege on 23/05/18.
//

// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "BitheapWrapper.hpp"
#include "BitHeap/BitHeap.hpp"

using namespace std;
namespace flopoco {




    BitheapWrapper::BitheapWrapper(Target* target, vector<int> wIn, int numberOfInputs, vector<bool> signs, vector<int> weights) :
            Operator(target) {

        useNumericStd();

        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName = "BitheapWrapper";

        // author
        setCopyrightString("Nicolai Fiege, 2018");

        // definition of the name of the operator
        ostringstream name;
        name << "BitheapWrapper_" << numberOfInputs;
        for (auto it : wIn) {
            name << "_" << it;
        }
        for (auto it : signs) {
            if (it == true) {
                name << "_s";
            } else {
                name << "_u";
            }
        }
        for (auto it : weights) {
            name << "_" << it;
        }
        setName(name.str());

        if(numberOfInputs != (int)signs.size())
        {
            THROWERROR("Size of sign vector doesn't match number of inputs");
        }

        if(numberOfInputs != (int)wIn.size())
        {
            THROWERROR("Size of word size vector doesn't match number of inputs");
        }

        if(numberOfInputs != (int)weights.size())
        {
            THROWERROR("Size of weight vector doesn't match number of inputs");
        }

        int maxWIn = 0;
        for(auto it : wIn)
        {
            if(it<=0) THROWERROR("Word Size must be > 0!");
            if(it>maxWIn) maxWIn = it;
        }

        for(auto it : weights)
        {
            if(it<0) THROWERROR("Weights must be > 0!");
        }

        int maxWeight = ceil(log2(numberOfInputs)) + maxWIn + 2;
        BitHeap* bh = new BitHeap(this, maxWeight, true, "Bit_Heap_Out");
        for(int i = 0; i < numberOfInputs; i++)
        {
            addInput("X"+to_string(i),wIn[i]);
            int weight = weights[i];
            if(signs[i]==true)
            {
                bh->addSignedBitVector(weight,"X"+to_string(i),wIn[i]);
            }
            else
            {
                bh->addUnsignedBitVector(weight,"X"+to_string(i),wIn[i]);
            }
        }
        //nextCycle(); // input register
        bh->useVariableColumnCompressors = false;
        bh->generateCompressorVHDL(false);
        string sumName = bh->getSumName();
        this->wOut = this->getSignalByName(sumName)->width()-1;

        addOutput("R",this->wOut);

        this->syncCycleFromSignal(sumName);
        nextCycle(); // output register

        // MSB of Bitheap is garbage => ignore it!
        this->vhdl << "R <= " << sumName << "(" << this->wOut-1 << " downto 0);" << endl;
    }

}//namespace flopoco