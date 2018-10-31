// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "FullConnectedCore.hpp"

#include "NeuralNetworks/Utility/BitheapWrapper.hpp"
#include "NeuralNetworks/Utility/Rounding.hpp"

using namespace std;
namespace flopoco {

    FullConnectedCore::FullConnectedCore(Target* target, multiplicationType multType, int wordSize, int fraction, int weightWordSize, int weightFraction, int secondInputWidth, char roundingMode)
            : Operator(target)
    {
        useNumericStd();
        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="FullConnectedCore";

        // author
        setCopyrightString("Nicolai Fiege, 2018");

        // definition of the name of the operator
        ostringstream name;
        name << "FullConnectedCore_" << wordSize << "_" << fraction << "_" << weightWordSize << "_" << weightFraction << "_" << secondInputWidth << "_" << ((int)roundingMode) << "_" << multType;
        setName(name.str());

        this->outputWordSize=-1;

        addInput("X",wordSize);
        addInput("Y",secondInputWidth);
        addInput("W",weightWordSize);

        // sync inputs
        this->vhdl << declare("X_reg",wordSize) << " <= X;" << endl;
        this->vhdl << declare("Y_reg",secondInputWidth) << " <= Y;" << endl;
        this->vhdl << declare("W_reg",weightWordSize) << " <= W;" << endl;
        nextCycle();

        int productSize = wordSize + weightWordSize;
        NeuralNetworkMultiplication* mult = new NeuralNetworkMultiplication(target,wordSize,weightWordSize,multType);
        addSubComponent(mult);
        inPortMap(mult,"X","X_reg");
        inPortMap(mult,"Y","W_reg");
        outPortMap(mult,"R","Product",true);
        vhdl << instance(mult,"Multiplication_instance");
        syncCycleFromSignal("Product");

        vector<int> wIns;
        wIns.push_back(productSize);
        wIns.push_back(secondInputWidth);
        vector<bool> signs(2,true);
        vector<int> weights(2,0);
        BitheapWrapper* bitHeap = new BitheapWrapper(target,wIns,wIns.size(),signs,weights);
        addSubComponent(bitHeap);
        inPortMap(bitHeap,"X0","Product");
        inPortMap(bitHeap,"X1","Y_reg");
        outPortMap(bitHeap,"R","AddResult",true);
        this->vhdl << instance(bitHeap,"BitHeap_instance");
        syncCycleFromSignal("AddResult");

        roundingTypeEnum round = roundingTypeEnum::saturation;
        if(roundingMode==0x00)
        {
            round = roundingTypeEnum::truncation;
        }

        Rounding* roundOp = new Rounding(target,getSignalByName("AddResult")->width(),fraction+weightFraction,secondInputWidth,fraction+weightFraction,round);
        addSubComponent(roundOp);
        inPortMap(roundOp,"X","AddResult");
        outPortMap(roundOp,"R","R_temp",true);
        this->vhdl << instance(roundOp,"Rounding_instance");

        addOutput("R",getSignalByName("R_temp")->width());
        this->vhdl << "R <= R_temp;" << endl;
    }
}