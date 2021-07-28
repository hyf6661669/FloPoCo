//
// Created by annika on 27.07.21.
//

#include "ConstantToBitheap.hpp"

using namespace std;
namespace flopoco {

    ConstantToBitheap::ConstantToBitheap(Operator *parentOp, Target *target, int constant) : Compressor(
            parentOp, target) {
        setCopyrightString("Annika Oeste");
        useNumericStd();


        //compressors are supposed to be combinatorial
        setCombinatorial();
        setShared();

        ostringstream name;
        name << "ConstantToBitheap_Const_" << abs(constant);
        setNameWithFreqAndUID(name.str());

        //heights = vector<int> newHeights();
        int sizeOut =  intlog2(constant);
        vector<int> _outHeights;
        for(int i = 0; i < sizeOut; i++){
            if (((1 << i) & constant) != 0) {
                _outHeights.push_back(1);
            } else {
                _outHeights.push_back(0);
            }
        }
        outHeights = _outHeights;
        //wOut = heights.size()+1;

        //setWordSizes();
        wIn = 0;
        wOut = outHeights.size();

        createInputsAndOutputs();

        vhdl << tab << "R <= STD_LOGIC_VECTOR(TO_UNSIGNED(" << constant << "," << sizeOut << "));" << endl;
    }

    ConstantToBitheap::~ConstantToBitheap() {
    }

    OperatorPtr ConstantToBitheap::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args) {
        //int wIn;
        int constant;

        //UserInterface::parseInt(args, "wIn", &wIn);
        UserInterface::parseInt(args, "constant", &constant);

//        vector<int> heights;
//        for(int i = 0; i < wIn; i++){
//            heights.push_back(1);
//        }

        return new ConstantToBitheap(parentOp, target, constant);
    }

    void ConstantToBitheap::registerFactory(){
        UserInterface::add("ConstantToBitheap", // name
                           "A compressor for adding a constant to the bit heap without inputs.",
                           "BasicInteger", // categories
                           "",
                           "wIn(int): input length, \
				            constant(int): constant to add, has to be within input length",
                           "",
                           Compressor::parseArguments
        ) ;
    }

    BasicConstantToBitheap::BasicConstantToBitheap(Operator* parentOp_, Target * target, int constant) : BasicCompressor(parentOp_, target, vector<int>(), 0, CompressorType::Constant, true) {
        //heights = _heights;
        _constant = constant;
    }

    Compressor* BasicConstantToBitheap::getCompressor(unsigned int middleLength) {
        return new ConstantToBitheap(parentOp, target, _constant);
    }
}