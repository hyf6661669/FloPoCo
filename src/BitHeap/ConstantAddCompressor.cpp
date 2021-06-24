//
// Created by Annika Oeste on 06.05.21.
//

#include "ConstantAddCompressor.hpp"

using namespace std;
namespace flopoco {

    ConstantAddCompressor::ConstantAddCompressor(Operator *parentOp, Target *target, vector<int> _heights, int constant) : Compressor(
            parentOp, target) {
        setCopyrightString("Annika Oeste");
        useNumericStd();


        //compressors are supposed to be combinatorial
        setCombinatorial();
        setShared();

        ostringstream name;
        name << "ConstantAddCompressor_Const_" << abs(constant);
        setNameWithFreqAndUID(name.str());

        heights = _heights;
        vector<int> _outHeights;
        for(int i = 0; i < heights.size(); i++){
            _outHeights.push_back(1);
        }
        outHeights = _outHeights;
        //wOut = heights.size()+1;

        //setWordSizes();
        wIn = heights.size();
        wOut = outHeights.size();

        createInputsAndOutputs();

        // create signal X as input concatenation
        bool isFirst = true;
        vhdl << tab << declare(
                "X", wIn, true) << tab << "<=";
        for(int i=heights.size()-1; i>=0; i--)
        {
            //no need to create a signal for columns of height 0
            if(heights[i] > 0)
            {
                if (isFirst) {
                    isFirst = false;
                    vhdl << tab << join("X",i);
                } else {
                    vhdl << tab << "&" << tab << join("X",i);
                }

            }
        }
        vhdl << ";" << endl;

        vhdl << tab << "R <= STD_LOGIC_VECTOR(SIGNED(X) + (" << constant << "));" << endl;
    }

    ConstantAddCompressor::~ConstantAddCompressor() {
    }

    OperatorPtr ConstantAddCompressor::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args) {
        int wIn;
        int constant;

        UserInterface::parseInt(args, "wIn", &wIn);
        UserInterface::parseInt(args, "constant", &constant);

        vector<int> heights;
        for(int i = 0; i < wIn; i++){
            heights.push_back(1);
        }

        return new ConstantAddCompressor(parentOp, target, heights, constant);
    }

    void ConstantAddCompressor::registerFactory(){
        UserInterface::add("ConstantAddCompressor", // name
                           "A compressor for addition with a constant.",
                           "BasicInteger", // categories
                           "",
                           "wIn(int): input length, \
				            constant(int): constant to add, has to be within input length",
                           "",
                           Compressor::parseArguments
        ) ;
    }

    BasicConstantAddCompressor::BasicConstantAddCompressor(Operator* parentOp_, Target * target, vector<int> _heights, int constant) : BasicCompressor(parentOp_, target, _heights, _heights.size(), CompressorType::Constant, true) {
        heights = _heights;
        _constant = constant;
    }

    Compressor* BasicConstantAddCompressor::getCompressor(unsigned int middleLength) {
        return new ConstantAddCompressor(parentOp, target, heights, _constant);
    }
}