//
// Created by nfiege on 18/05/18.
//

#include "ConvolutionalCoreWithWeightInputs.hpp"
#include "NeuralNetworks/Utility/BitheapWrapper.hpp"
#include "NeuralNetworks/Utility/Rounding.hpp"

namespace flopoco
{
    ConvolutionalCoreWithWeightInputs::ConvolutionalCoreWithWeightInputs(Target* target, multiplicationType multType, unsigned int wIn, unsigned int wF, unsigned int wWeightsIn, unsigned int wFWeights, unsigned int size, char roundingType, string id)
    : Operator(target), wordSize(wIn), fraction(wF), weightWordSize(wWeightsIn), weightFraction(wFWeights), outputWordSize(0), outputFraction(0)
    {
        // author
        setCopyrightString("Nicolai Fiege, 2017");

        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="ConvolutionalCoreWithWeightInputs";

        // definition of the name of the operator
        ostringstream name;
        name << "ConvolutionalCoreWithWeightInputs_" << id;
        setName(name.str());
        this->useNumericStd();

        for(unsigned int i=0; i<size*size; i++)
        {
            addInput("X"+to_string(i),wordSize);
            addInput("W"+to_string(i),weightWordSize);
            declare("W"+to_string(i)+"_reg",weightWordSize);
        }
        addInput("Weight_Enable",1);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////
        // set new weights, when DMA controller has new weights and the current feature finished calculating //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////
        this->vhdl << "process(clk)" << endl;
        this->vhdl << "begin" << endl;
        this->vhdl << tab << "if(rising_edge(clk)) then" << endl;
        this->vhdl << tab << tab << "if(rst='1') then" << endl;
        for(unsigned int i=0; i<size*size; i++)
        {
            this->vhdl << tab << tab << tab << "W" << i << "_reg <= (others => '0');" << endl;
        }
        this->vhdl << tab << tab << "elsif(Weight_Enable = \"1\") then" << endl;
        for(unsigned int i=0; i<size*size; i++)
        {
            this->vhdl << tab << tab << tab << "W" << i << "_reg <= W" << i << ";" << endl;
        }
        this->vhdl << tab << tab << "end if;" << endl;
        this->vhdl << tab << "end if;" << endl;
        this->vhdl << "end process;" << endl;

        ////////////////////////////////////////////
        // Calculate P = X * W_reg for each input //
        ////////////////////////////////////////////
        int productWidth = wordSize+weightWordSize;
        for(unsigned int i=0; i<size*size; i++)
        {
            NeuralNetworkMultiplication* mult = new NeuralNetworkMultiplication(target,wordSize,weightWordSize,multType);
            addSubComponent(mult);
            inPortMap(mult,"X","X"+to_string(i));
            inPortMap(mult,"Y","W"+to_string(i)+"_reg");
            outPortMap(mult,"R","P"+to_string(i),true);
            vhdl << instance(mult,"Multiplication_instance"+to_string(i));

            if(getSignalByName("P"+to_string(i))->width() != productWidth) THROWERROR("Bug in NeuralNetworkMultiplication! Wrong output word size!");
        }
        syncCycleFromSignal("P0"); // all multiplicators should have the same pipeline depth anyways, so just choose one of them


        //////////////////////////////
        // Calculate Sum of all P's //
        //////////////////////////////
        vector<int> wIns(size*size,productWidth);
        vector<bool> signs(size*size,true);
        vector<int> weights(size*size,0);
        BitheapWrapper* bit = new BitheapWrapper(target,wIns,size*size,signs,weights);
        addSubComponent(bit);
        for(unsigned int i=0; i<size*size; i++)
        {
            inPortMap(bit,"X"+to_string(i),"P"+to_string(i));
        }
        outPortMap(bit,"R","AddResult",true);
        this->vhdl << instance(bit,"Sum_of_all_Products");
        syncCycleFromSignal("AddResult");

        ///////////
        // Round //
        ///////////
        this->outputWordSize=this->getSignalByName("AddResult")->width();
        this->outputFraction=this->weightFraction+this->fraction;
        if(roundingType==0x00)
        {
            this->roundOutput(target,this->outputWordSize, this->weightFraction+this->fraction,"AddResult",roundingType);
            this->addOutput("R",this->outputWordSize);
            nextCycle();
            this->vhdl << "R <= R_t;" << endl;
        }
        else
        {
            this->addOutput("R",this->outputWordSize);
            nextCycle();
            this->vhdl << "R <= AddResult;" << endl;
        }

    }

    void ConvolutionalCoreWithWeightInputs::roundOutput(Target* target, unsigned int wordSizeFrom, unsigned int fractionFrom, string roundFromName, char roundingType)
    {
        flopoco::roundingTypeEnum rType;
        if(roundingType==0x00)
        {
            rType = flopoco::roundingTypeEnum::truncation;
        }
        else if(roundingType==0x01)
        {
            rType = flopoco::roundingTypeEnum::saturation;
        }
        else
        {
            cout << "ConvolutionalCoreSimple.roundOutput: atm only truncation/saturation are implemented. FloPoCo is rounding with truncation instead" << endl;
            rType = flopoco::roundingTypeEnum::truncation;
        }
        Rounding* roundOp = new Rounding(target,wordSizeFrom,fractionFrom,this->wordSize,this->fraction,rType);
        addSubComponent(roundOp);
        inPortMap(roundOp,"X","AddResult");
        outPortMap(roundOp,"R","R_t",true);
        this->vhdl << instance(roundOp,"Rounding_instance");

        this->outputWordSize=this->wordSize;
        this->outputFraction=this->fraction;
    }
}