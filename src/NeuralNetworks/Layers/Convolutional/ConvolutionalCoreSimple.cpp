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
#include "ConvolutionalCoreSimple.hpp"

#include "IntAddSubCmp/IntAdderTree.hpp"
#include "BitHeap/BitHeap.hpp"
#include "NeuralNetworks/Utility/Register.hpp"
#include "NeuralNetworks/Utility/Rounding.hpp"

using namespace std;
namespace flopoco {




    ConvolutionalCoreSimple::ConvolutionalCoreSimple(Target* target, unsigned int wordSize_, unsigned int fraction_, unsigned int weightWordSize_, unsigned int weightFraction_, unsigned int size_, vector <double> weights_, bool useBitHeap_, char roundingType_, string id_) :
        Operator(target), wordSize(wordSize_), fraction(fraction_), weightWordSize(weightWordSize_), weightFraction(weightFraction_), size(size_), weights(weights_), useBitHeap(useBitHeap_), roundingType(roundingType_), id(id_) {
        //throw error if fraction or weightFraction < 0
        if(fraction<0 || weightFraction<0)
        {
            stringstream e;
            e << "A fraction can't be smaller than 0";
            THROWERROR(e.str());
        }
        //throw error if the number of weights isn't correct
        unsigned int numberOfNeurons=size*size;
        if(weights.size()!=numberOfNeurons)
        {
            stringstream e;
            e << "The number of weights is not equal to the number of neurons in that layer";
            THROWERROR(e.str());
        }

        // author
        setCopyrightString("Nicolai Fiege, 2017");

		// definition of the source file name, used for info and error reporting using REPORT 
		srcFileName="ConvolutionalCoreSimple";

		// definition of the name of the operator
		ostringstream name;
        name << "ConvolutionalCoreSimple_" << id;
        setName(name.str());
        this->useNumericStd();

        //add inputs and output
        for(unsigned int i=0; i<numberOfNeurons; i++)
        {
            inputNames.push_back("X"+to_string(i));
            addInput(inputNames[i], wordSize);
        }

        //calculate products
        for(unsigned int i=0; i<numberOfNeurons; i++)
        {
            double scal = pow(2,weightFraction);
            int val = (int) round(weights[i] * scal);

            addConstant("const"+to_string(i),"signed", "to_signed("+to_string(val)+","+to_string(weightWordSize)+")");
            vhdl << tab << declare("prod"+to_string(i),wordSize+weightWordSize) << " <= std_logic_vector(signed("
                 << inputNames[i] << ")*const" << i << ");" << endl;
        }
        this->nextCycle();

        //sum everything up and assign to temporary output (bigger wordsize than actual output)
        unsigned int wordSizeTempOutput = wordSize+weightWordSize;
        unsigned int fractionTempOutput = fraction+weightFraction;

        if(useBitHeap==true)
        {
            // bit heap
            int maxWeightBitheap = ceil(log2(numberOfNeurons)) + wordSizeTempOutput + 2;
            BitHeap* bitH = new BitHeap(this, maxWeightBitheap, true, "AddProducts");
            for(unsigned int i=0; i<numberOfNeurons; i++)
            {
                bitH->addSignedBitVector(0,"prod"+to_string(i),wordSizeTempOutput);
            }
            bitH->useVariableColumnCompressors = false;
            bitH->generateCompressorVHDL(false);
            nextCycle();

            // set wordSizeTempOutput
            wordSizeTempOutput=bitH->getMaxWeight();

            // don't use MSB of bitheap because it's corrupted!
            vhdl << tab << declare("AddResult", wordSizeTempOutput) << " <= " << bitH->getSumName() << "(" << wordSizeTempOutput-1 << " downto 0)" << ";" << endl;


        }
        else
        {
            // calc wordsize of sum
            int numberOfStages = ceil(log2(numberOfNeurons+1));
            int wSum = wordSizeTempOutput + numberOfStages;
            // sign extend products to have needed word size
            for(int i=0;i<numberOfNeurons;i++)
            {
                string tmp="prod"+to_string(i)+"_resized";
                this->vhdl << tab << declare(tmp,wSum) << " <= std_logic_vector(resize(signed(prod" << i << ")," << wSum << "));" << endl;
            }

            if(numberOfNeurons==9)
            {
                // first stage
                this->vhdl << tab << declare("sum01",wSum) << " <= std_logic_vector(signed(prod0_resized)+signed(prod1_resized));" << endl;
                this->vhdl << tab << declare("sum23",wSum) << " <= std_logic_vector(signed(prod2_resized)+signed(prod3_resized));" << endl;
                this->vhdl << tab << declare("sum45",wSum) << " <= std_logic_vector(signed(prod4_resized)+signed(prod5_resized));" << endl;
                this->vhdl << tab << declare("sum67",wSum) << " <= std_logic_vector(signed(prod6_resized)+signed(prod7_resized));" << endl;
                nextCycle();
                // second stage
                this->vhdl << tab << declare("sum0123",wSum) << " <= std_logic_vector(signed(sum01)+signed(sum23));" << endl;
                this->vhdl << tab << declare("sum4567",wSum) << " <= std_logic_vector(signed(sum45)+signed(sum67));" << endl;
                nextCycle();
                // third stage
                this->vhdl << tab << declare("sum01234567",wSum) << " <= std_logic_vector(signed(sum0123)+signed(sum4567));" << endl;
                nextCycle();
                // fourth stage
                this->vhdl << tab << declare("AddResult",wSum) << " <= std_logic_vector(signed(sum01234567)+signed(prod8_resized));" << endl;
                nextCycle();

            }
            else
            {
                this->vhdl << tab << declare("AddResult", wSum) << " <= std_logic_vector(signed(prod0_resized)";
                for(unsigned int i=1; i<numberOfNeurons; i++)
                {
                    this->vhdl << " + signed(prod" << i << "_resized)";
                }
                this->vhdl << ");" << endl;
                nextCycle();
            }

        }

        if(roundingType==0x00 || roundingType ==0x01) //or other rounding types when they are implemented
        {
            this->outputWordSize=this->wordSize;
            this->outputFraction=this->fraction;
            roundOutput(target, wordSizeTempOutput, fractionTempOutput, roundingType);
            addOutput("R",this->outputWordSize);
            nextCycle();
            this->vhdl << tab << "R <= R_t;" << endl;
        }
        else
        {
            // don't round
            this->outputWordSize=wordSizeTempOutput;
            this->outputFraction=fractionTempOutput;
            addOutput("R",this->outputWordSize);
            nextCycle();
            this->vhdl << tab << "R <= AddResult;" << endl;
        }


    }

    void ConvolutionalCoreSimple::roundOutput(Target* target,unsigned int wordSizeFrom, unsigned int fractionFrom, char round)
    {
        flopoco::roundingTypeEnum rType;
        if(round==0x00)
        {
            rType = flopoco::roundingTypeEnum::truncation;
        }
        else if(round==0x01)
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
	
	
    unsigned int ConvolutionalCoreSimple::getOutputWordSize()
	{
		return this->outputWordSize;
	}
	
    unsigned int ConvolutionalCoreSimple::getOutputFraction()
	{
		return this->outputFraction;
	}

}//namespace flopoco
