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
#include "Register.hpp"

using namespace std;
namespace flopoco {




    ConvolutionalCoreSimple::ConvolutionalCoreSimple(Target* target, unsigned int wordSize_, unsigned int fraction_, unsigned int weightWordSize_, unsigned int weightFraction_, unsigned int size_, vector <double> weights_, bool useAdderTree_, string roundingType_, string id_) :
        Operator(target), wordSize(wordSize_), fraction(fraction_), weightWordSize(weightWordSize_), weightFraction(weightFraction_), size(size_), weights(weights_), useAdderTree(useAdderTree_), roundingType(roundingType_), id(id_) {
        //throw error if fraction or weightFraction < 0
        if(fraction<0 || weightFraction<0)
        {
            stringstream e;
            e << "A fraction can't be smaller than 0";
            THROWERROR(e);
        }
        //throw error if the number of weights isn't correct
        unsigned int numberOfNeurons=size*size;
        if(weights.size()!=numberOfNeurons)
        {
            stringstream e;
            e << "The number of weights is not equal to the number of neurons in that layer";
            THROWERROR(e);
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

        if(useAdderTree==true)
        {
            cout << "pipelinedepth before adder tree: " << this->getPipelineDepth() << endl;
            IntAdderTree* treeOp = new IntAdderTree(target, (wordSize+weightWordSize), numberOfNeurons, "bitheap",false);
            cout << "pipelinedepth of the addder tree: " << treeOp->getPipelineDepth() << endl;
            this->addSubComponent(treeOp);
            for(unsigned int i=0; i<numberOfNeurons; i++)
            {
                this->inPortMap(treeOp, "X"+to_string(i+1), "prod"+to_string(i));
            }
            this->outPortMap(treeOp, "Y", "AddResult", true);
            this->vhdl << this->instance(treeOp, treeOp->getName()+"_instance");
            wordSizeTempOutput=this->getSignalByName("AddResult")->width();

            syncCycleFromSignal("AddResult"); // to set pipeline stages
        }
        else
        {
            this->vhdl << tab << declare("AddResult", wordSizeTempOutput) << " <= std_logic_vector(signed(prod0)";
            for(unsigned int i=1; i<numberOfNeurons; i++)
            {
                this->vhdl << " + signed(prod" << i << ")";
            }
            this->vhdl << ");" << endl;

        }
        this->outputWordSize=wordSizeTempOutput;
        this->outputFraction=fractionTempOutput;

        if(roundingType=="Truncation") //or other rounding types when they are implemented
        {
            roundOutput(wordSizeTempOutput, fractionTempOutput, roundingType);
        }
        else
        {
            // don't round
            this->outputWordSize=wordSizeTempOutput;
            this->outputFraction=fractionTempOutput;
            declare("R_t",this->outputWordSize);
            this->vhdl << tab << "R_t <= AddResult;" << endl;
        }

		addOutput("R",this->outputWordSize);
        this->vhdl << tab << "R <= R_t;" << endl;

    }

    void ConvolutionalCoreSimple::roundOutput(unsigned int wordSizeFrom, unsigned int fractionFrom, string round)
    {
        if(round=="Truncation")
        {
            int MSBFrom = wordSizeFrom - fractionFrom;
            int MSBTo = wordSize-fraction;
            int from = wordSizeFrom-(MSBFrom-MSBTo)-1;
            int downto = fractionFrom-fraction;
            this->vhdl << tab << declare("roundingNumber",downto) << " <= AddResult(" << downto-1 << ")";
            if(downto>0)
            {
                if(downto>1)
                {
                    this->vhdl << " & (" << downto-2 << " downto 0 => '0')";
                }
                this->vhdl <<";" << endl;
                this->vhdl << tab << declare("R_beforeRound",wordSizeFrom) << " <= std_logic_vector(signed(AddResult)+"
                           << "signed(roundingNumber));" << endl;

                declare("R_t",wordSize);
                this->outputWordSize=wordSize;
                this->outputFraction=fraction;
                this->vhdl << tab << "R_t <= R_beforeRound(" << from << " downto " << downto << ");" << endl;
            }
            else if(downto==0)
            {
                this->vhdl << tab << "R_t <= AddResult(" << from << " downto " << downto << ");" << endl;
            }
            else
            {
                stringstream e;
                e << "Error while rounding, requested fraction is bigger than actual fraction! This should normally not happen";
                THROWERROR(e);
            }

        }
        else
        {
            cout << "ConvolutionalCoreSimple.roundOutput: atm only truncation is implemented. FloPoCo is rounding with this mode instead" << endl;
            roundOutput(wordSizeFrom, fractionFrom, "Truncation");
        }
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
