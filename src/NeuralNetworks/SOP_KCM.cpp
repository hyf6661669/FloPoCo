// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "SOP_KCM.hpp"
#include "../PrimitiveComponents/Xilinx/Xilinx_CFGLUT5.hpp"
#include "../BitHeap/BitHeap.hpp"
#include <math.h>

using namespace std;
namespace flopoco {




    SOP_KCM::SOP_KCM(Target* target,int inputWordSize,int constantWordSize, int parameterNo, int _defaultConstant) : Operator(target), input_bit_width(inputWordSize), Const_bit_width(constantWordSize), No_of_Products(parameterNo)
    {
        srcFileName="SOP_KCM";

		ostringstream name;
        name << "SOP_KCM_" + to_string(input_bit_width)+ "_" + to_string(Const_bit_width)+ "_" + to_string(No_of_Products);
		setName(name.str());
		// Copyright 
        setCopyrightString("UNIVERSITY of Kassel 2017");

        // MH debug int input_bit_width = 12;
        // MH debug int Const_bit_width = 4;
        // MH debug int No_of_Products = 2;



        signed_calculation = true;
        faithful_rounding = true;
        int defaultConstant = _defaultConstant;

        LUT_bit_width = 4;

        // abgeleitete variablen... MH Debug
        int LUT_per_stage = Const_bit_width+LUT_bit_width;
        //int No_of_stages = ceil((float)Const_bit_width / (float)LUT_bit_width);// Version mit der die ergebnisse produziert wurrden
	int No_of_stages = ceil((float)input_bit_width / (float)LUT_bit_width);
	

        //to prevent a overvlow, caused by mutiple additions.
        int additionalBitWidth = floor(log2((float)No_of_Products))+1;

        // gard bits for faithful rounding
        g =  floor(log2((float)(No_of_stages * No_of_Products)));

        output_bit_width = inputWordSize;



        addInput ("LUT_Config_clk_enable");

        vector<string> bitHeapStack;
        vector<unsigned int> bitHeapStackShifts;
        vector<int> bitHeapStackSize;
        bitHeapStackShifts.clear();
        bitHeapStack.clear();
        bitHeapStackSize.clear();

        int outputBits = Const_bit_width + input_bit_width + additionalBitWidth;
        int border =  outputBits - output_bit_width - g;



        addFullComment("Parameter:");
        addComment(join("input_bit_width=",input_bit_width));
        addComment(join("Const_bit_width=",Const_bit_width));
        addComment(join("output_bit_width=",output_bit_width));
        addComment(join("outputBits=",outputBits));
        addComment(join("border=",border));
        addComment(join("No_of_Products=",No_of_Products));
        addComment(join("No_of_stages=",No_of_stages));
        addComment(join("g=",g));

        nextCycle();

        for (int prNo = 0; prNo < No_of_Products; ++prNo)
        {
            string inputSignalName = join("X",prNo);
            addInput ( inputSignalName, input_bit_width);
            inputSignalName+="(";

            addFullComment(join("Product No ",prNo));
            string configurationStreamSignalName = "cdi_no_" + to_string(prNo);
            addInput(configurationStreamSignalName,ceil(((float)LUT_per_stage)/2)*No_of_stages);

            configurationStreamSignalName += "(";

            for(int stage = 0; stage < No_of_stages; ++stage)
            {
                addFullComment(join("Stage No ",stage));
                int Lut_start_Counter;
                if(faithful_rounding)
                {
                    Lut_start_Counter = max(0, border - stage * LUT_bit_width);
                }
                else
                {
                    Lut_start_Counter = 0;
                }

                if(Lut_start_Counter < LUT_per_stage)
                {
                    string outputSignalName = "ProductNo_" + to_string(prNo) + "_LUTstage_" + to_string(stage);
                    declare(outputSignalName, LUT_per_stage-Lut_start_Counter);
                    bitHeapStack.push_back(outputSignalName);
                    bitHeapStackShifts.push_back(stage*LUT_bit_width+Lut_start_Counter);
                    bitHeapStackSize.push_back(LUT_per_stage-Lut_start_Counter);

                    for(int LUT_No = Lut_start_Counter; LUT_No< LUT_per_stage; LUT_No += 2)
                    {

                        vhdl << std::endl;
                        Xilinx_CFGLUT5 *myCFGLUT = new Xilinx_CFGLUT5(target);
                        addToGlobalOpList(myCFGLUT);

                        myCFGLUT->setGeneric( "init", generateInitStringFor(defaultConstant,LUT_No) );

                        inPortMapCst(myCFGLUT, "CLK","clk");
                        inPortMap(myCFGLUT, "CE","LUT_Config_clk_enable");
                        inPortMap(myCFGLUT, "CDI",join(configurationStreamSignalName,to_string((int)(floor(LUT_No/2) + ceil(((float)LUT_per_stage)/2)*stage)),")"));

                        for(int i = 0; i <= LUT_bit_width; ++i)
                        {
                            if(i == LUT_bit_width)
			    {
                             inPortMapCst(myCFGLUT, join("I", to_string(i)),"'1'");// the highest bit is true to use the 5 input LUT as two 4 input Luts
			    }
			    else if((i+stage*LUT_bit_width) < input_bit_width)
			    {
			     inPortMap(myCFGLUT, join("I", to_string(i)),join(inputSignalName,to_string(i+stage*LUT_bit_width),")"));
			    }
                            else
			    {
                             inPortMapCst(myCFGLUT, join("I", to_string(i)),"'0'");
			    }
                        }

                        string outputSignalNameLUT = outputSignalName + "(" + to_string(LUT_No-Lut_start_Counter) +")";
                        outPortMap(myCFGLUT, "O6",outputSignalNameLUT,false);// MH switch 5 6
			
			
                        if(LUT_No+1 < LUT_per_stage)// if the nomber of LUts are odd
                        {
                            outputSignalNameLUT = outputSignalName + "(" + to_string(LUT_No+1-Lut_start_Counter) +")";
                            outPortMap(myCFGLUT, "O5",outputSignalNameLUT,false);// MH switch 5 6
                        }
                        else
                        {
                            outPortMap(myCFGLUT, "O5","open",false);// MH switch 5 6
                        }
                        //outPortMap(myCFGLUT, "CDO", "open",false);

			string debug_signalName = "CFGLUT_ProductNo_"+ to_string(prNo) + "_LUTstage_" + to_string(stage) + "_LUT_No_" +to_string(LUT_No)+ "_and_" + to_string(LUT_No+1);
			outPortMap(myCFGLUT, "CDO", debug_signalName,true);
			
                        string instanceName = "CFGLUT_inst_ProductNo_"+ to_string(prNo) + "_LUTstage_" + to_string(stage) + "_LUT_No_" +to_string(LUT_No)+ "_and_" + to_string(LUT_No+1);
                        vhdl << myCFGLUT->primitiveInstance(instanceName,this);

                    }
                }
            }
        }
        addFullComment("BitHeap structure");

        unsigned int maxWeight = Const_bit_width+input_bit_width+additionalBitWidth;
        BitHeap *myBitHeap = new BitHeap(this,maxWeight-1);
        myBitHeap->addConstantOneBit(border+g-1);
        for(unsigned int i=0; i < bitHeapStack.size(); ++i)
        {
            if(bitHeapStackSize[i] > 0)
            {
                if((signed_calculation) && (i == bitHeapStack.size()-1)) // last element is the MSB and is signed in case of signed computing
                {
                    myBitHeap->addSignedBitVector(bitHeapStackShifts[i],bitHeapStack[i],bitHeapStackSize[i]);
                }
                else
                {
                    myBitHeap->addUnsignedBitVector(bitHeapStackShifts[i],bitHeapStack[i],bitHeapStackSize[i]);
                }
            }
        }
        myBitHeap->generateCompressorVHDL();

        // declaring output
        addOutput("result",output_bit_width);

        nextCycle();
        vhdl << tab << "result <= " << myBitHeap->getSumName() << "(" <<  myBitHeap->getMaxWeight() << " downto "  << (myBitHeap->getMaxWeight() - output_bit_width+1) <<  ");" << std::endl;


    }

	
    //void ConvolutionalCore::emulate(TestCase * tc){}
    //void ConvolutionalCore::buildStandardTestCases(TestCaseList * tcl){}

    string SOP_KCM::generateInitStringFor(int weight, unsigned int LUTNo)
    {
        int LutLImit = 1 << this->LUT_bit_width; //this calculates 2^LUT_bit_width usaly (LUT_bit_width = 4) it is 16.
        vector<int>ProductList;
        string initString;

        for (unsigned int i=0; i < (1 << LUT_bit_width);i++)
        {
            ProductList.push_back(i*weight);
        }

        int bitMask;

        bitMask = 1 << (LUTNo+1);
        //for(int i = 0; i < LutLImit; i++)
        for(int i = LutLImit-1; i >= 0; i--)
        {
            if(ProductList[i] & bitMask)
            {
                initString += "1";
            }
            else
            {
                initString += "0";
            }
        }

        bitMask = 1 << (LUTNo);

        //for(int i = 0; i < LutLImit; i++)
        for(int i = LutLImit-1; i >= 0; i--)
        {
            if(ProductList[i] & bitMask)
            {
                initString += "1";
            }
            else
            {
                initString += "0";
            }
        }

        initString = "\"" + initString + "\"";
        return initString;
    }



    OperatorPtr SOP_KCM::parseArguments(Target *target, vector<string> &args)
    {
        int param0, param1, param2, param3;
        UserInterface::parseInt(args, "inputWordSize", &param0); // param0 has a default value, this method will recover it if it doesnt't find it in args,
        UserInterface::parseInt(args, "constantWordSize", &param1);
        UserInterface::parseInt(args, "no_of_products", &param2);
        UserInterface::parseInt(args, "defaultProduct", &param3);
        return new SOP_KCM(target, param0, param1, param2, param3);
	}
	
    void SOP_KCM::registerFactory()
    {
        UserInterface::add("SOP_KCM", // name
                                             "My first SOP_KCM.", // description, string
											 "NeuralNetworks", // category, from the list defined in UserInterface.cpp
											 "", //seeAlso
											 // Now comes the parameter description string.
											 // Respect its syntax because it will be used to generate the parser and the docs
											 // Syntax is: a semicolon-separated list of parameterDescription;
											 // where parameterDescription is parameterName (parameterType)[=defaultValue]: parameterDescriptionString 
                                             "inputWordSize(int): Defined the input word size; \
                                             constantWordSize(int): Defined the Constant word size; \
                                             no_of_products(int)=1: Nombre of products; \
                                             defaultProduct(int)=5: the intial product for all inputs",
											 // More documentation for the HTML pages. If you want to link to your blog, it is here.
											 "Feel free to experiment with its code, it will not break anything in FloPoCo. <br> Also see the developper manual in the doc/ directory of FloPoCo.",
                                             SOP_KCM::parseArguments
											 ) ;
	}

}//namespace
