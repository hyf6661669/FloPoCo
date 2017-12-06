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
#include "../PrimitiveComponents/Xilinx/Xilinx_CFGLUTShadow.hpp"
#include "../BitHeap/BitHeap.hpp"
#include <math.h>

using namespace std;
namespace flopoco {




    SOP_KCM::SOP_KCM(Target* target,int inputWordSize,int constantWordSize, int parameterNo, int _defaultConstant, bool useShadowLUTs, bool useFaithfulRounding) : Operator(target), input_bit_width(inputWordSize), Const_bit_width(constantWordSize), No_of_Products(parameterNo)
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
        faithful_rounding = useFaithfulRounding;
        allow_half_start_LUT = false;

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

        bool halfLUTusage_justO5 =false;// will be set wen a half used Lut is generated
        bool halfLUTusage_justO6 =false;// will be set wen a half used Lut is generated


        int outputBits = Const_bit_width + input_bit_width + additionalBitWidth;

        if(faithful_rounding)
        {
            output_bit_width = inputWordSize;
        }
        else
        {
            output_bit_width = outputBits;
        }

        addInput ("LUT_Config_clk_enable");

        vector<string> bitHeapStack;
        vector<unsigned int> bitHeapStackShifts;
        vector<int> bitHeapStackSize;
        bitHeapStackShifts.clear();
        bitHeapStack.clear();
        bitHeapStackSize.clear();



        int border =  outputBits - output_bit_width - g;;




        string differentLUTContent_output;
        int differentLUTContentCounter = 0;
        vector<int> differentLUTContent_counterLSB;
        vector<int> differentLUTContent_counterMSB;

        addFullComment("Parameter:");
        addComment(join("useFaithfulRounding=",useFaithfulRounding));
        addComment(join("useShadowLUTs=",useShadowLUTs));
        addComment(join("input_bit_width=",input_bit_width));
        addComment(join("Const_bit_width=",Const_bit_width));
        addComment(join("output_bit_width=",output_bit_width));
        addComment(join("outputBits=",outputBits));
        addComment(join("border=",border));
        addComment(join("additionalBitWidth=",additionalBitWidth));
        addComment(join("No_of_Products=",No_of_Products));
        addComment(join("No_of_stages=",No_of_stages));
        addComment(join("g=",g));

        nextCycle();

        for (int prNo = 0; prNo < No_of_Products; ++prNo)
        {
            differentLUTContent_output = "\n --Lut Content for Product " + to_string(prNo);
            differentLUTContent_counterLSB.clear();
            differentLUTContent_counterLSB.resize(LUT_per_stage);
            differentLUTContent_counterMSB.clear();
            differentLUTContent_counterMSB.resize(LUT_per_stage);


            string inputSignalName = join("X",prNo);
            addInput ( inputSignalName, input_bit_width);
            inputSignalName+="(";

            addFullComment(join("Product No ",prNo));

            //string configurationStreamSignalName = "cdi_no_" + to_string(prNo);
            //addInput(configurationStreamSignalName,ceil(((float)LUT_per_stage)/2)*No_of_stages);

            string configurationStreamSignalLSBName = "cdi_LSB_no_" + to_string(prNo);
            string configurationStreamSignalMSBName = "cdi_MSB_no_" + to_string(prNo);


            {
                cdi_bit_width = LUT_per_stage+1; // because of the last LUT which is just half used.

                addInput(configurationStreamSignalLSBName,cdi_bit_width);
                addInput(configurationStreamSignalMSBName,cdi_bit_width);
            }

            configurationStreamSignalLSBName += "(";
            configurationStreamSignalMSBName += "(";

            for(int stage = 0; stage < No_of_stages; ++stage)
            {
                string configurationStreamSignalName;
                if(stage == No_of_stages-1)
                {
                    configurationStreamSignalName = configurationStreamSignalMSBName; //the configuration stream is different for the las stage
                }
                else
                {
                    configurationStreamSignalName = configurationStreamSignalLSBName;// all other stages for the same constant can share the same stream.
                }

                addFullComment(join("Stage No ",stage));
                int Lut_start_Counter = 0;
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

                    int Lut_start_Counter_correction_Value=0;
                    if(((Lut_start_Counter % 2) != 0) && (allow_half_start_LUT == false))
                    {
                        Lut_start_Counter_correction_Value = 1;
                        differentLUTContent_counterLSB.resize(LUT_per_stage+1);
                        differentLUTContent_counterMSB.resize(LUT_per_stage+1);
                    }
                    for(int LUT_No = Lut_start_Counter-Lut_start_Counter_correction_Value; LUT_No< LUT_per_stage; LUT_No += 2)
                    {
                        bool switch_O5_and_O6 = false;
                        vhdl << std::endl;

                        Operator *myCFGLUT;
                        if(useShadowLUTs)
                        {
                            myCFGLUT = new Xilinx_CFGLUTShadow(target);
                        }
                        else
                        {
                            myCFGLUT = new Xilinx_CFGLUT5(target);
                            inPortMapCst(myCFGLUT, "clk","clk");
                            if(stage == No_of_stages-1)// if it is the MSB Lut stage
                            {
                                ((Xilinx_CFGLUT5*) myCFGLUT)->setGeneric( "init", generateInitStringFor(defaultConstant,LUT_No,1) );
                            }
                            else
                            {
                                ((Xilinx_CFGLUT5*) myCFGLUT)->setGeneric( "init", generateInitStringFor(defaultConstant,LUT_No,0) );
                            }
                        }
                        addToGlobalOpList(myCFGLUT);


                        inPortMap(myCFGLUT, "ce","LUT_Config_clk_enable");
                        for(int i = 0; i <= LUT_bit_width; ++i)
                        {
                            if(i == LUT_bit_width)
                            {
                                inPortMapCst(myCFGLUT, join("i", to_string(i)),"'1'");// the highest bit is true to use the 5 input LUT as two 4 input Luts
                            }
                            else if((i+stage*LUT_bit_width) < input_bit_width)
                            {
                                inPortMap(myCFGLUT, join("i", to_string(i)),join(inputSignalName,to_string(i+stage*LUT_bit_width),")"));
                            }
                            else
                            {
                                inPortMapCst(myCFGLUT, join("i", to_string(i)),"'0'");
                            }
                        }

                        string outputSignalNameLUT;
                        if((Lut_start_Counter_correction_Value == 0) || (LUT_No > Lut_start_Counter))
                        {
                            outputSignalNameLUT = outputSignalName + "(" + to_string(LUT_No-Lut_start_Counter) +")";
                            outPortMap(myCFGLUT, "o6",outputSignalNameLUT,false);// MH switch 5 6
                            if(stage < (No_of_stages-1))
                            {
                                differentLUTContent_counterLSB[LUT_No]++;
                            }
                            else
                            {
                                differentLUTContent_counterMSB[LUT_No]++;
                            }
                        }
                        else
                        { // switch the LUT output Bits O5 and O6 because using just o5 is not possible... so youst o6 is used and o5 stay unused...
                            switch_O5_and_O6 = true;
                            //halfLUTusage_justO5 = true;

                            outputSignalNameLUT = outputSignalName + "(" + to_string(LUT_No+1-Lut_start_Counter) +")";
                            outPortMap(myCFGLUT, "o6",outputSignalNameLUT,false);
                            if(stage < (No_of_stages-1))
                            {
                                differentLUTContent_counterLSB[(LUT_No)]++;
                            }
                            else
                            {
                                differentLUTContent_counterMSB[(LUT_No)]++;
                            }
                        }
			
			
                        if((switch_O5_and_O6 == false) && (LUT_No+1 < LUT_per_stage))//if the lut content wasnt switched and the Lut is realy used
                        {
                            outputSignalNameLUT = outputSignalName + "(" + to_string(LUT_No+1-Lut_start_Counter) +")";
                            outPortMap(myCFGLUT, "o5",outputSignalNameLUT,false);
                            if(stage < (No_of_stages-1))
                            {
                                differentLUTContent_counterLSB[(LUT_No+1)]++;
                            }
                            else
                            {
                                differentLUTContent_counterMSB[(LUT_No+1)]++;
                            }
                        }
                        else
                        {
                            halfLUTusage_justO6 = true;
                            outPortMap(myCFGLUT, "o5","open",false);
                        }
                        //outPortMap(myCFGLUT, "CDO", "open",false);

                        if(switch_O5_and_O6) // if the output Ports are switched a also modyfied configuration stream is nacessary to compute te correct output.
                        {
                            inPortMap(myCFGLUT, "cdi",join(configurationStreamSignalName,to_string(LUT_No),")"));
                        }
                        else
                        {
                            inPortMap(myCFGLUT, "cdi",join(configurationStreamSignalName,to_string(LUT_No+1),")"));
                        }

                        string debug_signalName = "CFGLUT_ProductNo_"+ to_string(prNo) + "_LUTstage_" + to_string(stage) + "_LUT_No_" +to_string(LUT_No)+ "_and_" + to_string(LUT_No+1);
                        outPortMap(myCFGLUT, "cdo", debug_signalName,true);
			
                        string instanceName = "CFGLUT_inst_ProductNo_"+ to_string(prNo) + "_LUTstage_" + to_string(stage) + "_LUT_No_" +to_string(LUT_No)+ "_and_" + to_string(LUT_No+1);

                        if(useShadowLUTs)
                        {
                            vhdl << instance(myCFGLUT,instanceName);
                        }
                        else
                        {
                            vhdl << ((Xilinx_CFGLUT5*) myCFGLUT)->primitiveInstance(instanceName,this);
                        }
                    }
                }
            }

            differentLUTContent_output += "\n--Product No: "+ to_string(prNo) +" differentLUTContent_counterLSB:";
            for(int i=differentLUTContent_counterLSB.size()-1; i >=0 ;--i)
            {
                differentLUTContent_output += to_string(differentLUTContent_counterLSB[i]) + " ";
                if(i%2==0)
                    differentLUTContent_output += "| ";
            }
            differentLUTContent_output += "\n--Product No: "+ to_string(prNo) +" differentLUTContent_counterMSB:";
            for(int i=differentLUTContent_counterMSB.size()-1; i >=0 ;--i)
            {
                differentLUTContent_output += to_string(differentLUTContent_counterMSB[i]) + " ";
                if(i%2==0)
                    differentLUTContent_output += "| ";
            }

            for(unsigned int i=0; i < differentLUTContent_counterMSB.size(); i += 2)
            {
                if ((differentLUTContent_counterLSB[i] != 0) || (differentLUTContent_counterLSB[i+1] != 0))
                {
                    differentLUTContentCounter++;
                }
                if ((differentLUTContent_counterMSB[i] != 0) || (differentLUTContent_counterMSB[i+1] != 0))
                {
                    differentLUTContentCounter++;
                }
            }

        }

        nextCycle();

        addFullComment("BitHeap structure");

        unsigned int maxWeight = Const_bit_width+input_bit_width+additionalBitWidth;
        BitHeap *myBitHeap = new BitHeap(this,maxWeight-1);

        if(faithful_rounding)
        {
            myBitHeap->addConstantOneBit(border+g-1);
        }

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
        vhdl << tab << "result <= " << myBitHeap->getSumName() << "(" <<  myBitHeap->getMaxWeight() << " downto "  << ((myBitHeap->getMaxWeight() - output_bit_width) + 1) <<  ");" << std::endl;

        vhdl << differentLUTContent_output <<std::endl;
        vhdl << "--" << name.str() << " differentLUTContentCounter:" << differentLUTContentCounter << " halfLUTusageO5:" << halfLUTusage_justO5 << " halfLUTusageO6:" << halfLUTusage_justO6 << std::endl;


    }

	
    //void ConvolutionalCore::emulate(TestCase * tc){}
    //void ConvolutionalCore::buildStandardTestCases(TestCaseList * tcl){}

    string SOP_KCM::generateInitStringFor(int weight, unsigned int LUTNo, bool MSBLUT)
    {
        int LutLImit = 1 << this->LUT_bit_width; //this calculates 2^LUT_bit_width usaly (LUT_bit_width = 4) it is 16.
        vector<int>ProductList;
        string initString;

        ProductList.resize(16);
        switch(MSBLUT)
        {
        case 0:
            ProductList[ 0]=(   0 * weight);
            ProductList[ 1]=(   1 * weight);
            ProductList[ 2]=(   2 * weight);
            ProductList[ 3]=(   3 * weight);
            ProductList[ 4]=(   4 * weight);
            ProductList[ 5]=(   5 * weight);
            ProductList[ 6]=(   6 * weight);
            ProductList[ 7]=(   7 * weight);
            ProductList[ 8]=(   8 * weight);
            ProductList[ 9]=(   9 * weight);
            ProductList[10]=(  10 * weight);
            ProductList[11]=(  11 * weight);
            ProductList[12]=(  12 * weight);
            ProductList[13]=(  13 * weight);
            ProductList[14]=(  14 * weight);
            ProductList[15]=(  15 * weight);
        break;
        case 1:
            ProductList[ 0]=(   0 * weight);
            ProductList[ 1]=(   1 * weight);
            ProductList[ 2]=(   2 * weight);
            ProductList[ 3]=(   3 * weight);
            ProductList[ 4]=(   4 * weight);
            ProductList[ 5]=(   5 * weight);
            ProductList[ 6]=(   6 * weight);
            ProductList[ 7]=(   7 * weight);
            ProductList[ 8]=(  -8 * weight);
            ProductList[ 9]=(  -7 * weight);
            ProductList[10]=(  -6 * weight);
            ProductList[11]=(  -5 * weight);
            ProductList[12]=(  -4 * weight);
            ProductList[13]=(  -3 * weight);
            ProductList[14]=(  -2 * weight);
            ProductList[15]=(  -1 * weight);
        }

        for(unsigned int i=0; i < ProductList.size(); i++)
        {
            std::cout << "MH: ProductList[" << i << "]: "<< ProductList[i] << std::endl;
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
        bool param4, param5;
        UserInterface::parseInt(args, "inputWordSize", &param0); // param0 has a default value, this method will recover it if it doesnt't find it in args,
        UserInterface::parseInt(args, "constantWordSize", &param1);
        UserInterface::parseInt(args, "no_of_products", &param2);
        UserInterface::parseInt(args, "defaultProduct", &param3);
        UserInterface::parseBoolean(args, "useShadowLUTs", &param4);
        UserInterface::parseBoolean(args, "useFaithfulRounding", &param5);
        //useShadowLUTs(bool)=false Switch to choose the CFGLUT implementation with ore without Shadow LUTs"
        return new SOP_KCM(target, param0, param1, param2, param3, param4, param5);
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
                                     no_of_products(int)=1: Number of products; \
                                     useShadowLUTs(bool)=false: Using two different LUts to configure one while using the other. The used one on the other can switched in one cycle; \
                                     useFaithfulRounding(bool)=true: reduce the output wordsize with facefoul rounding to the input wordsize; \
                                     defaultProduct(int)=32: the intial product for all inputs;",
                                     // More documentation for the HTML pages. If you want to link to your blog, it is here.
                                     "Feel free to experiment with its code, it will not break anything in FloPoCo. <br> Also see the developper manual in the doc/ directory of FloPoCo.",
                                     SOP_KCM::parseArguments
                                             ) ;
    }

    int SOP_KCM::get_cdi_bit_with()
    {
        return this->cdi_bit_width;
    }

    int SOP_KCM::get_input_bit_with()
    {
        return this->input_bit_width;
    }

    int SOP_KCM::get_output_bit_with()
    {
        return this->output_bit_width;
    }

}//namespace
