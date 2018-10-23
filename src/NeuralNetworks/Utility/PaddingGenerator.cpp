// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "PaddingGenerator.hpp"

#include "ModuloCounter.hpp"

using namespace std;
namespace flopoco {




    PaddingGenerator::PaddingGenerator(Target* target, unsigned int wordSize_, unsigned int windowSize_, unsigned int horizontalSize_, unsigned int verticalSize_, int padTop_, int strideH_, string padType_, bool genValidFinished_, bool buildForSerialCalculation_, int padBot_, int padLeft_, int padRight_, int strideV_) :
        Operator(target), wordSize(wordSize_), windowSize(windowSize_), horizontalSize(horizontalSize_), verticalSize(verticalSize_), padTop(padTop_), strideH(strideH_), padType(padType_), genValidFinished(genValidFinished_), buildForSerialCalculation(buildForSerialCalculation_), padBot(padBot_), padLeft(padLeft_), padRight(padRight_), strideV(strideV_) {

        // throw errors
        if(strideH<0)
        {
            stringstream e;
            e << "stride < 0 is not supported!";
            THROWERROR(e.str());
        }
        if(padTop<0)
        {
            stringstream e;
            e << "Padding < 0 is not supported!";
            THROWERROR(e.str());
        }
        if(windowSize==1)
        {
            stringstream e;
            e << "Window size=1 is not supported!";
            THROWERROR(e.str());
        }
        if(verticalSize<windowSize && padTop>0 && padBot!=0)
        {
            stringstream e;
            e << "This implementation of Padding doesn't work for this corner case (verticalSize<windowSize when padding top and bot): verticalSize=" << verticalSize << ", windowSize=" << windowSize << "!";
            THROWERROR(e.str());
        }

        // default values
        if(strideV<0)
        {
            strideV = strideH;
        }
        if(padBot<0)
        {
            padBot=((windowSize%2)==1?padTop:padTop-1);
        }
        if(padLeft<0)
        {
            padLeft=padTop;
        }
        if(padRight<0)
        {
            padRight=((windowSize%2)==1?padTop:padTop-1);
        }

        if(padRight+padLeft>windowSize || padTop+padBot>windowSize)
        {
            stringstream e;
            e << "Can't pad more than windowSize!";
            THROWERROR(e.str());
        }

        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="PaddingGenerator";

        // numeric std
        useNumericStd();

        // author
        setCopyrightString("Nicolai Fiege, 2017");

        // definition of the name of the operator
        ostringstream name;
        name << "PaddingGenerator_wordSize_" << wordSize << "_windowSize_" << windowSize << "_horizontalSize_" << horizontalSize << "_verticalSize_" << verticalSize << "_strideH_" << strideH << "_strideV_" << strideV << "_padTop_" << padTop << "_padBot_" << padBot << "_padLeft_" << padLeft << "_padRight_" << padRight << "_padType_" << padType << "_finished_valid_o_" << (genValidFinished==true?"true":"false") << ((buildForSerialCalculation==true)?("_for_serial"):("_for_parallel"));
        setName(name.str());

        // add in/out
        numberOfInputs=windowSize*windowSize;
        unsigned int inputCounter=0;
        for(unsigned int i=0; i<windowSize; i++)
        {
            vector <string> tmpvec;
            inputNames.push_back(tmpvec);
            outputNames.push_back(tmpvec);
            for(unsigned int j=0; j<windowSize; j++)
            {
                string tmp="X"+to_string(inputCounter);
                addInput(tmp, wordSize);
                inputNames[i].push_back(tmp);
                tmp="R"+to_string(inputCounter);
                addOutput(tmp, wordSize);
                outputNames[i].push_back(tmp);
                inputCounter++;
            }
        }
        if(buildForSerialCalculation==true)
        {
            addInput("newWeightsReady",1);
            declare("weightsArrivedAtConvCore",1);
            addOutput("readNewWeights",1);
            addOutput("convCoreWeightEnable",1);
        }

        addInput("validData_i",1);
        addInput("newStep",1);
        addOutput("getNewData",1);
        if(genValidFinished==true)
        {
            addOutput("finished",1);
            addOutput("validData_o",1);
        }


        // counters for horizontal and vertical position
        addConstant("hmax","integer",horizontalSize-1);
        addConstant("vmax","integer",verticalSize+padBot-(padRight>0?0:1));


        // calculate first and last valid output
        int firstValidOutputH;
        int firstValidOutputV;
        int lastValidOutputH;
        int lastValidOutputV;

        firstValidOutputH=(padLeft>0?windowSize-padLeft-1:windowSize-1);
        firstValidOutputV=(padTop>0?windowSize-padTop-1:windowSize-1);
        lastValidOutputH=(padRight>0?padRight-1:horizontalSize-1);
        lastValidOutputV=(padBot>0?(padRight>0?(verticalSize+padBot):(verticalSize+padBot-1)):(padRight>0?(verticalSize):(verticalSize-1))); //verticalSize+padBot

        declare("hcount",ceil(log2(horizontalSize+1)));
        declare("vcount",ceil(log2(verticalSize+padBot+2)));

        declare("waiting_for_newStep",1);

        // DEBUGGING START
        addOutput("hcount_debug",ceil(log2(horizontalSize+1)));
        addOutput("vcount_debug",ceil(log2(verticalSize+padBot+2)));
        this->vhdl << "hcount_debug <= hcount;" << endl;
        this->vhdl << "vcount_debug <= vcount;" << endl;
        // DEBUGGING END


        vhdl << "process(clk)" << endl;
        vhdl << "begin" << endl;
        vhdl << tab << "if(rising_edge(clk)) then" << endl;
        vhdl << tab << tab << "if(rst='1' or newStep=\"1\") then" << endl;
        vhdl << tab << tab << tab << "waiting_for_newStep <= \"0\";" << endl;
        vhdl << tab << tab << tab << "hcount <= (others => '0');" << endl;
        vhdl << tab << tab << tab << "vcount <= (others => '0');" << endl;
        vhdl << tab << tab << "else" << endl;
        vhdl << tab << tab << tab << "if(validData_i=\"1\") then" << endl;
        vhdl << tab << tab << tab << tab << "if(unsigned(hcount)<hmax) then" << endl;
        vhdl << tab << tab << tab << tab << tab << "hcount <= std_logic_vector(unsigned(hcount)+1);" << endl;
        vhdl << tab << tab << tab << tab << "else" << endl;
        vhdl << tab << tab << tab << tab << tab << "hcount <= (others => '0');" << endl;
        vhdl << tab << tab << tab << tab << tab << "if(unsigned(vcount)<vmax) then" << endl;
        vhdl << tab << tab << tab << tab << tab << tab << "vcount <= std_logic_vector(unsigned(vcount)+1);" << endl;
        vhdl << tab << tab << tab << tab << tab << "else" << endl;
        vhdl << tab << tab << tab << tab << tab << tab << "waiting_for_newStep <= \"1\";" << endl;
        vhdl << tab << tab << tab << tab << tab << "end if;" << endl;
        vhdl << tab << tab << tab << tab << "end if;" << endl;
        vhdl << tab << tab << tab << "end if;" << endl;
        vhdl << tab << tab << "end if;" << endl;
        vhdl << tab << "end if;" << endl;
        vhdl << "end process;" << endl;

        if(buildForSerialCalculation==true)
        {
            // tell the weight fetcher to read new weights from memory
            vhdl << "process(clk)" << endl;
            vhdl << "begin" << endl;
            vhdl << tab << "if(rising_edge(clk)) then" << endl;
            vhdl << tab << tab << "if(rst = '1' or newStep = \"1\") then" << endl;
            vhdl << tab << tab << tab << "weightsArrivedAtConvCore <= \"0\";" << endl;
            vhdl << tab << tab << tab << "readNewWeights <= \"1\";" << endl;
            vhdl << tab << tab << "elsif(weightsArrivedAtConvCore=\"0\" and newWeightsReady=\"1\") then" << endl;
            vhdl << tab << tab << tab << "weightsArrivedAtConvCore <= \"1\";" << endl;
            vhdl << tab << tab << tab << "readNewWeights <= \"0\";" << endl;
            vhdl << tab << tab << "end if;" << endl;
            vhdl << tab << "end if;" << endl;
            vhdl << "end process;" << endl;

            // tell the convolutional core to update its weights for the new calculation
            this->vhdl << "convCoreWeightEnable <= not weightsArrivedAtConvCore and newWeightsReady;" << endl;
        }


        this->vhdl << "getNewData <= \"1\" when " << declare("getNewData_guard",1) << "=\"1\" and waiting_for_newStep = \"0\" and (unsigned(vcount) < " << lastValidOutputV << " or ((unsigned(vcount) = " << lastValidOutputV << ") and (unsigned(hcount) < " << lastValidOutputH << ")))" << ((buildForSerialCalculation==true)?(" and (weightsArrivedAtConvCore=\"1\")"):("")) << " else \"0\";" << endl;
        vhdl << "process(clk)" << endl;
        vhdl << "begin" << endl;
        vhdl << tab << "if(rising_edge(clk)) then" << endl;
        vhdl << tab << tab << "if(rst='1' or newStep=\"1\") then" << endl;
        vhdl << tab << tab << tab << "getNewData_guard <= \"0\";" << endl;
        vhdl << tab << tab << "else" << endl;
        vhdl << tab << tab << tab << "getNewData_guard <= \"1\";" << endl;
        vhdl << tab << tab << "end if;" << endl;
        vhdl << tab << "end if;" << endl;
        vhdl << "end process;" << endl;

        if(genValidFinished==true)
        {
            if(strideH>1)
            {
                vhdl << "process(hcount)" << endl
                     << "begin" << endl
                     << tab << "if(unsigned(hcount)=" << firstValidOutputH-1 << ") then" << endl
                     << tab << tab << declare("modCounterHReset",1) << " <= \"1\";" << endl
                     << tab << "else" << endl
                     << tab << tab << "modCounterHReset <= \"0\";" << endl
                     << tab << "end if;" << endl
                     << "end process;" << endl;

                vhdl << declare("modCounterHTotalReset",1) << " <= modCounterHReset or newStep;" << endl;

                ModuloCounter* modCH = new ModuloCounter(target,strideH,true);
                this->addSubComponent(modCH);
                this->inPortMap(modCH,"enable","validData_i");
                this->inPortMap(modCH,"manualReset","modCounterHTotalReset");
                this->outPortMap(modCH,"counter","strideCounterH",true);
                vhdl << instance(modCH,"horizontalStrideCounter");
            }

            if(strideV>1)
            {

                vhdl << "process(vcount)" << endl
                     << "begin" << endl
                     << tab << "if(unsigned(vcount)=" << firstValidOutputV << ") then" << endl
                     << tab << tab << declare("modCounterVReset",1) << " <= \"1\";" << endl
                     << tab << "else" << endl
                     << tab << tab << "modCounterVReset <= \"0\";" << endl
                     << tab << "end if;" << endl
                     << "end process;" << endl;

                vhdl << "process(hcount)" << endl
                     << "begin" << endl
                     << tab << "if(unsigned(hcount)=" << firstValidOutputH-1 << ") then" << endl
                     << tab << tab << declare("strideCounterVEnable",1) << " <= \"1\";" << endl
                     << tab << "else" << endl
                     << tab << tab << "strideCounterVEnable <= \"0\";" << endl
                     << tab << "end if;" << endl
                     << "end process;" << endl;

                vhdl << declare("modCounterVTotalReset",1) << " <= modCounterVReset or newStep;" << endl;

                ModuloCounter* modCV = new ModuloCounter(target,strideV);
                this->addSubComponent(modCV);
                this->inPortMap(modCV,"enable","strideCounterVEnable");
                this->inPortMap(modCV,"manualReset","modCounterVTotalReset");
                this->outPortMap(modCV,"counter","strideCounterV",true);
                vhdl << instance(modCV,"verticalStrideCounter");
            }


            declare("finished_tmp",1);
            vhdl << "process(clk)" << endl;
            vhdl << "begin" << endl;
            vhdl << tab << "if(rising_edge(clk)) then" << endl;
            vhdl << tab << tab << "if(rst='1' or unsigned(vcount)<" << lastValidOutputV << ((lastValidOutputH!=0)?(" or unsigned(hcount)<" + to_string(lastValidOutputH)):("")) <<  ") then" << endl;
            vhdl << tab << tab << tab << "finished_tmp <= \"0\";" << endl;
            vhdl << tab << tab << "else" << endl;
            vhdl << tab << tab << tab << "finished_tmp <= \"1\";" << endl;
            vhdl << tab << tab << "end if;" << endl;
            vhdl << tab << "end if;" << endl;
            vhdl << "end process;" << endl;

            vhdl << "finished <= waiting_for_newStep or finished_tmp;" << endl;



            stringstream hcountIsInValidRange;
            hcountIsInValidRange << "(unsigned(hcount)>=" << firstValidOutputH;
            if(padRight>0)
                hcountIsInValidRange << " or unsigned(hcount)<=" << lastValidOutputH;
            else
                hcountIsInValidRange << " and unsigned(hcount)<=" << lastValidOutputH;
            hcountIsInValidRange << ")";

            vhdl << "process(hcount,vcount,validData_i" << ((strideH>1)?(",strideCounterH"):("")) << ((strideV>1)?(",strideCounterV"):("")) << ")" << endl;
            vhdl << "begin" << endl;
            vhdl << tab << "if(((unsigned(vcount)=" << firstValidOutputV << " and unsigned(hcount)>=" << firstValidOutputH << ") or (unsigned(vcount)>" << firstValidOutputV << " and " << "unsigned(vcount)<" << lastValidOutputV << " and (unsigned(hcount)>=" << firstValidOutputH << " or unsigned(hcount)<=" << lastValidOutputH << ")) or (unsigned(vcount)=" << lastValidOutputV << " and unsigned(hcount)<=" << lastValidOutputH << "))" << "and validData_i = \"1\"" << (strideH>1?" and unsigned(strideCounterH)=0":"") << (strideV>1?" and unsigned(strideCounterV)=0":"") << " and " << hcountIsInValidRange.str() << ") then" << endl;
            vhdl << tab<< tab << "validData_o <= \"1\";" << endl;
            vhdl << tab << "else" << endl;
            vhdl << tab << tab << "validData_o <= \"0\";" << endl;
            vhdl << tab << "end if;" << endl;
            vhdl << "end process;" << endl;
        }


        bool vpad = padTop>0 || padBot>0;
        bool hpad = padRight>0 || padLeft>0;
        if(vpad==true || hpad==true)
        {
            vhdl << "process(" << ((vpad==true)?"vcount":"") << ((vpad==true) && (hpad==true)?",":"") << ((hpad==true)?"hcount":"");
            for(unsigned int i=0; i<windowSize; i++)
            {
                for(unsigned int j=0; j<windowSize; j++)
                {
                    vhdl << "," << inputNames[i][j];
                }
            }
            vhdl << ")" << endl
                 << "begin" << endl;

            if(vpad==true)
            {
                vhdl << tab << tab << tab << "case to_integer(unsigned(vcount)) is" << endl;
            }

            // top padding
            if(padTop>0)
            {
                for(int vpadC=0; vpadC<(padRight>0?padTop+1:padTop);vpadC++)
                {
                    if(vpadC>0 || padLeft>0)
                    {
                        vhdl << tab << tab << tab << tab << "when " << firstValidOutputV+vpadC << " => " << endl; //vpad=0 => maximum vpadding
                    }
                    if(hpad==true)
                    {
                        vhdl << tab << tab << tab << tab << tab << "case to_integer(unsigned(hcount)) is" << endl;
                        if(padLeft>0)
                        {
                            // top left padding
                            unsigned int howMuchVPad=padTop-vpadC;
                            for(int hpadC=0; hpadC<padLeft;hpadC++)
                            {
                                vhdl << tab << tab << tab << tab << tab << tab << "when " << firstValidOutputH+hpadC << " => " << endl;
                                unsigned int howMuchHPad=padLeft-hpadC;
                                for(unsigned int y=0; y<windowSize; y++)
                                {
                                    for(unsigned int x=0; x<windowSize; x++)
                                    {
                                        vhdl << tab << tab << tab << tab << tab << tab << tab << outputNames[y][x] << " <= ";
                                        if(y<howMuchVPad && x<howMuchHPad)
                                        {
                                            string topLeftValidPixel = inputNames[howMuchVPad][howMuchHPad];
                                            vhdl << this->getPaddingValue(topLeftValidPixel);
                                        }
                                        else if(y<howMuchVPad)
                                        {
                                            string topValidPixel = inputNames[howMuchVPad][x];
                                            vhdl << this->getPaddingValue(topValidPixel);
                                        }
                                        else if(x<howMuchHPad)
                                        {
                                            string leftValidPixel = inputNames[y][howMuchHPad];
                                            vhdl << this->getPaddingValue(leftValidPixel);
                                        }
                                        else
                                        {
                                            vhdl << inputNames[y][x];
                                        }
                                        vhdl << ";" << endl;
                                    }
                                }
                            }
                        }
                        if(vpadC>0 && padRight>0)
                        {
                            // top right padding
                            unsigned int howMuchVPad=padTop+1-vpadC;
                            for(int hpadC=0; hpadC<padRight;hpadC++)
                            {
                                vhdl << tab << tab << tab << tab << tab << tab << "when " << padRight-hpadC-1 << " => " << endl;
                                unsigned int howMuchHPad=padRight-hpadC;
                                for(unsigned int y=0; y<windowSize; y++)
                                {
                                    for(unsigned int x=0; x<windowSize; x++)
                                    {
                                        vhdl << tab << tab << tab << tab << tab << tab << tab << outputNames[y][x] << " <= ";
                                        if(y<howMuchVPad && x>=windowSize-howMuchHPad)
                                        {
                                            string topRightValidPixel = inputNames[howMuchVPad][windowSize-1-howMuchHPad];
                                            vhdl << this->getPaddingValue(topRightValidPixel);
                                        }
                                        else if(y<howMuchVPad)
                                        {
                                            string topValidPixel = inputNames[howMuchVPad][x];
                                            vhdl << this->getPaddingValue(topValidPixel);
                                        }
                                        else if(x>=windowSize-howMuchHPad)
                                        {
                                            string rightValidPixel = inputNames[y][windowSize-1-howMuchHPad];
                                            vhdl << this->getPaddingValue(rightValidPixel);
                                        }
                                        else
                                        {
                                            vhdl << inputNames[y][x];
                                        }
                                        vhdl << ";" << endl;
                                    }
                                }
                            }
                        }

                        vhdl << tab << tab << tab << tab << tab << tab << "when others =>" << endl;
                    }

                    // purely top padding
                    unsigned int howMuchVPad=padTop-vpadC;
                    for(unsigned int y=0; y<windowSize; y++)
                    {
                        for(unsigned int x=0; x<windowSize; x++)
                        {
                            string topValidPixel = inputNames[howMuchVPad][x];
                            vhdl << tab << tab << tab << tab << tab << tab << tab << outputNames[y][x] << " <= ";
                            if(y<howMuchVPad)
                            {
                                vhdl << this->getPaddingValue(topValidPixel);
                            }
                            else
                            {
                                vhdl << inputNames[y][x];
                            }
                            vhdl << ";" << endl;
                        }
                    }

                    if(hpad==true)
                    {
                        vhdl << tab << tab << tab << tab << tab << "end case;" << endl;
                    }
                }
            }

            // bot padding
            if(padBot>0)
            {
                for(int vpadC=0; vpadC<(padRight>0?padBot+1:padBot);vpadC++)
                {
                    if(padLeft>0 || vpadC>0)
                    {
                        vhdl << tab << tab << tab << tab << "when " << lastValidOutputV-vpadC << " => " << endl; //vpad=0 => maximum vpadding
                    }
                    if(hpad==true)
                    {
                        vhdl << tab << tab << tab << tab << tab << "case to_integer(unsigned(hcount)) is" << endl;
                        if(vpadC>0 && padLeft>0)
                        {
                            // bot left padding
                            unsigned int howMuchVPad=padBot+1-vpadC;
                            for(int hpadC=0; hpadC<padLeft;hpadC++)
                            {
                                vhdl << tab << tab << tab << tab << tab << tab << "when " << firstValidOutputH+hpadC << " => " << endl;
                                unsigned int howMuchHPad=padLeft-hpadC;
                                for(unsigned int y=0; y<windowSize; y++)
                                {
                                    for(unsigned int x=0; x<windowSize; x++)
                                    {
                                        vhdl << tab << tab << tab << tab << tab << tab << tab << outputNames[y][x] << " <= ";
                                        if(y>=windowSize-howMuchVPad && x<howMuchHPad)
                                        {
                                            string botLeftValidPixel = inputNames[windowSize-1-howMuchVPad][howMuchHPad];
                                            vhdl << this->getPaddingValue(botLeftValidPixel);
                                        }
                                        else if(y>=windowSize-howMuchVPad)
                                        {
                                            string botValidPixel = inputNames[windowSize-1-howMuchVPad][x];
                                            vhdl << this->getPaddingValue(botValidPixel);
                                        }
                                        else if(x<howMuchHPad)
                                        {
                                            string leftValidPixel = inputNames[y][howMuchHPad];
                                            vhdl << this->getPaddingValue(leftValidPixel);
                                        }
                                        else
                                        {
                                            vhdl << inputNames[y][x];
                                        }
                                        vhdl << ";" << endl;
                                    }
                                }
                            }
                        }
                        if(vpadC<padBot && padRight>0)
                        {
                            // bot right padding
                            unsigned int howMuchVPad=padBot-vpadC;
                            for(int hpadC=0; hpadC<padRight;hpadC++)
                            {
                                vhdl << tab << tab << tab << tab << tab << tab << "when " << padLeft-hpadC-1 << " => " << endl;
                                unsigned int howMuchHPad=padRight-hpadC;
                                for(unsigned int y=0; y<windowSize; y++)
                                {
                                    for(unsigned int x=0; x<windowSize; x++)
                                    {
                                        vhdl << tab << tab << tab << tab << tab << tab << tab << outputNames[y][x] << " <= ";
                                        if(y>=windowSize-howMuchVPad && x>=windowSize-howMuchHPad)
                                        {
                                            string botRightValidPixel = inputNames[windowSize-1-howMuchVPad][windowSize-1-howMuchHPad];
                                            vhdl << this->getPaddingValue(botRightValidPixel);
                                        }
                                        else if(y>=windowSize-howMuchVPad)
                                        {
                                            string botValidPixel = inputNames[windowSize-1-howMuchVPad][x];
                                            vhdl << this->getPaddingValue(botValidPixel);
                                        }
                                        else if(x>=windowSize-howMuchHPad)
                                        {
                                            string rightValidPixel = inputNames[y][windowSize-1-howMuchHPad];
                                            vhdl << this->getPaddingValue(rightValidPixel);
                                        }
                                        else
                                        {
                                            vhdl << inputNames[y][x];
                                        }
                                        vhdl << ";" << endl;
                                    }
                                }
                            }
                        }
                        if(vpadC==padBot && padRight>0)
                        {
                            // purely right padding
                            for(int hpadC=0; hpadC<padRight;hpadC++)
                            {
                                vhdl << tab << tab << tab << tab << tab << tab << "when " << padRight-hpadC-1 << " => " << endl;
                                unsigned int howMuchHPad=padRight-hpadC;
                                for(unsigned int y=0; y<windowSize; y++)
                                {
                                    for(unsigned int x=0; x<windowSize; x++)
                                    {
                                        vhdl << tab << tab << tab << tab << tab << tab << tab << outputNames[y][x] << " <= ";
                                        if(x>=windowSize-howMuchHPad)
                                        {
                                            string rightValidPixel = inputNames[y][windowSize-1-howMuchHPad];
                                            vhdl << this->getPaddingValue(rightValidPixel);
                                        }
                                        else
                                        {
                                            vhdl << inputNames[y][x];
                                        }
                                        vhdl << ";" << endl;
                                    }
                                }
                            }
                        }

                        vhdl << tab << tab << tab << tab << tab << tab << "when others =>" << endl;
                    }

                    // purely bot padding
                    unsigned int howMuchVPad=padBot-vpadC+1; ///////
                    for(unsigned int y=0; y<windowSize; y++)
                    {
                        for(unsigned int x=0; x<windowSize; x++)
                        {
                            string botValidPixel = inputNames[windowSize-1-howMuchVPad][x];
                            vhdl << tab << tab << tab << tab << tab << tab << tab << outputNames[y][x] << " <= ";
                            if(y>=windowSize-howMuchVPad)
                            {
                                vhdl << this->getPaddingValue(botValidPixel);
                            }
                            else
                            {
                                vhdl << inputNames[y][x];
                            }
                            vhdl << ";" << endl;
                        }
                    }

                    if(hpad==true)
                    {
                        vhdl << tab << tab << tab << tab << tab << "end case;" << endl;
                    }
                }
            }

            // purely left/right
            if(vpad==true)
            {
                vhdl << tab << tab << tab << tab << "when others =>" << endl;
            }
            if(hpad==true)
            {
                vhdl << tab << tab << tab << tab << tab << "case to_integer(unsigned(hcount)) is" << endl;

                // purely left
                for(int hpadC=0; hpadC<padLeft; hpadC++)
                {
                    vhdl << tab << tab << tab << tab << tab << tab << "when " << firstValidOutputH+hpadC << " =>" << endl;
                    unsigned int howMuchHPad = padLeft-hpadC;
                    for(unsigned int y=0; y<windowSize; y++)
                    {
                        string leftValidPixel = inputNames[y][howMuchHPad];
                        for(unsigned int x=0; x<windowSize; x++)
                        {
                            vhdl << tab << tab << tab << tab << tab << tab << tab << outputNames[y][x] << " <= ";
                            if(x<howMuchHPad)
                            {
                                vhdl << this->getPaddingValue(leftValidPixel);
                            }
                            else
                            {
                                vhdl << inputNames[y][x];
                            }
                            vhdl << ";" << endl;
                        }
                    }
                }

                // purely right
                for(int hpadC=0; hpadC<padRight; hpadC++)
                {
                    vhdl << tab << tab << tab << tab << tab << tab << "when " << padRight-hpadC-1 << " =>" << endl;
                    unsigned int howMuchHPad = padRight-hpadC;
                    for(unsigned int y=0; y<windowSize; y++)
                    {
                        string rightValidPixel = inputNames[y][windowSize-1-howMuchHPad];
                        for(unsigned int x=0; x<windowSize; x++)
                        {
                            vhdl << tab << tab << tab << tab << tab << tab << tab << outputNames[y][x] << " <= ";
                            if(x>=windowSize-howMuchHPad)
                            {
                                vhdl << this->getPaddingValue(rightValidPixel);
                            }
                            else
                            {
                                vhdl << inputNames[y][x];
                            }
                            vhdl << ";" << endl;
                        }
                    }
                }
                vhdl << tab << tab << tab << tab << tab << tab << "when others =>" << endl;
            }

            // no padding at all
            for(unsigned int x=0; x<this->inputNames.size(); x++)
            {
                for(unsigned int y=0; y<this->inputNames[x].size(); y++)
                {
                    vhdl << tab << tab << tab << tab << tab << tab << tab << outputNames[x][y] << " <= " << inputNames[x][y] << ";" << endl;
                }
            }

            if(hpad==true)
            {
                vhdl << tab << tab << tab << tab << tab << "end case;" << endl;
            }
            if(vpad==true)
            {
                vhdl << tab << tab << tab << "end case;" << endl;
            }
            vhdl << "end process;" << endl;
        }
        else
        {
            // this component is only used to generate flag signals
            for(unsigned int x=0; x<this->inputNames.size(); x++)
            {
                for(unsigned int y=0; y<this->inputNames[x].size(); y++)
                {
                    vhdl << outputNames[x][y] << " <= " << inputNames[x][y] << ";" << endl;
                }
            }
        }
    }

    string PaddingGenerator::getPaddingValue(string name)
    {
        if(padType=="Zero")
        {
            return "(others => '0')";
        }
        else if(padType=="Value")
        {
            return name;
        }
        else
        {
            stringstream e;
            e << "Invalid Padding Type!";
            THROWERROR(e.str());
        }
    }

}//namespace flopoco
