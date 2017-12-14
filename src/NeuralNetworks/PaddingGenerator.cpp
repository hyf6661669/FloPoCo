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




    PaddingGenerator::PaddingGenerator(Target* target, unsigned int wordSize_, unsigned int windowSize_, unsigned int horizontalSize_, unsigned int verticalSize_, int padTop_, unsigned int stride_, string padType_, int padBot_, int padLeft_, int padRight_, bool genValidFinished_) :
        Operator(target), wordSize(wordSize_), windowSize(windowSize_), horizontalSize(horizontalSize_), verticalSize(verticalSize_), stride(stride_), padTop(padTop_), padBot(padBot_), padLeft(padLeft_), padRight(padRight_), padType(padType_), genValidFinished(genValidFinished_) {

        cout << "### Padding Generator Constructor!" << endl;
        // throw errors
        if(padTop<0)
        {
            stringstream e;
            e << "Padding < 0 is not supported!";
            THROWERROR(e);
        }
        cout << "### 1" << endl;
        if(windowSize==1)
        {
            stringstream e;
            e << "Window size=1 is not supported!";
            THROWERROR(e);
        }
        cout << "### 2" << endl;

        // default values
        if(padBot<0)
        {
            padBot=padTop;
        }
        cout << "### 3" << endl;
        if(padLeft<0)
        {
            padLeft=padTop;
        }
        cout << "### 4" << endl;
        if(padRight<0)
        {
            padRight=padTop;
        }

        cout << "### 5" << endl;
        cout << "### padRight=" << padRight << endl;
        cout << "### padLeft=" << padLeft << endl;
        cout << "### padTop=" << padTop << endl;
        cout << "### padBot=" << padBot << endl;
        cout << "### windowSize=" << windowSize << endl;
        cout << "### (double)windowSize=" << (double)windowSize << endl;
        cout << "### ((double)windowSize)/2)=" << ((double)windowSize)/2 << endl;
        cout << "### floor(((double)windowSize)/2)=" << floor(((double)windowSize)/2) << endl;
        if(padRight+padLeft>((double)windowSize) || padTop+padBot>((double)windowSize))
        {
            stringstream e;
            e << "Can't pad more than windowSize/2!";
            THROWERROR(e);
        }

        cout << "### 6" << endl;
        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="PaddingGenerator";

        // numeric std
        useNumericStd();

        // author
        setCopyrightString("Nicolai Fiege, 2017");

        cout << "### 7" << endl;
        // definition of the name of the operator
        ostringstream name;
        name << "PaddingGenerator_wordSize_" << wordSize << "_windowSize_" << windowSize << "_horizontalSize_" << horizontalSize << "_verticalSize_" << verticalSize << "_stride_" << stride << "_padTop_" << padTop << "_padBot_" << padBot << "_padLeft_" << padLeft << "_padRight_" << padRight << "_padType_" << padType << "_finishedSignal_" << (genValidFinished==true?"true":"false");
        setName(name.str());

        cout << "### Before adding in/out" << endl;
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
        addConstant("vmax","integer",verticalSize+padBot);


        // calculate first and last valid output
        int firstValidOutputH;
        int firstValidOutputV;
        int lastValidOutputH;
        int lastValidOutputV;

        firstValidOutputH=(padLeft>0?windowSize-padLeft-1:windowSize-1);
        firstValidOutputV=(padTop>0?windowSize-padTop-1:windowSize-1);
        lastValidOutputH=(padRight>0?padRight-1:horizontalSize-1);
        lastValidOutputV=(padBot>0?(padRight>0?(verticalSize+padBot):(verticalSize+padBot-1)):(padRight>0?(verticalSize):(verticalSize-1))); //verticalSize+padBot

        declare("hcount",ceil(log2(horizontalSize)));
        declare("vcount",ceil(log2(verticalSize+padBot+1)));
        vhdl << "process(clk)" << endl
             << "begin" << endl
             << tab << "if(rising_edge(clk)) then" << endl
             << tab << tab << "if(rst='1' or newStep=\"1\") then" << endl
             << tab << tab << tab << "hcount <= (others => '0');" << endl
             << tab << tab << tab << "vcount <= (others => '0');" << endl
             << tab << tab << "else" << endl
             << tab << tab << tab << "if(validData_i=\"1\") then" << endl
             << tab << tab << tab << tab << "if(unsigned(hcount)<hmax) then" << endl
             << tab << tab << tab << tab << tab << "hcount <= std_logic_vector(unsigned(hcount)+1);" << endl
             << tab << tab << tab << tab << "else" << endl
             << tab << tab << tab << tab << tab << "hcount <= (others => '0');" << endl
             << tab << tab << tab << tab << tab << "if(unsigned(vcount)<vmax) then" << endl
             << tab << tab << tab << tab << tab << tab << "vcount <= std_logic_vector(unsigned(vcount)+1);" << endl
             << tab << tab << tab << tab << tab << "else" << endl
             << tab << tab << tab << tab << tab << tab << "vcount <= (others => '0');" << endl
             << tab << tab << tab << tab << tab << "end if;" << endl
             << tab << tab << tab << tab << "end if;" << endl
             << tab << tab << tab << "end if;" << endl
             << tab << tab << "end if;" << endl
             << tab << "end if;" << endl
             << "end process;" << endl;


        // assign outputs
        vhdl << "process(vcount,hcount)" << endl
             << "begin" << endl
             << tab << "if(unsigned(vcount) < " << lastValidOutputV << " or ((unsigned(vcount) = " << lastValidOutputV << ") and (unsigned(hcount) < " << lastValidOutputH << "))) then" << endl
             << tab << tab << "getNewData <= \"1\";" << endl
             << tab << "else" << endl
             << tab << tab << "getNewData <= \"0\";" << endl
             << tab << "end if;" << endl
             << "end process;" << endl;

        if(genValidFinished==true)
        {
            vhdl << "process(hcount)" << endl
                 << "begin" << endl
                 << tab << "if(unsigned(hcount)=hmax) then" << endl
                 << tab << tab << declare("modCounterReset",1) << " <= \"1\";" << endl
                 << tab << "else" << endl
                 << tab << tab << "modCounterReset <= \"0\";" << endl
                 << tab << "end if;" << endl
                 << "end process;" << endl;

            vhdl << "process(vcount,hcount)" << endl
                 << "begin" << endl
                 << tab << "if(unsigned(vcount)=" << lastValidOutputV << " and unsigned(hcount)=" << lastValidOutputH << ") then" << endl
                 << tab << tab << "finished <= \"1\";" << endl
                 << tab << "else" << endl
                 << tab << tab << "finished <= \"0\";" << endl
                 << tab << "end if;" << endl
                 << "end process;" << endl;

            vhdl << declare("modCounterTotalReset",1) << " <= modCounterReset or newStep;" << endl;

            ModuloCounter* modC = new ModuloCounter(target,stride);
            this->inPortMap(modC,"enable","validData_i");
            this->inPortMap(modC,"manualReset","modCounterTotalReset");
            this->outPortMap(modC,"counter","strideCounter",true);


            vhdl << "process(vcount,hcount,strideCounter)" << endl
                 << "begin" << endl
                 << tab << "if(((unsigned(vcount)=" << firstValidOutputV << " and unsigned(hcount)>=" << firstValidOutputH << ") or (unsigned(vcount)>" << firstValidOutputV << " and " << "unsigned(vcount)<" << lastValidOutputV << " and (unsigned(hcount)>=" << firstValidOutputH << " or unsigned(hcount)<=" << lastValidOutputH << ")) or (unsigned(vcount)=" << lastValidOutputV << " and unsigned(hcount)<=" << lastValidOutputH << ")) and unsigned(strideCounter)=0) then" << endl
                 << tab << tab << "validData_o <= \"1\";" << endl
                 << tab << "else" << endl
                 << tab << tab << "validData_o <= \"0\";" << endl
                 << tab << "end if;" << endl
                 << "end process;" << endl;
        }


        cout << "### Before the padding-Process!" << endl;
        bool vpad = padTop>0 || padBot>0;
        bool hpad = padRight>0 || padLeft>0;
        if(vpad==true || hpad==true)
        {
            cout << "### PADDING IS TRU!" << endl;
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
                vhdl << tab << "case to_integer(unsigned(vcount)) is" << endl;
            }


            cout << "### Before Top padding" << endl;
            // top padding
            if(padTop>0)
            {
                for(int vpadC=0; vpadC<(padRight>0?padTop+1:padTop);vpadC++)
                {
                    if(vpadC>0 || padLeft>0)
                    {
                        vhdl << tab << tab << "when " << firstValidOutputV+vpadC << " => " << endl; //vpad=0 => maximum vpadding
                    }
                    if(hpad==true)
                    {
                        vhdl << tab << tab << tab << "case to_integer(unsigned(hcount)) is" << endl;
                        if(padLeft>0) // vpadC<padTop &&
                        {
                            // top left padding
                            unsigned int howMuchVPad=padTop-vpadC;
                            for(int hpadC=0; hpadC<padLeft;hpadC++)
                            {
                                vhdl << tab << tab << tab << tab << "when " << firstValidOutputH+hpadC << " => " << endl;
                                unsigned int howMuchHPad=padLeft-hpadC;
                                for(unsigned int y=0; y<windowSize; y++)
                                {
                                    for(unsigned int x=0; x<windowSize; x++)
                                    {
                                        vhdl << tab << tab << tab << tab << tab << outputNames[y][x] << " <= ";
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
                                vhdl << tab << tab << tab << tab << "when " << padRight-hpadC-1 << " => " << endl;
                                unsigned int howMuchHPad=padRight-hpadC;
                                for(unsigned int y=0; y<windowSize; y++)
                                {
                                    for(unsigned int x=0; x<windowSize; x++)
                                    {
                                        vhdl << tab << tab << tab << tab << tab << outputNames[y][x] << " <= ";
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

                        vhdl << tab << tab << tab << tab << "when others =>" << endl;
                    }

                    // purely top padding
                    unsigned int howMuchVPad=padTop-vpadC;
                    for(unsigned int y=0; y<windowSize; y++)
                    {
                        for(unsigned int x=0; x<windowSize; x++)
                        {
                            string topValidPixel = inputNames[howMuchVPad][x];
                            vhdl << tab << tab << tab << tab << tab << outputNames[y][x] << " <= ";
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
                        vhdl << tab << tab << tab << "end case;" << endl;
                    }
                }
            }

            cout << "### Before Bot padding" << endl;
            // bot padding
            if(padBot>0)
            {
                for(int vpadC=0; vpadC<(padRight>0?padBot+1:padBot);vpadC++)
                {
                    if(padLeft>0 || vpadC>0)
                    {
                        vhdl << tab << tab << "when " << lastValidOutputV-vpadC << " => " << endl; //vpad=0 => maximum vpadding
                    }
                    if(hpad==true)
                    {
                        vhdl << tab << tab << tab << "case to_integer(unsigned(hcount)) is" << endl;
                        if(vpadC>0 && padLeft>0)
                        {
                            // bot left padding
                            unsigned int howMuchVPad=padBot+1-vpadC;
                            for(int hpadC=0; hpadC<padLeft;hpadC++)
                            {
                                vhdl << tab << tab << tab << tab << "when " << firstValidOutputH+hpadC << " => " << endl;
                                unsigned int howMuchHPad=padLeft-hpadC;
                                for(unsigned int y=0; y<windowSize; y++)
                                {
                                    for(unsigned int x=0; x<windowSize; x++)
                                    {
                                        vhdl << tab << tab << tab << tab << tab << outputNames[y][x] << " <= ";
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
                                vhdl << tab << tab << tab << tab << "when " << padLeft-hpadC-1 << " => " << endl;
                                unsigned int howMuchHPad=padRight-hpadC;
                                for(unsigned int y=0; y<windowSize; y++)
                                {
                                    for(unsigned int x=0; x<windowSize; x++)
                                    {
                                        vhdl << tab << tab << tab << tab << tab << outputNames[y][x] << " <= ";
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

                        vhdl << tab << tab << tab << tab << "when others =>" << endl;
                    }

                    // purely bot padding
                    unsigned int howMuchVPad=padBot-vpadC+1; ///////
                    for(unsigned int y=0; y<windowSize; y++)
                    {
                        for(unsigned int x=0; x<windowSize; x++)
                        {
                            string botValidPixel = inputNames[windowSize-1-howMuchVPad][x];
                            vhdl << tab << tab << tab << tab << tab << outputNames[y][x] << " <= ";
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
                        vhdl << tab << tab << tab << "end case;" << endl;
                    }
                }
            }

            cout << "### Before Purely Left/Right padding" << endl;
            // purely left/right
            if(vpad==true)
            {
                vhdl << tab << tab << "when others =>" << endl;
            }
            if(hpad==true)
            {
                vhdl << tab << tab << tab << "case to_integer(unsigned(hcount)) is" << endl;

                // purely left
                for(int hpadC=0; hpadC<padLeft; hpadC++)
                {
                    vhdl << tab << tab << tab << tab << "when " << firstValidOutputH+hpadC << " =>" << endl;
                    unsigned int howMuchHPad = padLeft-hpadC;
                    for(unsigned int y=0; y<windowSize; y++)
                    {
                        string leftValidPixel = inputNames[y][howMuchHPad];
                        for(unsigned int x=0; x<windowSize; x++)
                        {
                            vhdl << tab << tab << tab << tab << tab << outputNames[y][x] << " <= ";
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
                    vhdl << tab << tab << tab << tab << "when " << padRight-hpadC-1 << " =>" << endl;
                    unsigned int howMuchHPad = padRight-hpadC;
                    for(unsigned int y=0; y<windowSize; y++)
                    {
                        string rightValidPixel = inputNames[y][windowSize-1-howMuchHPad];
                        for(unsigned int x=0; x<windowSize; x++)
                        {
                            vhdl << tab << tab << tab << tab << tab << outputNames[y][x] << " <= ";
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
                vhdl << tab << tab << tab << tab << "when others =>" << endl;
            }


            cout << "### Before No padding" << endl;
            // no padding at all
            for(unsigned int x=0; x<this->inputNames.size(); x++)
            {
                for(unsigned int y=0; y<this->inputNames[x].size(); y++)
                {
                    vhdl << tab << tab << tab << tab << tab << outputNames[x][y] << " <= " << inputNames[x][y] << ";" << endl;
                }
            }

            if(hpad==true)
            {
                vhdl << tab << tab << tab << "end case;" << endl;
            }
            if(vpad==true)
            {
                vhdl << tab << "end case;" << endl;
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
                    vhdl << inputNames[x][y] << " <= " << outputNames[x][y] << ";" << endl;
                }
            }
        }
        cout << "### Padding Constructor end" << endl;
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
            THROWERROR(e);
        }
    }

}//namespace flopoco
