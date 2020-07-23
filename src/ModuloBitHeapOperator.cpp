//
// Created by Annika Oeste on 20.06.20.
//
#include <iostream>
#include <sstream>
#include <bitset>
#include "gmp.h"
#include "mpfr.h"

#include "ModuloBitHeapOperator.hpp"

using namespace std;
namespace flopoco {


    ModuloBitHeapOperator::ModuloBitHeapOperator(Target* target, int wIn_, int wOut_, int modulo_) : Operator(target), wIn(wIn_), wOut(wOut_), modulo(modulo_) {
        srcFileName="ModuloBitHeapOperator";

        // definition of the name of the operator
        ostringstream name;
        name << "ModuloBitHeapOperator_" << wIn << "_" << wOut;
        setName(name.str()); // See also setNameWithFrequencyAndUID()
        // Copyright
        setCopyrightString("Annika Oeste, 2020");

        REPORT(INFO,"Declaration of ModuloBitHeapOperator \n");
        REPORT(DETAILED, "this operator has received three parameters: wIn " << wIn << " ,wOut " << wOut << " modulo " << modulo);        // modulo

        // debug message for developer
        REPORT(DEBUG,"debug of ModuloBitHeapOperator");

        // declaring inputs
        addInput ("X" , wIn, true);
        addConstant("M", "positive", modulo);
        // declaring output
        addOutput("S" , wOut);

        addFullComment("Start of vhdl generation");
        int bhStage = 0;
        // declares every bit in input as a signal
        for (int i = 0; i < wIn; ++i) {
            ostringstream bitName;
            bitName << "B_" << i << "_" << bhStage;
            vhdl << tab << declare(
                    bitName.str(), 1, false) << tab << "<= X" << of(i) << ";" << endl;

        }

        stage = 0;
        range.resize(1);
        range[stage].resize(2);
        range[stage][0] = 0;
        BitHeap *bitHeapTmp;

        bits.resize(15,vector<vector<unsigned long long> >(100,vector<unsigned long long>(64)));

        // initialize bits
        for (int i = 0; i < bits.size(); ++i) {
            for (int j = 0; j < bits[i].size(); ++j) {
                for (int k = 0; k < bits[i][j].size(); ++k) {
                    bits[i][j][k] = -1;
                }
            }
            }

        // fill bits for input
        for (int i = 0; i < wIn; ++i) {
            REPORT(DEBUG, "for loop i : " << i);
            bits[0][i][0] = i;
            range[stage][1] = range[stage][1] + (1 << i);
        }

        while (stage < bits.size() - 1) {
            int bitHeapHeight = getMaxHeight();
            REPORT(INFO, "bitHeapHeight: " << bitHeapHeight << " range: " << range[stage][1]);
            range.resize(stage+2);
            range[stage+1].resize(2);

            if (bitHeapHeight == 1 && range[stage][1] > (2*modulo)) {
                REPORT(DEBUG, "bitheapHeight == 1");
                applyPseudoCompressors();
            } else if (bitHeapHeight >= 2 || (bitHeapHeight == 1 && range[stage][1] < (2*modulo))) {
                REPORT(DEBUG, "bitheapHeight >= 2");
                ostringstream bhName;
                bhName << "modheap_" << bhStage;
                bitHeapTmp = new BitHeap(this, wOut-1, 0, bhName.str(), 1);

                // consider bits
                for (int j = 0; j < bits[stage].size(); ++j) {
                    for (int k = 0; k < bits[stage][j].size(); ++k) {
                        if (bits[stage][j][k] != -1) {
                            ostringstream signalName;

                            int bitOrigin = bits[stage][j][k] - (1 << (48));
                            if (bitOrigin - (1 << (32)) > 0) {
                                bitOrigin = bitOrigin - (1 << (32));
                                signalName << "B_" << bitOrigin << "_" << (bhStage);
                                REPORT(INFO,"bitheap " << bitHeapTmp->getName() << " signal use: " << signalName.str() << " bit value: " << bitOrigin);
                                bitHeapTmp->subtractSignal(signalName.str(), j);
                            } else {
                                signalName << "B_" << bitOrigin << "_" << (bhStage);
                                REPORT(INFO,"bitheap " << bitHeapTmp->getName() << " signal use: " << signalName.str() << " bit value: " << bitOrigin);
                                bitHeapTmp->addSignal(signalName.str(), j);
                            }
                        }
                    }
                }
                bitHeapTmp->startCompression();

                range[stage+1] = range[stage];
                REPORT(DEBUG, "range in next stage is: " << range[stage+1][1]);
                if (-modulo < range[stage+1][0] && range[stage+1][1] < (modulo*2)) {
                    REPORT(DEBUG, "modulo range reached in bitheap: break");
                    break;
                }

                // add sum bits to bits in next stage
                vector<vector<Bit*>> bhBits = bitHeapTmp->getBits();

                // add bits with weight to get correct int
                // this int then has to be put in the bits vector
                int bitHeapSum = 0;
                for (int i = 0; i < bhBits.size(); ++i) {
                    for (int j = 0; j < bhBits[i].size(); ++j) {
                        bitHeapSum += (1 << bhBits[i][j]->weight);
                    }
                }
                string sumString = bitset<10>(bitHeapSum).to_string();
                int msbSum = 0;
                for (int i = 0; i < sumString.size(); ++i) {
                    if ((sumString[i] - '0') == 1) {
                        msbSum = sumString.size() - 1 - i;
                        break;
                    }
                }
                for (int i = 0; i <= msbSum; ++i) {
                    bits[stage+1][i][0] = i;
                }

                for (int i = 0; i < wOut; ++i) {
                    ostringstream bitName;
                    bitName << "B_" << i << "_" << (bhStage+1);
                    vhdl << tab << declare(
                            bitName.str(), 1, false) << tab << "<=" << bitHeapTmp->getSumName() << of(i) << ";" << endl;
                }
                bhStage++;
            }
            stage++;
        }
        vhdl << tab << "S <= " << bitHeapTmp->getSumName(wOut-1,0) << " when ";
        vhdl << tab << bitHeapTmp->getSumName(wOut-1,0) << " < M else " << bitHeapTmp->getSumName(wOut-1,0) << " - M;" << endl;
        addFullComment("End of vhdl generation");
    }

    void ModuloBitHeapOperator::applyPseudoCompressors() {
        REPORT(INFO,"applyPseudoCompressors");
        int smsb = getMSBInStage();
        int mmsb = getModuloMSB();
        vector<vector<int>> rangeContrib(smsb+1,vector<int>(5));

        // fills the first two rows with the modulus and the congruent negative modulus
        for (int i = 0; i <= smsb; ++i) {
            rangeContrib[i][0] = (1 << (i)) % modulo;
            rangeContrib[i][1] = rangeContrib[i][0] - modulo;
            if (i == smsb && range[stage][0] > 0) {
                rangeContrib[i][0] = (-1 << (i)) % modulo;
                int power = (-1 << (i));
                //printf("negative power for i " + i + " : " + power);
                REPORT(INFO, "negative power for i " << i << " : " << power);
                rangeContrib[i][1] = rangeContrib[i][0] - modulo;
            }
        }

        int maxBhHeight = INT_MAX;
        int minBitsToBitHeap = INT_MAX;
        int rangeReached = 0;;
        vector<vector<int>> ranges(2,vector<int>(2));
        vector<vector<vector<unsigned long long>>> tmpBitHeap(smsb+1,vector<vector<unsigned long long> >(64,vector<unsigned long long>(2)));

        // initialize tmpBitHeap with -1
        for (int i = 0; i < tmpBitHeap.size(); ++i) {
            for (int j = 0; j < tmpBitHeap[i].size(); ++j) {
                for (int k = 0; k < tmpBitHeap[i][j].size(); ++k) {
                    tmpBitHeap[i][j][k] = -1;
                }
            }
        }

        for (int combi = 0; combi <= (1 << (smsb))-1; ++combi) {
            int combiMin = 0;
            int combiMax = 0;
            int bhHeight = 0;
            int bitsToBitHeap = 0;

            // calculate range max
            for (int i = 0; i <= smsb ; ++i) {
                if ((combi & (1 << (i - 1))) != 0 ) {
                    rangeContrib[i][2] = rangeContrib[i][0];
                    combiMax = combiMax + rangeContrib[i][0];
                } else {
                    rangeContrib[i][2] = rangeContrib[i][1];
                    combiMin = combiMin + rangeContrib[i][1];
                }
            }

            vector<vector<unsigned long long>> combiBitHeap(reqBitsForRange(combiMin, combiMax),vector<unsigned long long>(64));
            // initialize combiBitHeap with -1
            for (int i = 0; i < combiBitHeap.size(); ++i) {
                for (int j = 0; j < combiBitHeap[i].size(); ++j) {
                    combiBitHeap[i][j] = -1;
                }
            }
            vector<int> bhHeights(combiBitHeap.size());
            int constVec = 0;

            // for every column on bitheap
            for (int i = 0; i <= smsb; ++i) {
                int leadingZeroWeight = getLeadingZero(rangeContrib[i][2]);
                for (int j = 0; j < combiBitHeap.size(); ++j) {
                    if ((rangeContrib[i][2] & (1 << (j))) != 0 ) {
                        int bit = 0;
                        while (combiBitHeap[j][bit] != -1) {
                            bit++;
                        }
                        bitsToBitHeap++;
                        bhHeights[j] = bhHeights[j] + 1;
                        if (bhHeight < bhHeights[j]) {
                            bhHeight = bhHeights[j];
                        }
                        if (j <= leadingZeroWeight || 0 <= rangeContrib[i][2]) {
                            combiBitHeap[j][bit] = i + (1 << (48));
                        } else {
                            combiBitHeap[j][bit] = i + (1 << (48)) + (1 << (32));
                            constVec = (constVec + (1 << (combiBitHeap.size())) - 1 -
                                    (1 << (leadingZeroWeight)) -1) & (1 << (smsb));
                            break;
                        }
                    }
                }
            }
            for (int i = 0; i < combiBitHeap.size(); ++i) {
                if (constVec & (1 << (i)) != 0) {
                    int bit = 0;
                    while (combiBitHeap[i][bit] != -1) {
                        bit++;
                    }
                    combiBitHeap[i][bit] = i + (1 << (48)) + 5 * (1 << (32));
                    bitsToBitHeap++;
                    bhHeights[i] = bhHeights[i] + 1;
                    if (bhHeight < bhHeights[i]) {
                        bhHeight = bhHeights[i];
                    }
                }
            }
            // final stage of pseudocompression
            if (-modulo < combiMin && combiMax < modulo) {
                rangeReached = 1;
                ranges[1] = {combiMin, combiMax};
                for (int i = 0; i < rangeContrib.size(); ++i) {
                    rangeContrib[i][4] = rangeContrib[i][2];
                }
                for (int i = 0; i < combiBitHeap.size(); ++i) {
                    for (int j = 0; j < combiBitHeap[i].size(); ++j) {
                        tmpBitHeap[i][j][1] = combiBitHeap[i][j];
                    }
                }
            }

            // fill bits in temporary bitheap
            if (bitsToBitHeap < minBitsToBitHeap && bhHeight <= maxBhHeight) {
                minBitsToBitHeap = bitsToBitHeap;
                ranges[0] = {combiMin, combiMax};
                for (int i = 0; i < rangeContrib.size(); ++i) {
                    rangeContrib[i][3] = rangeContrib[i][2];
                }
                maxBhHeight = bhHeight;
                for (int i = 0; i < combiBitHeap.size(); ++i) {
                    for (int j = 0; j < combiBitHeap[i].size(); ++j) {
                        tmpBitHeap[i][j][0] = combiBitHeap[i][j];
                    }
                }
            }
        }

        // fill bits from tmp bitheap in bits
        range[stage+1] = ranges[rangeReached];
        int reqCols = reqBitsForRange(ranges[rangeReached][0], ranges[rangeReached][1]);
        for (int i = 0; i < tmpBitHeap.size(); ++i) {
            for (int j = 0; j < reqCols; ++j) {
                bits[stage+1][i][j] = tmpBitHeap[i][j][rangeReached];
            }
        }
    }

    int ModuloBitHeapOperator::getMSBInStage() {
        int smsb = 0;
        for (int w = 0; w < 100; ++w) {
            if (bits[stage][w][0] != -1) {
                smsb = w;
            }
        }
        return smsb;
    }

    int ModuloBitHeapOperator::getMaxHeight() {
        int height = 0;
        for (int c = 0; c < bits[stage].size(); ++c) {
            int currentHeight = 0;
            for (int b = 0; b < bits[stage][c].size(); ++b) {
                if (bits[stage][c][b] != -1) {
                    currentHeight++;
                }
            }
            if (currentHeight > height) {
                height = currentHeight;
            }
        }
        return height;
    }

    int ModuloBitHeapOperator::getModuloMSB() {
        int mmsb = 1;
        for (int w = 0; w < 31; ++w) {
            if (modulo & (1 << (w)) != 0) {
                mmsb = w;
            }
        }
        return mmsb;
    }

    int ModuloBitHeapOperator::reqBitsForRange(int min, int max) {
        int bit = 0;
        while((max >> bit) != 0) {
            bit++;
        }
        if (min < 0) {
            int negBit = 0;
            while((min >> negBit) != -1) {
                negBit++;
            }
            if (bit < negBit) {
                bit = negBit + 1;
            } else {
                bit++;
            }
        }
        return bit;
    }

    int ModuloBitHeapOperator::getLeadingZero(int value) {
        int msz = 0;
        for (int w = 0; w < 32; ++w) {
            if (value & (1 << (w)) == 0) {
                msz = w;
            }
        }
        return msz;
    }

    void ModuloBitHeapOperator::emulate(TestCase * tc) {
        /* This function will be used when the TestBench command is used in the command line
           we have to provide a complete and correct emulation of the operator, in order to compare correct output generated by this function with the test input generated by the vhdl code */
        /* first we are going to format the entries */
        mpz_class sx = tc->getInputValue("X");

        /* then we are going to manipulate our bit vectors in order to get the correct output*/
        mpz_class sr;
        sr = sx % modulo;

        /* at the end, we indicate to the TestCase object what is the expected
           output corresponding to the inputs */
        tc->addExpectedOutput("S",sr);
    }

    void ModuloBitHeapOperator::buildStandardTestCases(TestCaseList * tcl) {
        TestCase *tc;
        tc = new TestCase(this);
        tc->addInput("X", mpz_class(15));
        emulate(tc);
        tcl->add(tc);
    }

    OperatorPtr ModuloBitHeapOperator::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args) {
        int wIn, wOut, modulo;
        UserInterface::parseInt(args, "wIn", &wIn);
        UserInterface::parseInt(args, "wOut", &wOut);
        UserInterface::parseInt(args, "modulo", &modulo);
        return new ModuloBitHeapOperator(target, wIn, wOut, modulo);
    }

    void ModuloBitHeapOperator::registerFactory(){
        UserInterface::add("ModuloBitHeapOperator", // name
                           "An Operator to test the modulo calculation in combination with a bitHeap", // description, string
                           "Miscellaneous", // category, from the list defined in UserInterface.cpp
                           "", //seeAlso
                // Now comes the parameter description string.
                // Respect its syntax because it will be used to generate the parser and the docs
                // Syntax is: a semicolon-separated list of parameterDescription;
                // where parameterDescription is parameterName (parameterType)[=defaultValue]: parameterDescriptionString
                           "wIn(int)=16: A first parameter - the input size; \
                            wOut(int): A second parameter - the output size; \
                        modulo(int): modulo",
                // More documentation for the HTML pages. If you want to link to your blog, it is here.
                           "See the developer manual in the doc/ directory of FloPoCo.",
                           ModuloBitHeapOperator::parseArguments
        ) ;
    }
}//namespace
