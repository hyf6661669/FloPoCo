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


    ModuloBitHeapOperator::ModuloBitHeapOperator(OperatorPtr parentOp, Target* target, int wIn_, int modulo_, string method_)
    : Operator(parentOp,  target), wIn(wIn_), modulo(modulo_), method(method_) {
        srcFileName="ModuloBitHeapOperator";

        // definition of the name of the operator
        ostringstream name;
        name << "ModuloBitHeapOperator_" << wIn << "_" << modulo;
        setName(name.str()); // See also setNameWithFrequencyAndUID()
        // Copyright
        setCopyrightString("Annika Oeste, 2020");
        useNumericStd();

        REPORT(INFO,"Declaration of ModuloBitHeapOperator \n");
        REPORT(DETAILED, "this operator has received three parameters: wIn " << wIn << " modulo " << modulo << " method " << method);

        // debug message for developer
        REPORT(DEBUG,"debug of ModuloBitHeapOperator");

        // declaring inputs
        addInput ("X" , wIn, true);
        addConstant("M", "positive", modulo);
        // declaring output
        int outPutSize = getModuloMSB() + 1;
        addOutput("S" , outPutSize);

        int bhStage = 0;
        // declares every bit in input as a signal
        for (int i = 0; i < wIn; ++i) {
            ostringstream bitName;
            bitName << "B_" << i << "_" << bhStage;
            vhdl << tab << declare(
                    bitName.str(), 1, true) << tab << "<= X" << range(i, i) << ";" << endl;

        }

        stage = 0;
        rangeBH.resize(1);
        rangeBH[stage].resize(2);
        rangeBH[stage][0] = 0;
        BitHeap *bitHeapTmp;
        int bitHeapReqBits;
        int bitHeapUsedBits;

        bits.resize(30,vector<vector<long long> >(100,vector<long long>(64)));

        // initialize bits
        for (int i = 0; i < bits.size(); ++i) {
            for (int j = 0; j < bits[i].size(); ++j) {
                for (int k = 0; k < bits[i][j].size(); ++k) {
                    bits[i][j][k] = -1;
                }
            }
        }

        long long oneLL = static_cast<long long>(1);
        // fill bits for input
        for (int i = 0; i < wIn; ++i) {
            bits[0][i][0] = i;
            rangeBH[stage][1] = rangeBH[stage][1] + (oneLL << i);
        }

        while (stage < bits.size() - 1) {
            int bitHeapHeight = getMaxHeight();
            REPORT(INFO, "bitHeapHeight: " << bitHeapHeight << " upper: " << rangeBH[stage][1] << " lower: " << rangeBH[stage][0]);
            rangeBH.resize(stage+2);
            rangeBH[stage+1].resize(2);

            if (bitHeapHeight == 1 && (rangeBH[stage][1] >= modulo || rangeBH[stage][0] < -modulo)) {
                REPORT(INFO, "bitheapHeight == 1");
                if (method == "compl") {
                    applyPseudoCompressors();
                } else if (method == "lMinBits") {
                    applyPseudoCompressorsMinBits();
                    if (checkForRepetition()) {
                        for (int i = 0; i < bits[stage+1].size(); ++i) {
                            for (int j = 0; j < bits[stage+1][i].size(); ++j) {
                                bits[stage+1][i][j] = -1;
                            }
                        }
                        applyPseudoCompressorsMinRange();
                    }
                } else if (method == "minRange") {
                    applyPseudoCompressorsMinRange();
                } else if (method == "minRangeWeight") {
                    applyPseudoCompressorsMinRangeWeighted();
                    if (checkForRepetition()) {
                        for (int i = 0; i < bits[stage+1].size(); ++i) {
                            for (int j = 0; j < bits[stage+1][i].size(); ++j) {
                                bits[stage+1][i][j] = -1;
                            }
                        }
                        applyPseudoCompressorsMinRange();
                    }
                } else {
                    REPORT(INFO, "invalid method name " << method);
                    break;
                }
            } else if (bitHeapHeight >= 2 || (bitHeapHeight == 1 && (rangeBH[stage][1] < modulo || rangeBH[stage][0] >= -modulo))) {
                REPORT(INFO, "bitheapHeight >= 2");
                ostringstream bhName;
                bhName << "modheap_" << bhStage;
                bitHeapReqBits = reqBitsForRange(rangeBH[stage][0], rangeBH[stage][1]);
                bitHeapUsedBits = reqBitsForRange2Complement(rangeBH[stage][0], rangeBH[stage][1]);

                bitHeapTmp = new BitHeap(this, bitHeapUsedBits-1, 0, bhName.str(), 1);
                oneLL = static_cast<long long>(1);

                // consider bits
                for (int j = 0; j < bits[stage].size(); ++j) {
                    for (int k = 0; k < bits[stage][j].size(); ++k) {
                        if (bits[stage][j][k] != -1) {
                            ostringstream signalName;

                            long long bitOrigin = bits[stage][j][k] - (oneLL << (48));

                            if (bitOrigin - (oneLL << (32)) >= 0) {
                                bitOrigin = bitOrigin - (oneLL << (32));
                                signalName << "B_" << bitOrigin << "_" << (bhStage);
                                REPORT(INFO,"bitheap --- " << bitHeapTmp->getName() << " signal use: " << signalName.str() << " bit value: " << bitOrigin << " weight: " << j);

                                // for use with bits as std_logic signals
//                                ostringstream signalNameNeg;
//                                signalNameNeg << "B_" << bitOrigin << "_" << bhStage << "N";
//                                vhdl << tab << declare(
//                                        signalNameNeg.str(), 1, true) << tab << "<= \"1\" when " << signalName.str() << " = '1' else \"0\";" << endl;
                                bitHeapTmp->subtractSignal(signalName.str(), j);
                            } else {
                                signalName << "B_" << bitOrigin << "_" << (bhStage);
                                REPORT(INFO,"bitheap +++ " << bitHeapTmp->getName() << " signal use: " << signalName.str() << " bit value: " << bitOrigin << " weight: " << j);
                                bitHeapTmp->addSignal(signalName.str(), j);
                            }
                        }
                    }
                }
                bitHeapTmp->startCompression();

                rangeBH[stage+1] = rangeBH[stage];
                REPORT(DEBUG, "rangeBH in next stage is: " << rangeBH[stage+1][1]);
                if (-modulo <= rangeBH[stage+1][0] && rangeBH[stage+1][1] < modulo) {
                    REPORT(INFO, "modulo rangeBH reached in bitheap: break");
                    break;
                }

                for (int i = 0; i < bitHeapReqBits; ++i) {
                    bits[stage+1][i][0] = i;
                }

                for (int i = 0; i < bitHeapReqBits; ++i) {
                    ostringstream bitName;
                    bitName << "B_" << i << "_" << (bhStage+1);
                    vhdl << tab << declare(
                            bitName.str(), 1, true) << tab << "<=" << bitHeapTmp->getSumName() << range(i,i) << ";" << endl;
                }
                //delete bitHeapTmp;
                bhStage++;
            }
            stage++;
        }
//        vhdl << tab << declare(
//                "T", bitHeapReqBits, false) << tab << "<= STD_LOGIC_VECTOR(SIGNED(" << bitHeapTmp->getSumName(bitHeapReqBits-1,0) << ") + M);" << endl;

        vhdl << tab << declare(
                "STemp", bitHeapUsedBits, false) << tab << "<= STD_LOGIC_VECTOR(SIGNED(" << bitHeapTmp->getSumName(bitHeapUsedBits-1,0) << ") + M);" << endl;


        vhdl << tab << "S <= " << bitHeapTmp->getSumName(outPutSize-1,0) << " when ";
        vhdl << bitHeapTmp->getSumName() << of(bitHeapUsedBits-1) << " = '0'";
        vhdl << tab << "else STemp" << range(outPutSize-1,0) << ";" << endl;
        addFullComment("End of vhdl generation");
        delete bitHeapTmp;
    }

    void ModuloBitHeapOperator::applyPseudoCompressors() {
        REPORT(INFO,"applyPseudoCompressors complete search");
        int smsb = getMSBInStage();
        vector<vector<long long>> rangeContrib(smsb+1,vector<long long>(5));
        long long oneLL = static_cast<long long>(1);

        // fills the first two rows with the modulus and the congruent negative modulus
        for (int i = 0; i <= smsb; ++i) {
            rangeContrib[i][0] = (oneLL << (i)) % modulo;
            rangeContrib[i][1] = rangeContrib[i][0] - modulo;
            if (i == smsb && rangeBH[stage][0] < 0) {
                rangeContrib[i][0] = (((-oneLL << (i)) % modulo) + modulo) % modulo;
                rangeContrib[i][1] = rangeContrib[i][0] - modulo;
            }
        }

        int maxBhHeight = INT_MAX;
        int minBitsToBitHeap = INT_MAX;
        int minBHRange = INT_MAX;
        int rangeReached = 0;
        vector<vector<long long>> ranges(2,vector<long long>(2));
        vector<vector<vector<long long>>> tmpBitHeap(smsb+1,vector<vector<long long> >(64,vector<long long>(2)));

        // initialize tmpBitHeap with -1
        for (int i = 0; i < tmpBitHeap.size(); ++i) {
            for (int j = 0; j < tmpBitHeap[i].size(); ++j) {
                for (int k = 0; k < tmpBitHeap[i][j].size(); ++k) {
                    tmpBitHeap[i][j][k] = -1;
                }
            }
        }

        for (int combi = 0; combi <= (oneLL << (smsb+1))-1; ++combi) {
            long long combiMin = 0;
            long long combiMax = 0;
            int bhHeight = 0;
            int bitsToBitHeap = 0;
            // calculate rangeBH max and min
            for (int i = 0; i <= smsb ; ++i) {
                if ((combi & (oneLL << (i))) != 0 ) {
                    rangeContrib[i][2] = rangeContrib[i][0];
                    combiMax = combiMax + rangeContrib[i][0];
                } else {
                    rangeContrib[i][2] = rangeContrib[i][1];
                    combiMin = combiMin + rangeContrib[i][1];
                }
            }

            vector<vector<long long>> combiBitHeap(reqBitsForRange(combiMin, combiMax),vector<long long>(64));
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
                    if ((rangeContrib[i][2] & (oneLL << (j))) != 0 ) {
                        int bit = 0;
                        while (combiBitHeap[j][bit] != -1) {
                            bit++;
                        }
                        bitsToBitHeap++;
                        bhHeights[j] = bhHeights[j] + 1;
                        if (bhHeight < bhHeights[j]) {
                            bhHeight = bhHeights[j];
                        }
                        if (j < leadingZeroWeight+1 || 0 <= rangeContrib[i][2]) {
                            combiBitHeap[j][bit] = i + (oneLL << (48));
                        } else {
                            combiBitHeap[j][bit] = i + (oneLL << (48)) + (oneLL << (32));
                            constVec = (constVec + (oneLL << (combiBitHeap.size())) - 1 -
                                    ((oneLL << (leadingZeroWeight+1)) - 1));
                            break;
                        }
                    }
                }
            }

            for (int i = 0; i < combiBitHeap.size(); ++i) {
                if ((constVec & (oneLL << (i))) != 0) {
                    bitsToBitHeap++;
                    bhHeights[i] = bhHeights[i] + 1;
                    if (bhHeight < bhHeights[i]) {
                        bhHeight = bhHeights[i];
                    }
                }
            }
            // final stage of pseudocompression
            if (-modulo <= combiMin && combiMax < modulo) {
                REPORT(INFO, "rangeReached");
                rangeReached = 1;
                ranges[1] = {combiMin, combiMax};
                for (int i = 0; i < rangeContrib.size(); ++i) {
                    rangeContrib[i][4] = rangeContrib[i][2];
                }
                tmpBitHeap.resize(combiBitHeap.size());
                for (int i = 0; i < tmpBitHeap.size(); ++i) {
                    tmpBitHeap[i].resize(64);
                    for (int j = 0; j < tmpBitHeap[i].size(); ++j) {
                        tmpBitHeap[i][j].resize(2);
                    }
                }
                for (int i = 0; i < combiBitHeap.size(); ++i) {
                    for (int j = 0; j < combiBitHeap[i].size(); ++j) {
                        tmpBitHeap[i][j][1] = combiBitHeap[i][j];
                    }
                }
                break;
            }
            // fill bits in temporary bitheap
            if (bitsToBitHeap < minBitsToBitHeap && bhHeight <= maxBhHeight
                && (combiMax - combiMin) <= minBHRange) {
                minBitsToBitHeap = bitsToBitHeap;
                minBHRange = combiMax - combiMin;
                ranges[0] = {combiMin, combiMax};
                for (int i = 0; i < rangeContrib.size(); ++i) {
                    rangeContrib[i][3] = rangeContrib[i][2];
                }
                maxBhHeight = bhHeight;
                tmpBitHeap.resize(combiBitHeap.size());
                for (int i = 0; i < tmpBitHeap.size(); ++i) {
                    tmpBitHeap[i].resize(64);
                    for (int j = 0; j < tmpBitHeap[i].size(); ++j) {
                        tmpBitHeap[i][j].resize(2);
                    }
                }
                for (int i = 0; i < combiBitHeap.size(); ++i) {
                    for (int j = 0; j < combiBitHeap[i].size(); ++j) {
                        tmpBitHeap[i][j][0] = combiBitHeap[i][j];
                    }
                }
            }
        }
        // fill bits from tmp bitheap in bits
        rangeBH[stage+1] = ranges[rangeReached];
        int reqCols = reqBitsForRange(ranges[rangeReached][0], ranges[rangeReached][1]);
        for (int i = 0; i < reqCols; ++i) {
            for (int j = 0; j < tmpBitHeap[i].size(); ++j) {
                bits[stage+1][i][j] = tmpBitHeap[i][j][rangeReached];
            }
        }
    }

    void ModuloBitHeapOperator::applyPseudoCompressorsMinBits() {
        REPORT(INFO,"applyPseudoCompressors min bits");
        int smsb = getMSBInStage();
        vector<vector<long long>> rangeContrib(smsb+1,vector<long long>(5));
        long long oneLL = static_cast<long long>(1);
        long long combiMin = 0;
        long long combiMax = 0;

        // fills the first two rows with the modulus and the congruent negative modulus
        for (int i = 0; i <= smsb; ++i) {
            rangeContrib[i][0] = (oneLL << (i)) % modulo;
            rangeContrib[i][1] = rangeContrib[i][0] - modulo;
            if (i == smsb && rangeBH[stage][0] < 0) {
                rangeContrib[i][0] = (((-oneLL << (i)) % modulo) + modulo) % modulo;
                rangeContrib[i][1] = rangeContrib[i][0] - modulo;
            }

            if (countOnes(rangeContrib[i][0]) < countOnes(rangeContrib[i][1])) {
                rangeContrib[i][2] = rangeContrib[i][0];
            } else {
                rangeContrib[i][2] = rangeContrib[i][1];
            }

            if (rangeContrib[i][2] >= 0) {
                combiMax += rangeContrib[i][2];
            } else {
                combiMin += rangeContrib[i][2];
            }
        }
        int bhHeight = 0;
        int bitsToBitHeap = 0;

        vector<vector<long long>> combiBitHeap(reqBitsForRange(combiMin, combiMax),vector<long long>(64));
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
                if ((rangeContrib[i][2] & (oneLL << (j))) != 0 ) {
                    int bit = 0;
                    while (combiBitHeap[j][bit] != -1) {
                        bit++;
                    }
                    bitsToBitHeap++;
                    bhHeights[j] = bhHeights[j] + 1;
                    if (bhHeight < bhHeights[j]) {
                        bhHeight = bhHeights[j];
                    }
                    if (j < leadingZeroWeight+1 || 0 <= rangeContrib[i][2]) {
                        combiBitHeap[j][bit] = i + (oneLL << (48));
                    } else {
                        combiBitHeap[j][bit] = i + (oneLL << (48)) + (oneLL << (32));
                        constVec = (constVec + (oneLL << (combiBitHeap.size())) - 1 -
                                    ((oneLL << (leadingZeroWeight+1)) - 1));
                        break;
                    }
                }
            }
        }
        for (int i = 0; i < combiBitHeap.size(); ++i) {
            if ((constVec & (oneLL << (i))) != 0) {
                bitsToBitHeap++;
                bhHeights[i] = bhHeights[i] + 1;
                if (bhHeight < bhHeights[i]) {
                    bhHeight = bhHeights[i];
                }
            }
        }

        // fill bits from tmp bitheap in bits
        rangeBH[stage+1] = {combiMin, combiMax};
        int reqCols = reqBitsForRange(combiMin, combiMax);

        for (int i = 0; i < reqCols; ++i) {
            for (int j = 0; j < combiBitHeap[i].size(); ++j) {
                bits[stage+1][i][j] = combiBitHeap[i][j];
            }
        }
    }

    void ModuloBitHeapOperator::applyPseudoCompressorsMinRange() {
        REPORT(INFO,"applyPseudoCompressors minimize range");
        int smsb = getMSBInStage();
        vector<vector<long long>> rangeContrib(smsb+1,vector<long long>(5));
        long long oneLL = static_cast<long long>(1);
        long long combiMin = 0;
        long long combiMax = 0;

        // fills the first two rows with the modulus and the congruent negative modulus
        for (int i = 0; i <= smsb; ++i) {
            rangeContrib[i][0] = (oneLL << (i)) % modulo;
            rangeContrib[i][1] = rangeContrib[i][0] - modulo;
            if (i == smsb && rangeBH[stage][0] < 0) {
                rangeContrib[i][0] = (((-oneLL << (i)) % modulo) + modulo) % modulo;
                rangeContrib[i][1] = rangeContrib[i][0] - modulo;
            }

            // min range condition
            if (abs(rangeContrib[i][0]) < abs(rangeContrib[i][1])) {
                rangeContrib[i][2] = rangeContrib[i][0];
            } else {
                rangeContrib[i][2] = rangeContrib[i][1];
            }

            if (rangeContrib[i][2] >= 0) {
                combiMax += rangeContrib[i][2];
            } else {
                combiMin += rangeContrib[i][2];
            }
        }
        int bhHeight = 0;
        int bitsToBitHeap = 0;

        vector<vector<long long>> combiBitHeap(reqBitsForRange(combiMin, combiMax),vector<long long>(64));
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
                if ((rangeContrib[i][2] & (oneLL << (j))) != 0 ) {
                    int bit = 0;
                    while (combiBitHeap[j][bit] != -1) {
                        bit++;
                    }
                    bitsToBitHeap++;
                    bhHeights[j] = bhHeights[j] + 1;
                    if (bhHeight < bhHeights[j]) {
                        bhHeight = bhHeights[j];
                    }
                    if (j < leadingZeroWeight+1 || 0 <= rangeContrib[i][2]) {
                        combiBitHeap[j][bit] = i + (oneLL << (48));
                    } else {
                        combiBitHeap[j][bit] = i + (oneLL << (48)) + (oneLL << (32));
                        constVec = (constVec + (oneLL << (combiBitHeap.size())) - 1 -
                                    ((oneLL << (leadingZeroWeight+1)) - 1));
                        break;
                    }
                }
            }
        }
        for (int i = 0; i < combiBitHeap.size(); ++i) {
            if ((constVec & (oneLL << (i))) != 0) {
                bitsToBitHeap++;
                bhHeights[i] = bhHeights[i] + 1;
                if (bhHeight < bhHeights[i]) {
                    bhHeight = bhHeights[i];
                }
            }
        }

        // fill bits from tmp bitheap in bits
        rangeBH[stage+1] = {combiMin, combiMax};
        int reqCols = reqBitsForRange(combiMin, combiMax);

        for (int i = 0; i < reqCols; ++i) {
            for (int j = 0; j < combiBitHeap[i].size(); ++j) {
                bits[stage+1][i][j] = combiBitHeap[i][j];
            }
        }
    }

    void ModuloBitHeapOperator::applyPseudoCompressorsMinRangeWeighted() {
        REPORT(INFO,"applyPseudoCompressors minimize range weighted bits");
        int smsb = getMSBInStage();
        vector<vector<long long>> rangeContrib(smsb+1,vector<long long>(5));
        long long oneLL = static_cast<long long>(1);
        long long combiMin = 0;
        long long combiMax = 0;

        // fills the first two rows with the modulus and the congruent negative modulus
        for (int i = 0; i <= smsb; ++i) {
            rangeContrib[i][0] = (oneLL << (i)) % modulo;
            rangeContrib[i][1] = rangeContrib[i][0] - modulo;
            if (i == smsb && rangeBH[stage][0] < 0) {
                rangeContrib[i][0] = (((-oneLL << (i)) % modulo) + modulo) % modulo;
                rangeContrib[i][1] = rangeContrib[i][0] - modulo;
            }

            // min range + weight condition
            if ((countOnes(rangeContrib[i][0]) + abs(rangeContrib[i][0]))
                < (countOnes(rangeContrib[i][1]) + abs(rangeContrib[i][1]))) {
                rangeContrib[i][2] = rangeContrib[i][0];
            } else {
                rangeContrib[i][2] = rangeContrib[i][1];
            }

            if (rangeContrib[i][2] >= 0) {
                combiMax += rangeContrib[i][2];
            } else {
                combiMin += rangeContrib[i][2];
            }
        }
        int bhHeight = 0;
        int bitsToBitHeap = 0;

        vector<vector<long long>> combiBitHeap(reqBitsForRange(combiMin, combiMax),vector<long long>(64));
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
                if ((rangeContrib[i][2] & (oneLL << (j))) != 0 ) {
                    int bit = 0;
                    while (combiBitHeap[j][bit] != -1) {
                        bit++;
                    }
                    bitsToBitHeap++;
                    bhHeights[j] = bhHeights[j] + 1;
                    if (bhHeight < bhHeights[j]) {
                        bhHeight = bhHeights[j];
                    }
                    if (j < leadingZeroWeight+1 || 0 <= rangeContrib[i][2]) {
                        combiBitHeap[j][bit] = i + (oneLL << (48));
                    } else {
                        combiBitHeap[j][bit] = i + (oneLL << (48)) + (oneLL << (32));
                        constVec = (constVec + (oneLL << (combiBitHeap.size())) - 1 -
                                    ((oneLL << (leadingZeroWeight+1)) - 1));
                        break;
                    }
                }
            }
        }
        for (int i = 0; i < combiBitHeap.size(); ++i) {
            if ((constVec & (oneLL << (i))) != 0) {
                bitsToBitHeap++;
                bhHeights[i] = bhHeights[i] + 1;
                if (bhHeight < bhHeights[i]) {
                    bhHeight = bhHeights[i];
                }
            }
        }

        // fill bits from tmp bitheap in bits
        rangeBH[stage+1] = {combiMin, combiMax};
        int reqCols = reqBitsForRange(combiMin, combiMax);

        for (int i = 0; i < reqCols; ++i) {
            for (int j = 0; j < combiBitHeap[i].size(); ++j) {
                bits[stage+1][i][j] = combiBitHeap[i][j];
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
        long long oneLL = static_cast<long long>(1);
        int mmsb = 0;
        for (int w = 0; w < 31; ++w) {
            if ((modulo & (oneLL << (w))) != 0) {
                mmsb = w;
            }
        }
        return mmsb;
    }

    int ModuloBitHeapOperator::reqBitsForRange(int min, int max) {
        REPORT(DEBUG, "reqBitsForRange min: " << min << " max " << max);
        int bit = 0;
        if (max > 0) {
            while((max >> bit) != 0) {
                bit++;
            }
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

    int ModuloBitHeapOperator::reqBitsForRange2Complement(int min, int max) {
        REPORT(DEBUG, "reqBitsForRange 2 complement min: " << min << " max " << max);
        int bit = 0;
        if (max > 0) {
            while((max >> bit) != 0) {
                bit++;
            }
        }
        if (min < 0) {
            int negBit = 0;
            int currentMin = min;
            while(currentMin <= -1) {
                currentMin = currentMin / 2;
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

    int ModuloBitHeapOperator::getLeadingZero(long long value) {
        long long oneLL = static_cast<long long>(1);
        int msz = -1;
        for (int w = 0; w < 32; ++w) {
            if ((value & (oneLL << (w))) == 0) {
                msz = w;
            }
        }
        return msz;
    }

    int ModuloBitHeapOperator::countOnes(long long value) {
        long long oneLL = static_cast<long long>(1);
        int bits = 0;
        int limit = 63;
        if (value < 0) {
            limit = getLeadingZero(value);
        }

        for (int i = 0; i <= limit; ++i) {
            if ((value & (oneLL << i)) == (oneLL << i)) {
                bits++;
            }
        }
        return bits;
    }

    bool ModuloBitHeapOperator::checkForRepetition() {
        bool same = false;
        for (int s = 0; s <= stage; ++s) {
            int c = 0;
            same = true;
            while (c < bits[stage].size() && bits[s][c][0] != -1) {
                if (bits[stage+1][c][0] != bits[s][c][0]) {
                    same = false;
                }
                c++;
            }
            if (same) {
                return same;
            }
        }
        return same;
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

        if (wIn >= 10) {
            tc = new TestCase(this);
            tc->addInput("X", mpz_class(634));
            emulate(tc);
            tcl->add(tc);
        }
        if (wIn >= 8) {
            tc = new TestCase(this);
            tc->addInput("X", mpz_class(222));
            emulate(tc);
            tcl->add(tc);
        }
        if (wIn >= 7) {
            tc = new TestCase(this);
            tc->addInput("X", mpz_class(120));
            emulate(tc);
            tcl->add(tc);
        }
        if (wIn >= 6) {
            tc = new TestCase(this);
            tc->addInput("X", mpz_class(33));
            emulate(tc);
            tcl->add(tc);

            tc = new TestCase(this);
            tc->addInput("X", mpz_class(40));
            emulate(tc);
            tcl->add(tc);
        }
        if (wIn >= 4) {
            tc = new TestCase(this);
            tc->addInput("X", mpz_class(15));
            emulate(tc);
            tcl->add(tc);
        }
    }

    TestList ModuloBitHeapOperator::unitTest(int index)
    {
        // the static list of mandatory tests
        TestList testStateList;
        vector<pair<string,string>> paramList;
        string methods[] = {"compl","lMinBits","minRange","minRangeWeight"};

        if(index==-1) 	{ // The unit tests
            for (int i = 0; i < 4; ++i) {
                for (int wIn=6; wIn<=12; wIn++) {
                    for(int m=2; m<20; m++) {
                        paramList.push_back(make_pair("wIn",to_string(wIn)));
                        paramList.push_back(make_pair("modulo",to_string(m)));
                        paramList.push_back(make_pair("method",methods[i]));
                        paramList.push_back(make_pair("TestBench n=",to_string(100)));
                        testStateList.push_back(paramList);

                        paramList.clear();
                    }
                }
            }
        }
        return testStateList;
    }

    OperatorPtr ModuloBitHeapOperator::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args) {
        int wIn, modulo;
        string method;
        UserInterface::parseInt(args, "wIn", &wIn);
        UserInterface::parseInt(args, "modulo", &modulo);
        UserInterface::parseString(args, "method", &method);
        return new ModuloBitHeapOperator(parentOp, target, wIn, modulo, method);
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
                            modulo(int): modulo; \
                            method(string)=compl: The method used for the pseudo compression."
                           " 4 available: 'compl' - for the complete search, 'lMinBits' - for the local bit minimization"
                           " 'minRange' - for the per bit range minimization, 'minRangeWeight' - for a combination of local bit and range minimization",
                // More documentation for the HTML pages. If you want to link to your blog, it is here.
                           "See the developer manual in the doc/ directory of FloPoCo.",
                           ModuloBitHeapOperator::parseArguments,
                           ModuloBitHeapOperator::unitTest
        ) ;
    }
}//namespace
