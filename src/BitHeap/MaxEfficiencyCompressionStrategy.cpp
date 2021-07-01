
#include "MaxEfficiencyCompressionStrategy.hpp"
#include "ModReduction/PseudoCompressor.hpp"
#include "ModReduction/RowAdder.hpp"
#include "BitHeap/ConstantAddCompressor.hpp"
//#include "CompressionStrategy.hpp"
//#include "BitHeap/BitHeap.hpp"


using namespace std;

namespace flopoco{


	MaxEfficiencyCompressionStrategy::MaxEfficiencyCompressionStrategy(BitHeap* bitheap) : CompressionStrategy(bitheap)
	{
		lowerBounds.resize(1);
		lowerBounds[0] = 0.0;
	}




	void MaxEfficiencyCompressionStrategy::compressionAlgorithm()
	{
		REPORT(DEBUG, "compressionAlgorithm is maxEfficiency");
		// add row adder
        BasicCompressor *rowAdder = new BasicRowAdder(bitheap->getOp(), bitheap->getOp()->getTarget(), 2);
        //possibleCompressors.push_back(rowAdder);

		//for the maxEfficiency algorithm, the compressors should be ordered by efficiency
		orderCompressorsByCompressionEfficiency();

		//adds the Bits to stages and columns
		orderBitsByColumnAndStage();

		//populates bitAmount. on this simple structure the maxEfficiency algorithm is working
		fillBitAmounts();

		//prints out how the inputbits of the bitheap looks like
		printBitAmounts();

		//new solution
		solution = BitHeapSolution();
		solution.setSolutionStatus(BitheapSolutionStatus::HEURISTIC_PARTIAL);

		//generates the compressor tree. Works only on bitAmount, compressors will be put into solution
		maxEfficiencyAlgorithm();

		//reports the area in LUT-equivalents
        printSolutionStatistics();

		//here the VHDL-Code for the compressors as well as the bits->compressors->bits are being written.
		applyAllCompressorsFromSolution();
	}

	void MaxEfficiencyCompressionStrategy::maxEfficiencyAlgorithm(){
		unsigned int s = 0;
		logicalStages = 0;
		bool adderReached = false;
//        int moduloRangeMax = 0;
//        int moduloRangeMin = 0;

		vector<mpz_class> currentRanges(bitheap->width);
		vector<bool> currentRangesInvertedBits(bitheap->width);

		negativeSignExtension = 0;
		needToApplyNegativeSignExtension = false;

		// set to false when single pseudo comps don't improve the range
		bool shouldUseSinglePseudoComps = true;
		bool useSingleBitSubtraction = false;

		// range case differentiation
		mpz_class currentPossibleMaxRange = 0;
		// se vector
		mpz_class extraRangeToSubtract = 0;
		// single bit
		mpz_class currentRangeMaxRange = 0;
		// single bit and msb cases. Max range before remainderExtension is subtracted
		mpz_class currentOnlyPositiveMaxRange = 0;

		string compressionMode = bitheap->mode;

        mpz_class oneMpz = mpz_class(1);
        for (int i = 0; i < bitheap->width; ++i) {
            //moduloRangeMax += (oneLL << i);
            currentRanges[i] = (oneMpz << i);
            currentPossibleMaxRange += (oneMpz << i);
            currentRangeMaxRange += (oneMpz << i);
            currentOnlyPositiveMaxRange += (oneMpz << i);
        }

        if (compressionMode.find("msbcases") != string::npos) {
            if (bitheap->maxInput != -1) {
                currentPossibleMaxRange = bitheap->maxInput;
            }
        }

		while(!adderReached){
            bool reachedModuloRange = true;

            mpz_class moduloRangeMax = 0;
            mpz_class moduloRangeMin = 0;
            int remainderExtension = 0;


            if (computeModulo) {
                cerr << "negativeSignExtension " << negativeSignExtension << endl;
                remainderExtension = negativeSignExtension % bitheap->modulus;

                for (int i = 0; i < currentRanges.size(); ++i) {
                    //cerr << "current ranges " << i << " is " << currentRanges[i] << endl;
                    if (currentRanges[i] >= 0) {
                        moduloRangeMax += currentRanges[i];
                    } else {
                        moduloRangeMin += currentRanges[i];
                    }

                }
                cerr << "current ranges max " << moduloRangeMax << " min " << moduloRangeMin << endl;

                if (compressionMode.find("msbcases") != string::npos) {
                    //moduloRangeMax = getMaxRangeForMaxValue(currentPossibleMaxRange, currentRanges);
                    moduloRangeMax = getMaxRangeForMaxValue(currentPossibleMaxRange, currentRanges, currentRangesInvertedBits);
                    //vector<bool> pseudoCompSet(bitheap->width, true);
                    //moduloRangeMax = getMaxRangeForStage(currentPossibleMaxRange, currentRanges, bitAmount[s], pseudoCompSet, currentRangesInvertedBits);
                }

                if (compressionMode.find("singlebit") != string::npos) {
                    currentRangeMaxRange = moduloRangeMax;
                    if (useSingleBitSubtraction) {
                        moduloRangeMax -= extraRangeToSubtract;
                    }

                    //cerr << "extraRangeToSubtract" << extraRangeToSubtract << endl;
                }

                currentOnlyPositiveMaxRange = moduloRangeMax;

                if (compressionMode.find("sevector") != string::npos) {
                    moduloRangeMax += remainderExtension;
                    moduloRangeMin += remainderExtension;
                }

                cerr << "modMin: " << moduloRangeMin << " modMax: " << moduloRangeMax << endl;
                if (compressionMode.find("pos") != string::npos) {
                    reachedModuloRange = moduloRangeMax < bitheap->modulus*2;
                } else {
                    reachedModuloRange = (moduloRangeMin >= -bitheap->modulus && moduloRangeMax < bitheap->modulus);
                }
            }

            if(checkAlgorithmReachedAdder(2, s) && reachedModuloRange && !needToApplyNegativeSignExtension){
                cerr << "reached Adder and range: break" << endl;
                break;
            }

			//make sure there is the stage s+1 with the same amount of columns as s
			while(bitAmount.size() <= s + 1){
				bitAmount.resize(bitAmount.size() + 1);
				bitAmount[bitAmount.size() - 1].resize(bitAmount[bitAmount.size() - 2].size(), 0);
			}
            //unsigned int requiredBitsForRange = reqBitsForRange2Complement(moduloRangeMin, moduloRangeMax);
            unsigned int requiredBitsForRange;

            if (compressionMode.find("sevector") != string::npos) {
                requiredBitsForRange = reqBitsForRange2Complement(0, currentOnlyPositiveMaxRange);
            } else {
                requiredBitsForRange = reqBitsForRange2Complement(moduloRangeMin, moduloRangeMax);
            }


            if (negativeSignExtension < 0 && shouldPlacePseudoCompressors(bitAmount[s]) && reachedModuloRange) {
                vector<int> compInput(bitheap->width, 1);
                BasicCompressor* compressor = nullptr;
                compressor = new BasicConstantAddCompressor(bitheap->getOp(), bitheap->getOp()->getTarget(), compInput, remainderExtension);
                placeCompressor(s, 0, compressor);

                needToApplyNegativeSignExtension =  false;
            } else {
                // pseudo compression
                if (shouldPlacePseudoCompressors(bitAmount[s]) && computeModulo) {
                    cerr << "pseudo compression" << endl;
                    bool useNegativeMSBValue = false;

                    if (compressionMode.find("sevector") == string::npos) {
                        if (moduloRangeMin < 0) {
                            useNegativeMSBValue = true;
                        }

                        extraRangeToSubtract = 0;
                    }


                    currentRangesInvertedBits.clear();
                    for(unsigned int c = 0; c < bitAmount[s].size(); c++){
                        pair<int,bool> resultPlacedComp = placePseudoCompressor(s, c, requiredBitsForRange, true, useNegativeMSBValue);
                        int rangeChange = resultPlacedComp.first;

                        if (c < currentRanges.size()) {
                            currentRanges[c] = rangeChange;

                            if (resultPlacedComp.second) {
                                currentRangesInvertedBits.push_back(true);
                            } else {
                                currentRangesInvertedBits.push_back(false);
                            }
                        }
                    }

                    if (compressionMode.find("msbcases") != string::npos) {
                        currentPossibleMaxRange = currentOnlyPositiveMaxRange;
                    }
                    if (compressionMode.find("singlebit") != string::npos) {
                        shouldUseSinglePseudoComps = true;
                        useSingleBitSubtraction = false;
                    }

                    mpz_class afterPseudoCompMax = 0;
                    mpz_class afterPseudoCompMin = 0;
                    for (int i = 0; i < currentRanges.size(); ++i) {
                        //cerr << "current ranges " << i << " is " << currentRanges[i] << endl;
                        if (currentRanges[i] >= 0) {
                            afterPseudoCompMax += currentRanges[i];
                        } else {
                            afterPseudoCompMin += currentRanges[i];
                        }
                    }

                    if (afterPseudoCompMax == moduloRangeMax && afterPseudoCompMin == moduloRangeMin) {
                        if (bitheap->pseudoCompMode == "lMinBits") {
                            bitheap->pseudoCompMode = "minRange";
                        }
                    }

                } else {
                    cerr << "normal compression" << endl;

                    // place PseudoCompressors where column height = 1
                    // TODO: leave reachedModuloRange here?
                    if (compressionMode.find("singlebit") != string::npos && !reachedModuloRange && shouldUseSinglePseudoComps) {

                        vector<int> bitDistributionStage = bitAmount[s];
                        int bitsInStage = 0;

                        for (int i = 0; i < bitDistributionStage.size(); ++i) {
                            for (int j = 0; j < bitDistributionStage[i]; ++j) {
                                bitsInStage++;
                            }
                        }

                        // for more bits the range calculation takes very long
                        if (bitsInStage <= 30) {
                            bool useNegativeMSBValue = false;

                            if (compressionMode.find("sevector") == string::npos) {
                                if (moduloRangeMin < 0) {
                                    useNegativeMSBValue = true;
                                }
                            }

                            vector<bool> pseudoCompSet;
                            vector<bool> invertedRangeBits;

                            bool pseudoCompWasSet = false;
                            bool pseudoCompChangedRange = false;

                            for(unsigned int c = 0; c < bitAmount[s].size(); c++){
                                if (bitAmount[s][c] == 1 && c < bitheap->width) {
                                    pair<int,bool> resultPlacedComp = placePseudoCompressor(s, c, requiredBitsForRange, false,  useNegativeMSBValue);
                                    cerr << "placed single bit pseudo comp " << resultPlacedComp.first << " at " << c << endl;

                                    pseudoCompSet.push_back(true);
                                    pseudoCompWasSet = true;
                                    if (resultPlacedComp.second) {
                                        invertedRangeBits.push_back(true);
                                    } else {
                                        invertedRangeBits.push_back(false);
                                    }

                                    if (resultPlacedComp.first != mpz_class(1) << c) {
                                        pseudoCompChangedRange = true;
                                    }
                                } else {
                                    pseudoCompSet.push_back(false);
                                    invertedRangeBits.push_back(false);
                                }
                            }
                            if (pseudoCompWasSet) {
                                mpz_class newRangeForStage = getMaxRangeForStage(currentOnlyPositiveMaxRange, currentRanges, bitDistributionStage, pseudoCompSet, invertedRangeBits);
                                extraRangeToSubtract = currentRangeMaxRange - newRangeForStage;
                                //cerr << "single bit currentRangeMaxRange " << currentRangeMaxRange << endl;
                                //cerr << "single bit newRangeForStage " << newRangeForStage << endl;
                                //cerr << "single bit currentOnlyPositiveMaxRange " << currentOnlyPositiveMaxRange << endl;

                                if (currentOnlyPositiveMaxRange == newRangeForStage && pseudoCompChangedRange) {
                                    //cerr << "shouldUseSinglePseudoComps false" << endl;
                                    shouldUseSinglePseudoComps = false;
                                }
                                useSingleBitSubtraction = true;
                            }
                        }
                    }

                    bool found = true;
                    while(found){
                        found = false;

                        double achievedEfficiencyBest = -1.0;
                        BasicCompressor* compressor = nullptr;
                        unsigned int column = 0;
                        unsigned int middleLengthBest = 0;

                        for(unsigned int e = 0; e < possibleCompressors.size(); e++){
                            if (possibleCompressors[e]->type != CompressorType::Pseudo) {
                                BasicCompressor* currentCompressor = possibleCompressors[e];
                                REPORT(DEBUG, "compressor is " << currentCompressor->getStringOfIO());
                                vector<bool> used;
                                used.resize(bitAmount[s].size(), false);

                                unsigned int columnsAlreadyChecked = 0;
                                //check if the achievedEfficiency is better than the maximal efficiency possible by this compressor. If true, it's not necessary to check this and the following compressors. Therefore return.
                                while(columnsAlreadyChecked < bitAmount[s].size() && !((found == true) && currentCompressor->getEfficiency() - achievedEfficiencyBest < 0.0001)){

                                    unsigned int currentMaxColumn = 0;
                                    int currentSize = 0;
                                    for(unsigned int c = 0; c < bitAmount[s].size(); c++){
                                        if(!used[c] && bitAmount[s][c] > currentSize){
                                            currentMaxColumn = c;
                                            currentSize = bitAmount[s][c];
                                        }
                                    }
                                    used[currentMaxColumn] = true;

                                    double achievedEfficiencyCurrent;
                                    unsigned int middleLengthCurrent = 0;
                                    if(currentCompressor->type != CompressorType::Variable){
                                        achievedEfficiencyCurrent = getCompressionEfficiency(s, currentMaxColumn, currentCompressor);
                                    } else {
                                        double variableEfficiencyBest = -1.0;
                                        for (int m = 0; m < bitAmount[s].size() - 2; m++) {
                                            double currentVariableEfficiency = getCompressionEfficiency(s, currentMaxColumn, currentCompressor, m);
                                            if (currentVariableEfficiency >= variableEfficiencyBest) {
                                                variableEfficiencyBest = currentVariableEfficiency;
                                                middleLengthCurrent = m;
                                            }
                                        }
                                        achievedEfficiencyCurrent = variableEfficiencyBest;
                                    }

                                REPORT(FULL, "checked " << currentCompressor->getStringOfIO() << " in stage " << s << " and column " << currentMaxColumn << " with an efficiency of " << achievedEfficiencyCurrent);

                                    float lowerBound;
                                    if(s < lowerBounds.size())
                                        lowerBound = lowerBounds[s];
                                    else
                                        lowerBound = 0.0;

                                    if(achievedEfficiencyCurrent > (achievedEfficiencyBest + 0.0001) && achievedEfficiencyCurrent > (lowerBound - 0.0001)){
                                        achievedEfficiencyBest = achievedEfficiencyCurrent;
                                        compressor = currentCompressor;
                                        found = true;
                                        column = currentMaxColumn;
                                        middleLengthBest = middleLengthCurrent;
                                    }
                                    columnsAlreadyChecked++;
                                }
                            }
                        }
                        if(found){

                            REPORT(DETAILED, "placed compressor " << compressor->getStringOfIO() << " in stage " << s << " and column " << column);
                            cerr << "placed compressor " << compressor->getStringOfIO() << " in stage " << s << " and column " << column << endl;
                            REPORT(DETAILED, "efficiency is " << achievedEfficiencyBest);
                            placeCompressor(s, column, compressor, middleLengthBest);
                        }
                    }
                    logicalStages++;
                }
            }



			//finished one stage. bring the remaining bits in bitAmount to the new stage
			for(unsigned int c = 0; c < bitAmount[s].size(); c++){
				if(bitAmount[s][c] > 0){
					bitAmount[s + 1][c] += bitAmount[s][c];
					bitAmount[s][c] = 0;
				}
				solution.setEmptyInputsByRemainingBits(s, bitAmount[s]);
			}
			REPORT(DEBUG, "finished stage " << s);
			printBitAmounts();
			if (s > 100 && computeModulo) {
			    cerr << "break because stage limit reached" << endl;
			    break;
			}
			s++;
		}
	}

    int MaxEfficiencyCompressionStrategy::reqBitsForRange2Complement(mpz_class min, mpz_class max) {
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

    bool MaxEfficiencyCompressionStrategy::isRemainderMoreEfficient(int rem, int remToCompare) {
	    if (bitheap->mode.find("pos") != string::npos) {
            return rem >= 0;
	    } else if (bitheap->pseudoCompMode == "minRange") {
            return abs(rem) < abs(remToCompare);
	    } else if (bitheap->pseudoCompMode == "lMinBits"){
	        cerr << "For value " << rem << " contains ones " << countOnes(rem) << endl;
            cerr << "For value 2 " << remToCompare << " contains ones " << countOnes(remToCompare) << endl;
            return countOnes(rem) < countOnes(remToCompare);
	    } else {
            THROWERROR("No matching mode to decide which pseudo compressor is better");
	    }
	}

    int MaxEfficiencyCompressionStrategy::getLeadingZero(long long value) {
        long long oneLL = static_cast<long long>(1);
        int msz = -1;
        for (int w = 0; w < 32; ++w) {
            if ((value & (oneLL << (w))) == 0) {
                msz = w;
            }
        }
        return msz;
    }

    int MaxEfficiencyCompressionStrategy::countOnes(long long value) {
        long long oneLL = static_cast<long long>(1);
        int bits = 0;
        int limit = 63;
        if (value < 0 && bitheap->mode.find("sevector") != string::npos) {
            limit = getLeadingZero(value) + 1;
        }
        //cerr << "countOnes value " << value << " limit is " << limit << endl;

        for (int i = 0; i <= limit; ++i) {
            if ((value & (oneLL << i)) == (oneLL << i)) {
                bits++;
            }
        }
        return bits;
    }

    bool MaxEfficiencyCompressionStrategy::shouldPlacePseudoCompressors(vector<int> bitAmountStage) {
        int maxHeightBitAmount = *max_element(bitAmountStage.begin(), bitAmountStage.end());
        return maxHeightBitAmount == 1;
	}

    pair<int,bool> MaxEfficiencyCompressionStrategy::placePseudoCompressor(int s, int column, int requiredBitsForRange, bool allowDeletion, bool useNegativeMSBValue) {

        BasicCompressor* compressor = nullptr;
        bool found = false;
        bool invertedBitComp = false;

        if (useNegativeMSBValue && column == requiredBitsForRange - 1 && bitAmount[s][column] > 0) {
            cerr << "negative MSB comp" << endl;
            // make new compressor for negative MSB
            vector<int> compInput(requiredBitsForRange, 0);
            compInput[compInput.size()-1] = 1;
            int modulo = bitheap->modulus;
            int wIn = bitheap->width;

            int newRem = (((-1 << column) % modulo) + modulo) % modulo;
            int newReciprocal = newRem - modulo;

            if (abs(newRem) <= abs(newReciprocal)) {
                vector<int> compOutput;
                for(int j = 1; j < 1<<wIn; j <<= 1){
                    if(j&newRem){
                        compOutput.push_back(1);
                    } else {
                        compOutput.push_back(0);
                    }
                }

                compressor = new BasicPseudoCompressor(bitheap->getOp(), bitheap->getOp()->getTarget(), compInput, compOutput, newRem);
                found = true;
            } else {
                vector<int> compOutputRec;
                int ones_vector_start = 0, cnt = 1;
                for(int j = 1; j < 1<<wIn; j <<= 1){
                    if(j&newReciprocal){
                        compOutputRec.push_back(1);
                    } else {
                        compOutputRec.push_back(0);
                        ones_vector_start = cnt;
                    }
                    cnt++;
                }

                compressor = new BasicPseudoCompressor(bitheap->getOp(), bitheap->getOp()->getTarget(), compInput, compOutputRec, newReciprocal, ones_vector_start);
                found = true;
            }

        } else if (column >= requiredBitsForRange && bitAmount[s][column] > 0 && column < bitheap->width && allowDeletion) {
            cerr << "set 0 pseudo comp c is " << column << endl;
            vector<int> compInput(column+1, 0);
            compInput[compInput.size()-1] = 1;
            vector<int> compOutput(bitAmount[s].size(), 0);
            compOutput[bitAmount[s] .size()-1] = 1;
            compressor = new BasicPseudoCompressor(bitheap->getOp(), bitheap->getOp()->getTarget(), compInput, compOutput, 0);
            found = true;
        } else {
            int currentRange = INT_MAX;
            for (unsigned int i = 0; i < possibleCompressors.size(); ++i) {
                if (possibleCompressors[i]->type == CompressorType::Pseudo) {
                    if (possibleCompressors[i]->heights.size() - 1 == column && bitAmount[s][column] > 0 && column < requiredBitsForRange) {
                        if (isRemainderMoreEfficient(possibleCompressors[i]->range_change, currentRange)) {
                            if (possibleCompressors[i]->range_change < 0) {
                                currentRange = possibleCompressors[i]->range_change;
                                if (bitheap->mode.find("sevector") != string::npos) {
                                    compressor = createCompWithoutSignExtension(possibleCompressors[i]);
                                    invertedBitComp = true;
                                } else {
                                    compressor = possibleCompressors[i];
                                }

                                found = true;
                            } else {
                                currentRange = possibleCompressors[i]->range_change;
                                compressor = possibleCompressors[i];
                                found = true;
                            }
                        }
                    }
                }
            }
        }

        if(found){
            REPORT(DETAILED, "range change is " << compressor->range_change);
            placeCompressor(s, 0, compressor);
            if (bitheap->mode.find("sevector") != string::npos) {
                if (compressor->ones_vector_start < INT32_MAX) {
                    negativeSignExtension -= (1 << compressor->ones_vector_start);
                    if (!needToApplyNegativeSignExtension) {
                        needToApplyNegativeSignExtension = true;
                    }
                }
            }

            return make_pair(compressor->range_change, invertedBitComp);
        } else {
            return make_pair(0, false);
        }
	}

    BasicCompressor* MaxEfficiencyCompressionStrategy::createCompWithoutSignExtension(BasicCompressor* compressor) {
        vector<int> compOutput;
        for(int i = 0; i < compressor->outHeights.size(); i++){
            if (i == compressor->ones_vector_start){
                compOutput.push_back(1);
            } else if (i < compressor->ones_vector_start){
                compOutput.push_back(compressor->outHeights[i]);
            } else {
                compOutput.push_back(0);
            }
        }

        int rangeChange = 1 << compressor->ones_vector_start;
        cerr << "range change for ones start " << compressor->ones_vector_start << " is " << rangeChange << endl;

        BasicPseudoCompressor* pseudoCompressor = nullptr;
        pseudoCompressor = new BasicPseudoCompressor(bitheap->getOp(), bitheap->getOp()->getTarget(), compressor->heights, compOutput, rangeChange, compressor->ones_vector_start);
        pseudoCompressor->setHasExternalSignExtension(true);
        return pseudoCompressor;
	}

    mpz_class MaxEfficiencyCompressionStrategy::getMaxRangeForMaxValue(mpz_class maxValue, vector<mpz_class> currentRanges, vector<bool> currentRangesInvertedBits) {
	    mpz_class currentMaxRange = 0;
	    mpz_class currentOneMSBRange = 0;
        int maxInputSize = intlog2(maxValue);

        for (int i = maxInputSize-1; i >= 0; i--) {
            if ((maxValue & (mpz_class(1)) << i) != 0 || currentRangesInvertedBits[i]) {
                // case 0
                mpz_class msbZeroRange = currentOneMSBRange;
                for (int j = 0; j < currentRanges.size(); ++j) {
                    if (j < i) {
                        msbZeroRange += currentRanges[j];
                    }
                }
                if (currentRangesInvertedBits[i]) {
                    msbZeroRange += currentRanges[i];
                }

                if (msbZeroRange > currentMaxRange) {
                    currentMaxRange = msbZeroRange;
                }
                // case 1
                currentOneMSBRange += currentRanges[i];
            }
        }

        if (currentOneMSBRange > currentMaxRange) {
            currentMaxRange = currentOneMSBRange;
        }

        return currentMaxRange;
	}

    mpz_class MaxEfficiencyCompressionStrategy::getMaxRangeForStage(mpz_class maxValue, vector<mpz_class> currentRanges, vector<int> bitDistribution, vector<bool> setPseudoComps, vector<bool> invertedRangeBits) {
	    vector<RangeEntry> actualRanges;
        mpz_class oneMpz = mpz_class(1);

        for (int i = 0; i < bitDistribution.size(); ++i) {
            for (int j = 0; j < bitDistribution[i]; ++j) {
                RangeEntry rangeEntry;
                rangeEntry.isSet = false;

                if (invertedRangeBits[i]) {
                    rangeEntry.weight = 0;
                } else {
                    rangeEntry.weight = oneMpz << i;
                }

                if (setPseudoComps[i]) {
                    rangeEntry.range = currentRanges[i];
                } else {
                    rangeEntry.range = oneMpz << i;
                }
                actualRanges.push_back(rangeEntry);
            }
        }

        return maxRangeForPosition(actualRanges, actualRanges.size()-1, maxValue);
	}

	mpz_class MaxEfficiencyCompressionStrategy::maxRangeForPosition(vector<RangeEntry> actualRanges, int currentPosition, mpz_class maxValue) {
	    mpz_class maxWeightRange = 0;
        mpz_class minWeightRange = 0;

        for (int i = 0; i < actualRanges.size(); ++i) {
            if (i > currentPosition) {
                maxWeightRange += actualRanges[i].isSet ? actualRanges[i].weight : 0;
                minWeightRange += actualRanges[i].isSet ? actualRanges[i].weight : 0;
            } else {
                maxWeightRange += actualRanges[i].weight;
                minWeightRange += 0;
            }
        }

        // break 1: cannot be higher than maxValue
        if (maxWeightRange <= maxValue) {
            mpz_class rangeNow = 0;
            for (int i = 0; i < actualRanges.size(); ++i) {
                if (i > currentPosition) {
                    rangeNow += actualRanges[i].isSet ? actualRanges[i].range : 0;
                } else {
                    rangeNow += actualRanges[i].range;
                }
            }
            return rangeNow;
        }

        // break 2: maxValue not in range
        if (maxValue > maxWeightRange || maxValue < minWeightRange) {
            return -1;
        }


        mpz_class rangeZero = -1;
        mpz_class rangeOne = -1;
        int newPosition = currentPosition-1;

	    // case 0
        actualRanges[currentPosition].isSet = false;
        rangeZero = maxRangeForPosition(actualRanges, newPosition, maxValue);


	    // case 1
        actualRanges[currentPosition].isSet = true;
        rangeOne = maxRangeForPosition(actualRanges, newPosition, maxValue);

        return max(rangeZero, rangeOne);
	}
}
