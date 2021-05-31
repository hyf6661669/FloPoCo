
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
		bool adderReached = false;
//        int moduloRangeMax = 0;
//        int moduloRangeMin = 0;

		vector<int> currentRanges(bitheap->width);
		negativeSignExtension = 0;
		needToApplyNegativeSignExtension = false;
		int currentPossibleMaxRange = 0;
		int extraRangeToSubtract = 0;
		string compressionMode = bitheap->mode;

        long long oneLL = static_cast<long long>(1);
        for (int i = 0; i < bitheap->width; ++i) {
            //moduloRangeMax += (oneLL << i);
            currentRanges[i] = (oneLL << i);
            currentPossibleMaxRange += (oneLL << i);
        }

        if (compressionMode.find("msbcases") != string::npos) {
            if (bitheap->maxInput != -1) {
                currentPossibleMaxRange = bitheap->maxInput;
            }
        }

		while(!adderReached){
            bool reachedModuloRange = true;

            int moduloRangeMax = 0;
            int moduloRangeMin = 0;
            int remainderExtension = 0;


            if (computeModulo) {
                remainderExtension = negativeSignExtension % bitheap->modulus;

                for (int i = 0; i < currentRanges.size(); ++i) {
                    cerr << "current ranges " << i << " is " << currentRanges[i] << endl;
                    if (currentRanges[i] >= 0) {
                        moduloRangeMax += currentRanges[i];
                    } else {
                        moduloRangeMin += currentRanges[i];
                    }
                }

                if (compressionMode.find("msbcases") != string::npos) {
                    moduloRangeMax = getMaxRangeForMaxValue(currentPossibleMaxRange, currentRanges);
                }

                if (compressionMode.find("singlebit") != string::npos) {
                    moduloRangeMax -= extraRangeToSubtract;
                }

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
            unsigned int requiredBitsForRange = reqBitsForRange2Complement(moduloRangeMin, moduloRangeMax);


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


                    for(unsigned int c = 0; c < bitAmount[s].size(); c++){
                        int rangeChange = placePseudoCompressor(s, c, requiredBitsForRange, true, useNegativeMSBValue);
                        cerr << "set range change " << rangeChange << " for column " << c << endl;
                        if (c < currentRanges.size()) {
                            currentRanges[c] = rangeChange;
                        }
                    }
                    if (compressionMode.find("msbcases") != string::npos) {
                        currentPossibleMaxRange = moduloRangeMax;
                    }
                } else {
                    cerr << "normal compression" << endl;

                    // place PseudoCompressors where column height = 1
                    // TODO: leave reachedModuloRange here?
                    if (compressionMode.find("singlebit") != string::npos && !reachedModuloRange) {
                        cerr << "placing single bit pseudo comps" << endl;
                        bool useNegativeMSBValue = false;

                        if (compressionMode.find("sevector") == string::npos) {
                            if (moduloRangeMin < 0) {
                                useNegativeMSBValue = true;
                            }
                        }
                        for(unsigned int c = 0; c < bitAmount[s].size(); c++){
                            if (bitAmount[s][c] == 1) {
                                int rangeChange = placePseudoCompressor(s, c, requiredBitsForRange, false,  useNegativeMSBValue);
                                cerr << "set range change 2 " << rangeChange << " for column " << c << endl;
                                int bitWeight = 1 << c;
                                int newRangeChange = bitWeight - rangeChange;
                                cerr << "newRangeChange " << newRangeChange << endl;
                                extraRangeToSubtract += newRangeChange;
                                cerr << "extraRangeToSubtract " << extraRangeToSubtract << endl;
//                                while (newRangeChange > 0) {
//                                    for (int i = 0; i < currentRanges.size(); i++) {
//                                        if (newRangeChange > 0 && currentRanges[i] > 0) {
//                                            currentRanges[i] -= 1;
//                                            newRangeChange--;
//                                        }
//                                    }
//                                }
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
                            cerr << "placed compressor " << compressor->getStringOfIO() << " in stage " << s << " and column " << column << endl;
                            REPORT(DETAILED, "placed compressor " << compressor->getStringOfIO() << " in stage " << s << " and column " << column);
                            REPORT(DETAILED, "efficiency is " << achievedEfficiencyBest);
                            placeCompressor(s, column, compressor, middleLengthBest);
                        }
                    }
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
			if (s > 30 && computeModulo) {
			    cerr << "break because stage limit reached" << endl;
			    break;
			}
			s++;
		}
	}

    int MaxEfficiencyCompressionStrategy::reqBitsForRange2Complement(int min, int max) {
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
	    } else {
            return abs(rem) < abs(remToCompare);
	    }
	}

    bool MaxEfficiencyCompressionStrategy::shouldPlacePseudoCompressors(vector<int> bitAmountStage) {
        int maxHeightBitAmount = *max_element(bitAmountStage.begin(), bitAmountStage.end());
        return maxHeightBitAmount == 1;
	}

    int MaxEfficiencyCompressionStrategy::placePseudoCompressor(int s, int column, int requiredBitsForRange, bool allowDeletion, bool useNegativeMSBValue) {

        BasicCompressor* compressor = nullptr;
        bool found = false;

        if (useNegativeMSBValue && column == requiredBitsForRange - 1 && bitAmount[s][column] > 0) {
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

            return compressor->range_change;
        } else {
            return 0;
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

    int MaxEfficiencyCompressionStrategy::getMaxRangeForMaxValue(int maxValue, vector<int> currentRanges) {
	    int currentMaxRange = 0;
	    int currentOneMSBRange = 0;
        int maxInputSize = floor(log2(maxValue)+1);
        for (int i = maxInputSize-1; i >= 0; i--) {
            if ((maxValue & (1 << i)) != 0) {
                // case 0
                int msbZeroRange = currentOneMSBRange;
                for (int j = 0; j < currentRanges.size(); ++j) {
                    if (j < i) {
                        msbZeroRange += currentRanges[j];
                    }
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
}
