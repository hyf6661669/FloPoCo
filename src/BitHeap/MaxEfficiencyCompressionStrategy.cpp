
#include "MaxEfficiencyCompressionStrategy.hpp"
#include "ModReduction/PseudoCompressor.hpp"
#include "ModReduction/RowAdder.hpp"
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
        possibleCompressors.push_back(rowAdder);

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

		int moduloRangeMin = 0;
		int moduloRangeMax = 0;
        long long oneLL = static_cast<long long>(1);
        for (int i = 0; i < bitheap->width; ++i) {
            moduloRangeMax += (oneLL << i);
        }


		while(!adderReached){
            bool breakToWhile = false;
            bool reachedModuloRange = true;

            if (computeModulo) {
                cerr << "modMin: " << moduloRangeMin << " modMax: " << moduloRangeMax << endl;
                reachedModuloRange = (moduloRangeMin >= -bitheap->modulus && moduloRangeMax < bitheap->modulus);
            }


            if(checkAlgorithmReachedAdder(2, s) && reachedModuloRange){
                break;
            }

			//make sure there is the stage s+1 with the same amount of columns as s
			while(bitAmount.size() <= s + 1){
				bitAmount.resize(bitAmount.size() + 1);
				bitAmount[bitAmount.size() - 1].resize(bitAmount[bitAmount.size() - 2].size(), 0);
			}

			// check if mod is used and if stage height = 1
			// then use other compressors -> pseudocompressors
			int maxHeightBitAmount = *max_element(bitAmount[s].begin(), bitAmount[s].end());
			if (maxHeightBitAmount == 1 && computeModulo) {
			    bool useNegativeMSBValue = false;
			    unsigned int requiredBitsForRange = reqBitsForRange2Complement(moduloRangeMin, moduloRangeMax);
			    if (moduloRangeMin < 0) {
			        useNegativeMSBValue = true;
			    }
			    moduloRangeMax = 0;
			    moduloRangeMin = 0;

                bool found = true;
                while(found) {
                    found = false;
                    BasicCompressor* compressor = nullptr;

                    // go through columns
                    // check if allowed
                    // check if another is better
                    // check if there is a bit in this column
                    for(unsigned int c = 0; c < bitAmount[s].size(); c++){
                        found = false;
                        unsigned int currentRange = INT_MAX;

                        if (useNegativeMSBValue && c == requiredBitsForRange - 1 && bitAmount[s][c] > 0) {
                            // make new compressor for negative MSB
                            vector<int> compInput(requiredBitsForRange, 0);
                            compInput[compInput.size()-1] = 1;
                            int modulo = bitheap->modulus;
                            int wIn = bitheap->width;

                            int newRem = (((-1 << c) % modulo) + modulo) % modulo;
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

                                currentRange = abs(newRem);
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

                                currentRange = abs(newReciprocal);
                                compressor = new BasicPseudoCompressor(bitheap->getOp(), bitheap->getOp()->getTarget(), compInput, compOutputRec, newReciprocal, ones_vector_start);
                                found = true;
                            }

                        } else if (c >= requiredBitsForRange && bitAmount[s][c] > 0 && c < bitheap->width) {
                            cerr << "set 0 pseudo comp c is " << c << endl;
                            vector<int> compInput(c+1, 0);
                            compInput[compInput.size()-1] = 1;
                            cerr << "compInput size: " << compInput.size() << endl;
                            vector<int> compOutput(bitAmount[s].size(), 0);
                            compOutput[bitAmount[s] .size()-1] = 1;
                            compressor = new BasicPseudoCompressor(bitheap->getOp(), bitheap->getOp()->getTarget(), compInput, compOutput, 0);
                            found = true;
                        } else {
                            for (unsigned int i = 0; i < possibleCompressors.size(); ++i) {
                                if (possibleCompressors[i]->type == CompressorType::Pseudo) {
                                    if (possibleCompressors[i]->heights.size() - 1 == c && bitAmount[s][c] > 0 && c < requiredBitsForRange) {
                                        if (abs(possibleCompressors[i]->range_change) < currentRange) {
                                            currentRange = possibleCompressors[i]->range_change;
                                            compressor = possibleCompressors[i];
                                            found = true;
                                        }
                                    }
                                }
                            }
                        }

                        if(found){
                            REPORT(DETAILED, "range change is " << compressor->range_change);
                            placeCompressor(s, 0, compressor);
                            if (compressor->range_change >= 0) {
                                moduloRangeMax += compressor->range_change;
                            } else {
                                moduloRangeMin += compressor->range_change;
                            }
                            printBitAmounts();
                        }
                    }
                }
			} else {
			    cerr << "normal compression" << endl;

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
                        REPORT(DETAILED, "efficiency is " << achievedEfficiencyBest);
                        placeCompressor(s, column, compressor, middleLengthBest);
                    }
                }
			}

			if (breakToWhile) {
			    break;
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
}
