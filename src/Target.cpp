/*
  The abstract class that models different target FGPAs for FloPoCo. 
  Classes for real chips (in the Targets directory) inherit from this one.
 
  Author : Florent de Dinechin
 
  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  2008-2011.
  All rights reserved.

 */

#include <iostream>
#include <sstream>
#include <string>
#include "Target.hpp"
#include <math.h>


using namespace std;


namespace flopoco{

	extern int verbose;

	string Target::getID(){
		return id_;
	}

	void Target::setPipelined() {
		pipeline_=true;
	}

	void Target::setNotPipelined() {
		pipeline_=false;
	}

	bool Target::isPipelined() {
		return pipeline_;
	}

	int Target::lutInputs() {
		return lutInputs_;
	}

	double Target::frequency(){
		return frequency_;
	}


	double Target::frequencyMHz(){
		return floor(frequency_/1000000);
	}

	double Target::normalizedFrequency(){
		return frequencyMHz()/maxFrequencyMHz_;
	}
	void Target::setFrequency(double f){
		frequency_ = f;
	}

	void Target::setUseHardMultipliers(bool v){
		useHardMultipliers_ = v;  
	}

	bool Target::getUseHardMultipliers(){
		return useHardMultipliers_ ;
	} 

	/*----------------- Resource Estimation related items ------------*/
	
	//TODO: ensure that this function is implemented in the classes that
	//		extend Target, else it will cause failures
	//		for now it defaults to the Target default lutInputs_
	int Target::suggestLUTType(){
		
		return lutInputs();
	}
	
	//TODO: get a better approximation of the number of LUT used
	//		currently, the LUT can at most be split into 2 indep. functions
	int Target::suggestLUTCount(int count, int type){
		int lutType, increment;
		
		(type == 0) ? lutType = suggestLUTType() : lutType = type;
		if(lutSplitInputs(lutType)){
			return count/2;
		}else
			return count;
	}
	
	
	//FIXME: for now, just simply return the number of multipliers as a
	//result when the widths are smaller than the standard multiplier size
	int Target::suggestMultiplierCount(int count, int widthX, int widthY){
		double ratioX, ratioY;
		int sizeDSPx, sizeDSPy;
		
		getDSPWidths(&sizeDSPx, &sizeDSPy, true);
		ratioX = (double)widthX/sizeDSPx;
		ratioY = (double)widthY/sizeDSPy;
		return count*ceil(ratioX)*ceil(ratioY);
	}
	
	//TODO: create a data structure/ function that gives the multiplier 
	//		input width levels and the associated number of DSPs
	int Target::suggestDSPfromMultiplier(int count, int widthX, int widthY){
		double ratioX, ratioY;
		int sizeDSPx, sizeDSPy;
		
		if(dspSplitMult()){
			return count*getMultPerDSP(widthX, widthY);
		}else{
			getDSPWidths(&sizeDSPx, &sizeDSPy, true);
			ratioX = (double)widthX/sizeDSPx;
			ratioY = (double)widthY/sizeDSPy;
			return count*ceil(ratioX)*ceil(ratioY);
		}
	}
		
	//TODO: create a data structure/ function that gives the possible 
	//		word sizes 
	//TODO: take into account the memory type (RAM or ROM); depending on 
	//		the type, might be implemented through distributed memory
	int Target::suggestMemoryCount(int count, int size, int width, int type){
		int wordsPerBlock;
		
		wordsPerBlock = getWordsPerBlock(width);
		return ceil(count*(double)size/wordsPerBlock);
	}
	
	
	int Target::suggestSRLCount(int count, int width, int depth){
		int defaultShifterWidth;
		
		defaultShifterDepth = getSRLDepth(depth);
		return count*width*ceil((double)depth/defaultShifterDepth);
	}
	
	
	int Target::suggestLUTfromSRL(int count, int width, int depth){
		double srlPerLut;
		
		srlPerLut = getLUTPerSRL(depth);
		if(srlPerLut != 0){
			return ceil(count*width*srlPerLut);
		}else
			return 0;
	}
	
	
	int Target::suggestFFfromSRL(int count, int width, int depth){
		int ffPerSRL;
		
		ffPerSRL = getFFPerSRL(depth);
		if(ffPerSRL != 0){
			return ceil(count*width*ffPerSRL);
		}else
			return 0;
	}
	
	
	int Target::suggestRAMfromSRL(int count, int width, int depth){
		int ramPerSRL;
		
		ramPerSRL = getRAMPerSRL(width, depth);
		if(ramPerSRL != 0){
			return ceil(count*ramPerSRL);
		}else
			return 0;
	}
	
	//TODO: get a more accurate count of the number of multiplexers 
	//		needed; currently specific resources are not taken into account
	int Target::suggestMuxCount(int count, int nrInputs, int width){
		
		return count;
	}
	
	
	int Target::suggestLUTfromMux(int count, int nrInputs, int width){
		int stdInputs;
		
		stdInputs = 1;
		while(stdInputs<nrInputs)
			stdInputs *= 2;
		if(stdInputs == 1)
			return 0;
		return ceil(count*width*getLUTfromMux(nrInputs));
	}
	
	
	int Target::suggestLUTfromCounter(int count, int width){
		
		return ceil(count*width*getLUTPerCounter(width));
	}
	
	
	int Target::suggestFFfromCounter(int count, int width){
		
		return ceil(count*width*getFFPerCounter(width));
	}
	
	
	int Target::suggestLUTfromAccumulator(int count, int width, bool useDSP){
		
		return ceil(count*width*getLUTPerAccumulator(width, useDSP));
	}
	
	
	int Target::suggestFFfromAccumulator(int count, int width, bool useDSP){
		
		return ceil(count*width*getFFPerAccumulator(width, useDSP));
	}
	
	
	int Target::suggestDSPfromAccumulator(int count, int width, bool useDSP){
		
		if(useDSP)
			return ceil(count*width*getDSPPerAccumulator(width));
		else
			return 0;
	}
	
	
	int Target::suggestLUTfromDecoder(int count, int width){
		
		return ceil(count*width*getLutFromDecoder(width));
	}
	
	
	int Target::suggestFFfromDecoder(int count, int width){
		
		return ceil(count*width*getLutFromDecoder(width));
	}
	
	
	int Target::suggestRAMfromDecoder(int count, int width){
		
		return ceil(count*width*getRAMFromDecoder(width));
	}
	
	
	int Target::suggestLUTfromArithmeticOperator(int count, int nrInputs, int width){
		double lutsPerArithOp;
		
		lutsPerArithOp = (double)nrInputs/lutInputs();
		return ceil(count*width*lutsPerArithOp);
	}
	
	
	//TODO: find a better approximation for the resources
	//		currently just logic corresponding to the multiplexers
	int Target::suggestLUTfromFSM(int count, int nrStates, int nrTransitions){
		int lutCount = 0;
		
		lutCount += suggestLUTfromMux(count*nrStates*nrTransitions, 2, ceil(log2(nrStates)));
		return lutCount;
	}
	
	
	//TODO: find a better approximation for the resources
	//		currently just logic corresponding to the state register
	int Target::suggestFFfromFSM(int count, int nrStates, int nrTransitions){
		int ffCount = 0;
		
		ffCount += ceil(count*log2(nrStates));
		return ffCount;
	}
	
	
	//TODO: find a better approximation for the resources
	//		for now, RAM blocks are not used
	int Target::suggestRAMfromFSM(int count, int nrStates, int nrTransitions){
		
		return 0;
	}
	
	/*-------- Resource Estimation - target specific functions -------*/
	
	bool lutSplitInputs(int lutType){
		bool fpgaDefault;
		
		if(vendor_ == "Xilinx"){
			if(id_ == "Spartan3"){
				fpgaDefault = false;
			}else if(id_ == "Virtex4"){
				fpgaDefault = false;
			}else if(id_ == "Virtex5"){
				fpgaDefault = true;
			}else if(id_ == "Virtex6"){
				fpgaDefault = true;
			}
		}else if(vendor_.compare("Altera") == 0){
			if(id_ == "StratixII"){
				fpgaDefault = true;
			}else if(id_ == "StratixIII"){
				fpgaDefault = true;
			}else if(id_ == "StratixIV"){
				fpgaDefault = true;
			}
		}
		return fpgaDefault && (lutType <= ceil(lutInputs()/2.0));
	}
	
	
	bool dspSplitMult(){
		
		if(vendor_ == "Xilinx"){
			if(id_ == "Spartan3"){
				return false;
			}else if(id_ == "Virtex4"){
				return false;
			}else if(id_ == "Virtex5"){
				return false;
			}else if(id_ == "Virtex6"){
				return false;
			}
		}else if(vendor_.compare("Altera") == 0){
			if(id_ == "StratixII"){
				return true;
			}else if(id_ == "StratixIII"){
				return true;
			}else if(id_ == "StratixIV"){
				return true;
			}
		}
	}
	
	
	double getMultPerDSP(int widthX, int widthY){
		int defaultWidthX, defaultWidthY, maxWidth;
		
		getDSPWidths(&defaultWidthX, &defaultWidthY, true);
		(defaultWidthX>defaultWidthY) ? maxWidth = defaultWidthX : maxWidth = defaultWidthY;
		if(vendor_ == "Xilinx"){
			if(id_ == "Spartan3"){
				return ceil((double)widthX/defaultWidthX)*ceil((double)widthY/defaultWidthY);
			}else if(id_ == "Virtex4"){
				return ceil((double)widthX/defaultWidthX)*ceil((double)widthY/defaultWidthY);
			}else if(id_ == "Virtex5"){
				return ceil((double)widthX/defaultWidthX)*ceil((double)widthY/defaultWidthY);
			}else if(id_ == "Virtex6"){
				return ceil((double)widthX/defaultWidthX)*ceil((double)widthY/defaultWidthY);
			}
		}else if(vendor_.compare("Altera") == 0){
			if(id_ == "StratixII"){
				if(maxWidth<=9){
					return 1.0/8;
				}else if(maxWidth<=18){
					return 1.0/4;
				}else if(maxWidth<=36){
					return 1;
				}else{
					return ceil((double)widthX/defaultWidthX)*ceil((double)widthY/defaultWidthY);
				}
			}else if(id_ == "StratixIII"){
				if(maxWidth<=9){
					return 1.0/8;
				}else if(maxWidth<=12){
					return 1.0/6;
				}else if(maxWidth<=18){
					return 1.0/4;
				}else if(maxWidth<=36){
					return 1.0/2;
				}else{
					return ceil((double)widthX/defaultWidthX)*ceil((double)widthY/defaultWidthY);
				}
			}else if(id_ == "StratixIV"){
				if(maxWidth<=9){
					return 1.0/8;
				}else if(maxWidth<=12){
					return 1.0/6;
				}else if(maxWidth<=18){
					return 1.0/4;
				}else if(maxWidth<=36){
					return 1.0/2;
				}else{
					return ceil((double)widthX/defaultWidthX)*ceil((double)widthY/defaultWidthY);
				}
			}
		}
	}
	
	//TODO: give indication of the memory type being used, so as to get 
	//		better results (for example, Altera has several types of  
	//		RAM memory, and depending on the type, the size varies)
	//		for now, the function returns the best results, mixing 
	//		memory types
	int wordsPerBlock(int width){
		
		if(vendor_ == "Xilinx"){
			if(id_ == "Spartan3"){
				if(width<=1){
					return 16384;
				}else if(width<=2){
					return 8192;
				}else if(width<=4){
					return 4096;
				}else if(width<=8){
					return 2048;
				}else if(width<=16){
					return 1024;
				}else if(width<=32){
					return 512;
				}else if(width<=72){
					return 256;
				}
			}else if(id_ == "Virtex4"){
				if(width<=1){
					return 16384;
				}else if(width<=2){
					return 8192;
				}else if(width<=4){
					return 4096;
				}else if(width<=9){
					return 2048;
				}else if(width<=18){
					return 1024;
				}else if(width<=36){
					return 512;
				}
			}else if(id_ == "Virtex5"){
				if(width<=1){
					return 32768;
				}else if(width<=2){
					return 16384;
				}else if(width<=4){
					return 8192;
				}else if(width<=9){
					return 4096;
				}else if(width<=18){
					return 2048;
				}else if(width<=36){
					return 1024;
				}else if(width<=72){
					return 512;
				}
			}else if(id_ == "Virtex6"){
				if(width<=1){
					return 32768;
				}else if(width<=2){
					return 16384;
				}else if(width<=4){
					return 8192;
				}else if(width<=9){
					return 4096;
				}else if(width<=18){
					return 2048;
				}else if(width<=36){
					return 1024;
				}else if(width<=72){
					return 512;
				}
			}
		}else if(vendor_.compare("Altera") == 0){
			if(id_ == "StratixII"){
				if(width<=1){
					return 4096;
				}else if(width<=2){
					return 2048;
				}else if(width<=4){
					return 1024;
				}else if(width<=8){
					return 65536;
				}else if(width<=9){
					return 65536;
				}else if(width<=16){
					return 32768;
				}else if(width<=18){
					return 32768;
				}else if(width<=32){
					return 16384;
				}else if(width<=64){
					return 8192;
				}else if(width<=72){
					return 8192;
				}else if(width<=128){
					return 4096;
				}else if(width<=144){
					return 4096;
				}
			}else if(id_ == "StratixIII"){
				if(width<=1){
					return 8192;
				}else if(width<=2){
					return 4096;
				}else if(width<=4){
					return 2048;
				}else if(width<=8){
					return 16384;
				}else if(width<=9){
					return 16384;
				}else if(width<=16){
					return 8192;
				}else if(width<=18){
					return 8192;
				}else if(width<=32){
					return 4096;
				}else if(width<=36){
					return 4096;
				}else if(width<=64){
					return 2048;
				}else if(width<=72){
					return 2048;
				}
			}else if(id_ == "StratixIV"){
				if(width<=1){
					return 8192;
				}else if(width<=2){
					return 4096;
				}else if(width<=4){
					return 2048;
				}else if(width<=8){
					return 16384;
				}else if(width<=9){
					return 16384;
				}else if(width<=16){
					return 8192;
				}else if(width<=18){
					return 8192;
				}else if(width<=32){
					return 4096;
				}else if(width<=36){
					return 4096;
				}else if(width<=64){
					return 2048;
				}else if(width<=72){
					return 2048;
				}
			}
		}
	}
	
	//TODO: check validity of the shift register widths for the Altera
	//		families of devices
	int getSRLDepth(int depth){
		
		if(vendor_ == "Xilinx"){
			if(id_ == "Spartan3"){
				return ceil(depth/16.0);
			}else if(id_ == "Virtex4"){
				return ceil(depth/16.0);
			}else if(id_ == "Virtex5"){
				if(depth<=16)
					return ceil(depth/16.0);
				else
					return ceil(depth/32.0);
			}else if(id_ == "Virtex6"){
				if(depth<=16)
					return ceil(depth/16.0);
				else
					return ceil(depth/32.0);
			}
		}else if(vendor_.compare("Altera") == 0){
			if(id_ == "StratixII"){
				if(depth<=18)
					return ceil(depth/16.0);
				else if(depth<=36)
					return ceil(depth/36.0);
				else if(depth<=144)
					return ceil(depth/144.0);
			}else if(id_ == "StratixIII"){
				if(depth<=18)
					return ceil(depth/16.0);
				else if(depth<=36)
					return ceil(depth/36.0);
				else if(depth<=72)
					return ceil(depth/72.0);
			}else if(id_ == "StratixIV"){
				if(depth<=18)
					return ceil(depth/16.0);
				else if(depth<=36)
					return ceil(depth/36.0);
				else if(depth<=72)
					return ceil(depth/72.0);
			}
		}
	}
	
	
	double getLUTPerSRL(int depth){
		cerr << "Error: function getLUTPerSRL() not yet implemented" << endl;
		exit(1);
	}
	
	
	double getFFPerSRL(int depth){
		cerr << "Error: function getFFPerSRL() not yet implemented" << endl;
		exit(1);
	}
	
	
	double getRAMPerSRL(int depth){
		cerr << "Error: function getRAMPerSRL() not yet implemented" << endl;
		exit(1);
	}
	
	
	double getLUTFromMux(int nrInputs){
		cerr << "Error: function getLUTFromMux() not yet implemented" << endl;
		exit(1);
	}
	
	
	double getLUTPerCounter(int width){
		cerr << "Error: function getLUTPerCounter() not yet implemented" << endl;
		exit(1);
	}
	
	
	double getFFPerCounter(int width){
		cerr << "Error: function getFFPerCounter() not yet implemented" << endl;
		exit(1);
	}
	
	
	double getLUTPerAccumulator(int width, bool useDSP){
		cerr << "Error: function getLUTPerAccumulator() not yet implemented" << endl;
		exit(1);
	}
	
	
	double getFFPerAccumulator(int width, bool useDSP){
		cerr << "Error: function getFFPerAccumulator() not yet implemented" << endl;
		exit(1);
	}
	
	
	double getDSPPerAccumulator(int width){
		cerr << "Error: function getDSPPerAccumulator() not yet implemented" << endl;
		exit(1);
	}
	
	
	double getLUTPerDecoder(int width){
		cerr << "Error: function getLUTPerDecoder() not yet implemented" << endl;
		exit(1);
	}
	
	
	double getFFPerDecoder(int width){
		cerr << "Error: function getFFPerDecoder() not yet implemented" << endl;
		exit(1);
	}
	
	
	double getRAMPerDecoder(int width){
		cerr << "Error: function getRAMPerDecoder() not yet implemented" << endl;
		exit(1);
	}
	
	/*----------------------------------------------------------------*/

}















