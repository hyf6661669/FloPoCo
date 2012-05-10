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
	
	//TODO: get a better approximation of the number of LUT used
	//		currently, the LUT can at most be split into 2 indep. functions
	int Target::lutCount(int nrInputs, int count){
		int lutType;
		
		(nrInputs == 0) ? lutType = lutInputs() : lutType = nrInputs;
		if(lutSplitInputs(lutType)){
			return count/2;
		}else
			return count;
	}
	
	
	//FIXME: for now, just simply return the number of multipliers as a
	//result when the widths are smaller than the standard multiplier size
	int Target::multiplierCount(int widthX, int widthY, int count){
		double ratioX, ratioY;
		int sizeDSPx, sizeDSPy;
		
		getDSPWidths(sizeDSPx, sizeDSPy, true);
		ratioX = (double)widthX/sizeDSPx;
		ratioY = (double)widthY/sizeDSPy;
		return count*ceil(ratioX)*ceil(ratioY);
	}
	
	
	//TODO: create a data structure/ function that gives the multiplier 
	//		input width levels and the associated number of DSPs
	int Target::dspForMultiplier(int widthX, int widthY, int count){
		double ratioX, ratioY;
		int sizeDSPx, sizeDSPy;
		
		if(dspSplitMult()){
			return count*getMultPerDSP(widthX, widthY);
		}else{
			getDSPWidths(sizeDSPx, sizeDSPy, true);
			ratioX = (double)widthX/sizeDSPx;
			ratioY = (double)widthY/sizeDSPy;
			return count*ceil(ratioX)*ceil(ratioY);
		}
	}
	
	int Target::lutForMultiplier(int widthX, int widthY, int count){
		int maxWidth;
		
		(widthX>widthY) ? maxWidth=widthX : maxWidth=widthY;
		return ceil(count*maxWidth*getLUTPerMultiplier(widthX, widthY));
	}
	
	
	int Target::ffForMultiplier(int widthX, int widthY, int count){
		int maxWidth;
		
		(widthX>widthY) ? maxWidth=widthX : maxWidth=widthY;
		return ceil(count*maxWidth*getFFPerMultiplier(widthX, widthY));
	}
	
	int Target::lutForAdderSubtracter(int widthX, int widthY, int count){
		int maxWidth;
		
		(widthX>widthY) ? maxWidth=widthX : maxWidth=widthY;
		
		//method based on FloPoCo formulas
		double period;
		int alpha, beta, kappa;
		
		/*
		period = 1.0/frequency_;
		//alpha = 1 + floor((period-2.0*lutDelay())/carryPropagateDelay());
		suggestSubaddSize(alpha, maxWidth);
		kappa = maxWidth/alpha + 1;
		beta = maxWidth - (kappa - 1)*alpha;
		*/
		
		period = 1.0/frequency_;
		suggestSubaddSize(alpha, maxWidth);
		beta = (maxWidth % alpha == 0 ? alpha : maxWidth % alpha);
		kappa = (maxWidth % alpha == 0 ? maxWidth/alpha : maxWidth/alpha + 1);
		
		//TODO: decide formula to use based on the type of the adder
		if(kappa == 1){
			return maxWidth;
		}else if(kappa == 2){
			return alpha+2*beta;
		}else{
			return (4*kappa-8)*alpha+3*beta+kappa-3;
		}
		//
		
		
		//method based on approximations
		//return ceil(count*maxWidth*getLUTPerAdderSubtracter(widthX, widthY));
	}
	
	
	int Target::ffForAdderSubtracter(int widthX, int widthY, int count){
		int maxWidth;
		
		(widthX>widthY) ? maxWidth=widthX : maxWidth=widthY;
		
		//method based on FloPoCo formulas
		double period;
		int alpha, beta, kappa;
		
		/*
		period = 1.0/frequency_;
		alpha = 1 + floor((period-2*lutDelay())/carryPropagateDelay());
		kappa = maxWidth/alpha + 1;
		beta = maxWidth - (kappa - 1)*alpha;
		*/
		
		period = 1.0/frequency_;
		suggestSubaddSize(alpha, maxWidth);
		beta = (maxWidth % alpha == 0 ? alpha : maxWidth % alpha);
		kappa = (maxWidth % alpha == 0 ? maxWidth/alpha : maxWidth/alpha + 1);
		
		//TODO: decide formula to use based on the type of the adder
		if(kappa == 1)
			return maxWidth;
		else
			return (2*kappa-3)*alpha+beta+2*kappa-3;
		//
		
		//method based on approximations
		//return ceil(count*maxWidth*getFFPerAdderSubtracter(widthX, widthY));
	}
	
		
	//TODO: create a data structure/ function that gives the possible 
	//		word sizes 
	//TODO: take into account the memory type (RAM or ROM); depending on 
	//		the type, might be implemented through distributed memory or
	//		dedicated memory blocks
	int Target::memoryCount(int size, int width, int type, int count){
		int WpB;
		
		WpB = wordsPerBlock(width);
		return ceil(count*(double)size/WpB);
	}
	
	
	int Target::srlCount(int width, int depth, int count){
		int defaultShifterDepth;
		
		defaultShifterDepth = getSRLDepth(depth);
		return count*width*ceil((double)depth/defaultShifterDepth);
	}
	
	
	int Target::lutForSRL(int width, int depth, int count){
		double srlPerLut;
		
		srlPerLut = getLUTPerSRL(depth);
		if(srlPerLut != 0){
			return ceil(count*width*srlPerLut);
		}else
			return 0;
	}
	
	
	int Target::ffForSRL(int width, int depth, int count){
		int ffPerSRL;
		
		ffPerSRL = getFFPerSRL(depth);
		if(ffPerSRL != 0){
			return ceil(count*width*ffPerSRL);
		}else
			return 0;
	}
	
	
	int Target::ramForSRL(int width, int depth, int count){
		int ramPerSRL;
		
		ramPerSRL = getRAMPerSRL(depth);
		if(ramPerSRL != 0){
			return ceil(count*ramPerSRL);
		}else
			return 0;
	}
	
	//TODO: get a more accurate count of the number of multiplexers 
	//		needed; currently specific resources are not taken into account
	int Target::muxCount(int nrInputs, int width, int count){
		
		return count;
	}
	
	
	int Target::lutForMux(int nrInputs, int width, int count){
		int stdInputs;
		
		stdInputs = 1;
		while(stdInputs<nrInputs)
			stdInputs *= 2;
			
		if(stdInputs == 1){
			cout << "Warning: lutForMux(): trying to add multiplexer with a single input" << endl;
			return 0;
		}
		return ceil(count*width*getLUTFromMux(stdInputs));
	}
	
	//TODO: get estimations when using specific resources (like DSPs)
	//		involves also changes to getLUTPerCounter()
	int Target::lutForCounter(int width, int count){
		
		return ceil(count*width*getLUTPerCounter(width));
	}
	
	//TODO: get estimations when using specific resources (like DSPs)
	//		involves also changes to getLUTPerCounter()
	int Target::ffForCounter(int width, int count){
		
		return ceil(count*width*getFFPerCounter(width));
	}
	
	
	int Target::lutForAccumulator(int width, bool useDSP, int count){
		
		return ceil(count*width*getLUTPerAccumulator(width, useDSP));
	}
	
	
	int Target::ffForAccumulator(int width, bool useDSP, int count){
		
		return ceil(count*width*getFFPerAccumulator(width, useDSP));
	}
	
	
	int Target::dspForAccumulator(int width, bool useDSP, int count){
		
		if(useDSP)
			return ceil(count*width*getDSPPerAccumulator(width));
		else
			return 0;
	}
	
	
	int Target::lutForDecoder(int width, int count){
		
		return ceil(count*width*getLUTPerDecoder(width));
	}
	
	
	int Target::ffForDecoder(int width, int count){
		
		return ceil(count*width*getFFPerDecoder(width));
	}
	
	
	int Target::lutForArithmeticOperator(int nrInputs, int width, int count){
		double lutsPerArithOp;
		
		lutsPerArithOp = (double)nrInputs/lutInputs();
		return ceil(count*width*lutsPerArithOp);
	}
	
	
	//TODO: find a better approximation for the resources
	//		currently just logic corresponding to the multiplexers
	int Target::lutForFSM(int nrStates, int nrTransitions, int count){
		int lutCount = 0;
		
		lutCount += lutForMux(count*nrStates*nrTransitions, 2, ceil(log2(nrStates)));
		return lutCount;
	}
	
	
	//TODO: find a better approximation for the resources
	//		currently just logic corresponding to the state register
	int Target::ffForFSM(int nrStates, int nrTransitions, int count){
		int ffCount = 0;
		
		ffCount += ceil(count*log2(nrStates));
		return ffCount;
	}
	
	
	//TODO: find a better approximation for the resources
	//		for now, RAM blocks are not used
	int Target::ramForFSM(int nrStates, int nrTransitions, int count){
		
		return 0;
	}
	
	/*-------- Resource Estimation - target specific functions -------*/
	
	bool Target::lutSplitInputs(int nrInputs){
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
		
		return fpgaDefault && (nrInputs <= ceil(lutInputs()/2.0));
	}
	
	
	bool Target::dspSplitMult(){
		
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
		
		return false;
	}
	
	
	double Target::getMultPerDSP(int widthX, int widthY){
		int defaultWidthX, defaultWidthY, maxWidth;
		
		getDSPWidths(defaultWidthX, defaultWidthY, true);
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
		
		return 0;
	}
	
	double Target::getLUTPerMultiplier(int widthX, int widthY){
		int maxWidth;
		
		(widthX>widthY) ? maxWidth=widthX : maxWidth=widthY;
		if(vendor_ == "Xilinx"){
			if(id_ == "Spartan3"){
				if(maxWidth<=9){
					return 12.555; 
				}else if(maxWidth<=12){
					return 14.917;
				}else if(maxWidth<=20){
					if(widthX==18 && widthY==18){
						return 0.000;
					}else{
						return 3.600;
					}
				}else if(maxWidth<=25){
					if((widthX==25 && widthY==18) || (widthY==25 && widthX==18)){
						return 0.680;
					}else{
						return 14.958;
					}
				}else if(maxWidth<=35){
					return 2.486;
				}else if(maxWidth<=53){
					return 9.283;
				}
			}else if(id_ == "Virtex4"){
				if(maxWidth<=9){
					return 12.555; 
				}else if(maxWidth<=12){
					return 14.917;
				}else if(maxWidth<=20){
					if(widthX==18 && widthY==18){
						return 0.000;
					}else{
						return 3.600;
					}
				}else if(maxWidth<=25){
					if((widthX==25 && widthY==18) || (widthY==25 && widthX==18)){
						return 0.680;
					}else{
						return 14.958;
					}
				}else if(maxWidth<=35){
					return 2.486;
				}else if(maxWidth<=53){
					return 9.283;
				}
			}else if(id_ == "Virtex5"){
				if(maxWidth<=9){
					return 12.111; 
				}else if(maxWidth<=12){
					return 13.917;
				}else if(maxWidth<=20){
					if(widthX==18 && widthY==18){
						return 0.000;
					}else{
						return 1.200;
					}
				}else if(maxWidth<=25){
					if((widthX==25 && widthY==18) || (widthY==25 && widthX==18)){
						return 0.000;
					}else{
						return 7.083;
					}
				}else if(maxWidth<=35){
					return 1.000;
				}else if(maxWidth<=53){
					return 4.320;
				}
			}else if(id_ == "Virtex6"){
				if(maxWidth<=9){
					return 12.778;
				}else if(maxWidth<=12){
					return 14.667;
				}else if(maxWidth<=20){
					if(widthX==18 && widthY==18){
						return 0.000;
					}else{
						return 1.150;
					}
				}else if(maxWidth<=25){
					if((widthX==25 && widthY==18) || (widthY==25 && widthX==18)){
						return 0.000;
					}else{
						return 8.333;
					}
				}else if(maxWidth<=35){
					return 0.942;
				}else if(maxWidth<=53){
					return 2.509;
				}
			}
		}else if(vendor_.compare("Altera") == 0){
			if(id_ == "StratixII"){
				//can use just DSPs and no additional logic
				if((widthX==8 && widthY==8) || (widthX==16 && widthY==16) || (widthX==32 && widthY==32))
					return 0.000;
				else
					return 1.281;
			}else if(id_ == "StratixIII"){
				if((widthX==8 && widthY==8) || (widthX==16 && widthY==16) || (widthX==32 && widthY==32))
					return 0.000;
				else
					return 1.281;
			}else if(id_ == "StratixIV"){
				if((widthX==8 && widthY==8) || (widthX==16 && widthY==16) || (widthX==32 && widthY==32))
					return 0.000;
				else
					return 1.281;
			}
		}
		
		return 0;
	}
	
	double Target::getFFPerMultiplier(int widthX, int widthY){
		int maxWidth;
		
		(widthX>widthY) ? maxWidth=widthX : maxWidth=widthY;
		if(vendor_ == "Xilinx"){
			if(id_ == "Spartan3"){
				if(maxWidth<=9){
					return 11.555;
				}else if(maxWidth<=12){
					return 13.750;
				}else if(maxWidth<=20){
					if(widthX==18 && widthY==18){
						return 0.000;
					}else{
						return 3.750;
					}
				}else if(maxWidth<=25){
					if((widthX==25 && widthY==18) || (widthY==25 && widthX==18)){
						return 1.000;
					}else{
						return 3.000;
					}
				}else if(maxWidth<=35){
					return 5.000;
				}else if(maxWidth<=53){
					return 8.566;
				}
			}else if(id_ == "Virtex4"){
				if(maxWidth<=9){
					return 11.555;
				}else if(maxWidth<=12){
					return 13.750;
				}else if(maxWidth<=20){
					if(widthX==18 && widthY==18){
						return 0.000;
					}else{
						return 3.750;
					}
				}else if(maxWidth<=25){
					if((widthX==25 && widthY==18) || (widthY==25 && widthX==18)){
						return 1.000;
					}else{
						return 3.000;
					}
				}else if(maxWidth<=35){
					return 5.000;
				}else if(maxWidth<=53){
					return 8.566;
				}
			}else if(id_ == "Virtex5"){
				if(maxWidth<=9){
					return 12.222;
				}else if(maxWidth<=12){
					return 14.917;
				}else if(maxWidth<=20){
					if(widthX==18 && widthY==18){
						return 0.000;
					}else{
						return 1.200;
					}
				}else if(maxWidth<=25){
					if((widthX==25 && widthY==18) || (widthY==25 && widthX==18)){
						return 0.000;
					}else{
						return 9.833;
					}
				}else if(maxWidth<=35){
					return 2.485;
				}else if(maxWidth<=53){
					return 5.283;
				}
			}else if(id_ == "Virtex6"){
				if(maxWidth<=9){
					return 12.222;
				}else if(maxWidth<=12){
					return 14.917;
				}else if(maxWidth<=20){
					if(widthX==18 && widthY==18){
						return 0.000;
					}else{
						return 1.200;
					}
				}else if(maxWidth<=25){
					if((widthX==25 && widthY==18) || (widthY==25 && widthX==18)){
						return 0.000;
					}else{
						return 9.833;
					}
				}else if(maxWidth<=35){
					return 1.971;
				}else if(maxWidth<=53){
					return 5.207;
				}
			}
		}else if(vendor_.compare("Altera") == 0){
			if(id_ == "StratixII"){
				if((widthX==8 && widthY==8) || (widthX==16 && widthY==16) || (widthX==32 && widthY==32))
					return 0.000;
				else
					return 2.000;
			}else if(id_ == "StratixIII"){
				if((widthX==8 && widthY==8) || (widthX==16 && widthY==16) || (widthX==32 && widthY==32))
					return 0.000;
				else
					return 2.000;
			}else if(id_ == "StratixIV"){
				if((widthX==8 && widthY==8) || (widthX==16 && widthY==16) || (widthX==32 && widthY==32))
					return 0.000;
				else
					return 2.000;
			}
		}
		
		return 0;
	}
	
	
	//TODO: get statistics for Virtex6 (for now using Virtex5 data)
	//		verify statistics for Altera families		
	double Target::getLUTPerAdderSubtracter(int widthX, int widthY){
		int maxWidth;
		
		(widthX>widthY) ? maxWidth=widthX : maxWidth=widthY;
		if(vendor_ == "Xilinx"){
			if(id_ == "Spartan3"){
				if(isPipelined()){
					if(maxWidth<=8){
						return 3.875;
					}else if(maxWidth<=12){
						return 4.500;
					}else if(maxWidth<=16){
						return 4.750;
					}else if(maxWidth<=20){
						return 4.900;
					}else if(maxWidth<=24){
						return 4.917;
					}else if(maxWidth<=32){
						return 5.156;
					}else if(maxWidth<=48){
						return 7.167;
					}else if(maxWidth<=64){
						return 7.968;
					}else if(maxWidth<=96){
						return 9.833;
					}else if(maxWidth<=128){
						return 11.671;
					}else if(maxWidth<=256){
						return 18.179;
					}
				}else{
					return 1.000;
				}
			}else if(id_ == "Virtex4"){
				if(isPipelined()){
					if(maxWidth<=32){
						return 1.000;
					}else if(maxWidth<=48){
						return 1.312;
					}else if(maxWidth<=52){
						return 1.365;
					}else if(maxWidth<=64){
						return 1.484;
					}else if(maxWidth<=96){
						return 2.667;
					}else if(maxWidth<=128){
						return 2.757;
					}else if(maxWidth<=256){
						return 3.347;
					}
				}else{
					return 1.000;
				}
			}else if(id_ == "Virtex5"){
				if(isPipelined()){
					if(maxWidth<=32){
						return 1.000;
					}else if(maxWidth<=48){
						return 1.145;
					}else if(maxWidth<=52){
						return 1.211;
					}else if(maxWidth<=64){
						return 1.328;
					}else if(maxWidth<=96){
						return 2.562;
					}else if(maxWidth<=128){
						return 2.687;
					}else if(maxWidth<=256){
						return 3.214;
					}
				}else{
					return 1.000;
				}
			}else if(id_ == "Virtex6"){
				if(isPipelined()){
					if(maxWidth<=32){
						return 1.000;
					}else if(maxWidth<=48){
						return 1.145;
					}else if(maxWidth<=52){
						return 1.211;
					}else if(maxWidth<=64){
						return 1.328;
					}else if(maxWidth<=96){
						return 2.562;
					}else if(maxWidth<=128){
						return 2.687;
					}else if(maxWidth<=256){
						return 3.214;
					}
				}else{
					return 1.000;
				}
			}
		}else if(vendor_.compare("Altera") == 0){
			if(id_ == "StratixII"){
				if(isPipelined()){
					if(maxWidth<=24){
						return 1.208;
					}else if(maxWidth<=32){
						return 1.406;
					}else if(maxWidth<=48){
						return 1.645;
					}else if(maxWidth<=64){
						return 1.765;
					}else if(maxWidth<=96){
						return 2.062;
					}else if(maxWidth<=128){
						return 2.031;
					}else if(maxWidth<=256){
						return 2.164;
					}
				}else{
					return 1.000;
				}
			}else if(id_ == "StratixIII"){
				if(isPipelined()){
					if(maxWidth<=48){
						return 1.167;
					}else if(maxWidth<=64){
						return 1.375;
					}else if(maxWidth<=96){
						return 1.604;
					}else if(maxWidth<=128){
						return 1.718;
					}else if(maxWidth<=256){
						return 1.910;
					}
				}else{
					return 1.000;
				}
			}else if(id_ == "StratixIV"){
				if(isPipelined()){
					if(maxWidth<=96){
						return 1.187;
					}else if(maxWidth<=128){
						return 1.390;
					}else if(maxWidth<=256){
						return 1.710;
					}
				}else{
					return 1.000;
				}
			}
		}
		
		return 0;
	}
	
	
	//TODO: get statistics for Virtex6 (for now Virtex5 data)
	//		verify statistics for Altera families
	double Target::getFFPerAdderSubtracter(int widthX, int widthY){
		int maxWidth;
		
		(widthX>widthY) ? maxWidth=widthX : maxWidth=widthY;
		if(vendor_ == "Xilinx"){
			if(id_ == "Spartan3"){
				if(isPipelined()){
					if(maxWidth<=8){
						return 1.875;
					}else if(maxWidth<=12){
						return 2.583;
					}else if(maxWidth<=16){
						return 2.687;
					}else if(maxWidth<=20){
						return 2.950;
					}else if(maxWidth<=24){
						return 2.958;
					}else if(maxWidth<=32){
						return 3.750;
					}else if(maxWidth<=48){
						return 4.500;
					}else if(maxWidth<=64){
						return 4.875;
					}else if(maxWidth<=96){
						return 6.437;
					}else if(maxWidth<=128){
						return 6.828;
					}else if(maxWidth<=256){
						return 7.414;
					}
				}else{
					return 0.000;
				}
			}else if(id_ == "Virtex4"){
				if(isPipelined()){
					if(maxWidth<=32){
						return 1.000;
					}else if(maxWidth<=48){
						return 1.020;
					}else if(maxWidth<=52){
						return 1.019;
					}else if(maxWidth<=64){
						return 1.015;
					}else if(maxWidth<=96){
						return 1.010;
					}else if(maxWidth<=128){
						return 1.273;
					}else if(maxWidth<=256){
						return 2.167;
					}
				}else{
					return 0.000;
				}
			}else if(id_ == "Virtex5"){
				if(isPipelined()){
					if(maxWidth<=32){
						return 1.000;
					}else if(maxWidth<=48){
						return 1.020;
					}else if(maxWidth<=52){
						return 1.019;
					}else if(maxWidth<=64){
						return 1.015;
					}else if(maxWidth<=96){
						return 1.010;
					}else if(maxWidth<=128){
						return 1.376;
					}else if(maxWidth<=256){
						return 2.117;
					}
				}else{
					return 0.000;
				}
			}else if(id_ == "Virtex6"){
				if(isPipelined()){
					if(maxWidth<=32){
						return 1.000;
					}else if(maxWidth<=48){
						return 1.020;
					}else if(maxWidth<=52){
						return 1.019;
					}else if(maxWidth<=64){
						return 1.015;
					}else if(maxWidth<=96){
						return 1.010;
					}else if(maxWidth<=128){
						return 1.376;
					}else if(maxWidth<=256){
						return 2.117;
					}
				}else{
					return 0.000;
				}
			}
		}else if(vendor_.compare("Altera") == 0){
			if(id_ == "StratixII"){
				if(isPipelined()){
					if(maxWidth<=24){
						return 0.000;
					}else if(maxWidth<=32){
						return 1.031;
					}else if(maxWidth<=48){
						return 2.062;
					}else if(maxWidth<=64){
						return 3.093;
					}else if(maxWidth<=96){
						return 3.0114;
					}else if(maxWidth<=128){
						return 2.976;
					}else if(maxWidth<=256){
						return 2.593;
					}
				}else{
					return 0.000;
				}
			}else if(id_ == "StratixIII"){
				if(isPipelined()){
					if(maxWidth<=48){
						return 0.000;
					}else if(maxWidth<=64){
						return 1.015;
					}else if(maxWidth<=96){
						return 2.031;
					}else if(maxWidth<=128){
						return 3.046;
					}else if(maxWidth<=256){
						return 5.964;
					}
				}else{
					return 0.000;
				}
			}else if(id_ == "StratixIV"){
				if(isPipelined()){
					if(maxWidth<=96){
						return 0.000;
					}else if(maxWidth<=128){
						return 1.007;
					}else if(maxWidth<=256){
						return 3.023;
					}
				}else{
					return 0.000;
				}
			}
		}
		
		return 0;
	}
	
	//TODO: give indication of the memory type being used, so as to get 
	//		better results (for example, Altera has several types of  
	//		RAM memory, and depending on the type, the size varies)
	//		for now, the function returns the best results, mixing 
	//		memory types
	int Target::wordsPerBlock(int width){
		
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
		
		return 0;
	}
	
	//TODO: check validity of the shift register widths for the Altera
	//		families of devices
	int Target::getSRLDepth(int depth){
		
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
		
		return 0;
	}
	
	
	double Target::getLUTPerSRL(int depth){
		
		if(vendor_ == "Xilinx"){
			if(id_ == "Spartan3"){
				return depth;
			}else if(id_ == "Virtex4"){
				return depth;
			}else if(id_ == "Virtex5"){
				return depth/2;
			}else if(id_ == "Virtex6"){
				return depth/2;
			}
		}else if(vendor_.compare("Altera") == 0){
			if(id_ == "StratixII"){
				return depth/2;
			}else if(id_ == "StratixIII"){
				return depth/2;
			}else if(id_ == "StratixIV"){
				return depth/2;
			}
		}
		
		return 0;
	}
	
	
	double Target::getFFPerSRL(int depth){
		
		if(vendor_ == "Xilinx"){
			if(id_ == "Spartan3"){
				return depth/4;
			}else if(id_ == "Virtex4"){
				return depth/4;
			}else if(id_ == "Virtex5"){
				return depth/4;
			}else if(id_ == "Virtex6"){
				return depth/4;
			}
		}else if(vendor_.compare("Altera") == 0){
			if(id_ == "StratixII"){
				return depth/4;
			}else if(id_ == "StratixIII"){
				return depth/4;
			}else if(id_ == "StratixIV"){
				return depth/4;
			}
		}
		
		return 0;
	}
	
	
	double Target::getRAMPerSRL(int depth){
		
		if(vendor_ == "Xilinx"){
			if(id_ == "Spartan3"){
				return 0;
			}else if(id_ == "Virtex4"){
				return 0;
			}else if(id_ == "Virtex5"){
				return 0;
			}else if(id_ == "Virtex6"){
				return 0;
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
		
		return 0;
	}
	
	
	double Target::getLUTFromMux(int nrInputs){
		
		if(vendor_ == "Xilinx"){
			if(id_ == "Spartan3"){
				if(nrInputs<=32)
					return nrInputs/2;
				else{
					int totalMux = 0, levels = ceil(log2(nrInputs)), i;
					
					i = levels;
					while(i>=0){
						totalMux += nrInputs/32;
						nrInputs /= 32;
						i--;
					}
					return totalMux;
				}
			}else if(id_ == "Virtex4"){
				if(nrInputs==2){
					return nrInputs/4.0;
				}if(nrInputs<=32)
					return ceil(nrInputs/2.0);
				else{
					int totalMux = 0, levels = ceil(log2(nrInputs)), i;
					
					i = levels;
					while(i>=0){
						totalMux += nrInputs/32;
						nrInputs /= 32;
						i--;
					}
					return totalMux;
				}
			}else if(id_ == "Virtex5"){
				if(nrInputs==2){
					return nrInputs/4.0;
				}else if(nrInputs<=16){
					return ceil(nrInputs/4.0);
				}else{
					int totalMux = 0, levels = ceil(log2(nrInputs)), i;
					
					i = levels;
					while(i>=0){
						totalMux += nrInputs/16;
						nrInputs /= 16;
						i--;
					}
					return totalMux;
				}
			}else if(id_ == "Virtex6"){
				if(nrInputs==2){
					return nrInputs/4.0;
				}if(nrInputs<=16)
					return ceil(nrInputs/4.0);
				else{
					int totalMux = 0, levels = ceil(log2(nrInputs)), i;
					
					i = levels;
					while(i>=0){
						totalMux += nrInputs/16;
						nrInputs /= 16;
						i--;
					}
					return totalMux;
				}
			}
		}else if(vendor_.compare("Altera") == 0){
			if(id_ == "StratixII"){
				if(nrInputs<=8)
					return ceil(nrInputs/4.0);
				else{
					int totalMux = 0, levels = ceil(log2(nrInputs)), i;
					
					i = levels;
					while(i>=0){
						totalMux += nrInputs/8;
						nrInputs /= 8;
						i--;
					}
					return totalMux;
				}
			}else if(id_ == "StratixIII"){
				if(nrInputs<=8)
					return ceil(nrInputs/4.0);
				else{
					int totalMux = 0, levels = ceil(log2(nrInputs)), i;
					
					i = levels;
					while(i>=0){
						totalMux += nrInputs/8;
						nrInputs /= 8;
						i--;
					}
					return totalMux;
				}
			}else if(id_ == "StratixIV"){
				if(nrInputs<=8)
					return ceil(nrInputs/4.0);
				else{
					int totalMux = 0, levels = ceil(log2(nrInputs)), i;
					
					i = levels;
					while(i>=0){
						totalMux += nrInputs/8;
						nrInputs /= 8;
						i--;
					}
					return totalMux;
				}
			}
		}
		
		cerr << "Warning: getLutFromMux(): unknown type of target" << endl;
		return 0;
	}
	
	
	//TODO: get more precise estimations for Virtex4, Virtex6, StratixII
	double Target::getLUTPerCounter(int width){
		
		if(vendor_ == "Xilinx"){
			if(id_ == "Spartan3"){
				if(width<=16)
					return 2.722;
				else if(width<=64)
					return 2.000;
				else if(width<=256)
					return 2.560;
			}else if(id_ == "Virtex4"){
				if(width<=16)
					return 1.055;
				else if(width<=64)
					return 1.510;
				else if(width<=256)
					return 1.945;
			}else if(id_ == "Virtex5"){
				if(width<=16)
					return 1.055;
				else if(width<=64)
					return 1.510;
				else if(width<=256)
					return 1.945;
			}else if(id_ == "Virtex6"){
				if(width<=16)
					return 1.055;
				else if(width<=64)
					return 1.510;
				else if(width<=256)
					return 1.945;
			}
		}else if(vendor_.compare("Altera") == 0){
			if(id_ == "StratixII"){
				if(width<=4)
					return 2.250;
				else if(width<=8)
					return 1.125;
				else if(width<=16)
					return 1.063;
				else if(width<=24)
					return 1.041;
				else if(width<=32)
					return 1.031;
				else if(width<=64)
					return 1.015;
			}else if(id_ == "StratixIII"){
				if(width<=4)
					return 2.250;
				else if(width<=8)
					return 1.125;
				else if(width<=16)
					return 1.063;
				else if(width<=24)
					return 1.041;
				else if(width<=32)
					return 1.031;
				else if(width<=64)
					return 1.015;
			}else if(id_ == "StratixIV"){
				if(width<=4)
					return 2.250;
				else if(width<=8)
					return 1.125;
				else if(width<=16)
					return 1.063;
				else if(width<=24)
					return 1.041;
				else if(width<=32)
					return 1.031;
				else if(width<=64)
					return 1.015;
			}
		}
		
		return 0;
	}
	
	
	//TODO: get more precise estimations for Virtex4, Virtex6, StratixII
	double Target::getFFPerCounter(int width){
		
		if(vendor_ == "Xilinx"){
			if(id_ == "Spartan3"){
				if(width<=16)
					return 2.000;
				else if(width<=64)
					return 2.170;
				else if(width<=256)
					return 2.320;
			}else if(id_ == "Virtex4"){
				if(width<=16)
					return 2.000;
				else if(width<=64)
					return 2.170;
				else if(width<=256)
					return 2.320;
			}else if(id_ == "Virtex5"){
				if(width<=16)
					return 1.000;
				else if(width<=64)
					return 1.829;
				else if(width<=256)
					return 2.095;
			}else if(id_ == "Virtex6"){
				if(width<=16)
					return 1.000;
				else if(width<=64)
					return 1.829;
				else if(width<=256)
					return 2.095;
			}
		}else if(vendor_.compare("Altera") == 0){
			if(id_ == "StratixII"){
				if(width<=4)
					return 1.000;
				else if(width<=8)
					return 1.000;
				else if(width<=16)
					return 1.000;
				else if(width<=24)
					return 1.000;
				else if(width<=32)
					return 1.000;
				else if(width<=64)
					return 1.000;
			}else if(id_ == "StratixIII"){
				if(width<=4)
					return 1.000;
				else if(width<=8)
					return 1.000;
				else if(width<=16)
					return 1.000;
				else if(width<=24)
					return 1.000;
				else if(width<=32)
					return 1.000;
				else if(width<=64)
					return 1.000;
			}else if(id_ == "StratixIV"){
				if(width<=4)
					return 1.000;
				else if(width<=8)
					return 1.000;
				else if(width<=16)
					return 1.000;
				else if(width<=24)
					return 1.000;
				else if(width<=32)
					return 1.000;
				else if(width<=64)
					return 1.000;
			}
		}
		
		return 0;
	}
	
	
	//TODO: get more precise estimations for Virtex4, Virtex6, StratixII
	//		get more specific estimates for Altera (several possibilities)
	double Target::getLUTPerAccumulator(int width, bool useDSP){
		
		if(vendor_ == "Xilinx"){
			if(id_ == "Spartan3"){
				if(useDSP){
					if(width<=8)
						return 0.000;
					else if(width<=32)
						return 0.000;
					else if(width<=64)
						return 0.000;
					else if(width<=128)
						return 0.000;
				}else{
					if(width<=8)
						return 1.000;
					else if(width<=32)
						return 1.187;
					else if(width<=64)
						return 2.218;
					else if(width<=128)
						return 3.050;
				}
			}else if(id_ == "Virtex4"){
				if(useDSP){
					if(width<=8)
						return 0.000;
					else if(width<=32)
						return 0.000;
					else if(width<=64)
						return 0.000;
					else if(width<=128)
						return 0.000;
				}else{
					if(width<=8)
						return 1.000;
					else if(width<=32)
						return 1.125;
					else if(width<=64)
						return 1.843;
					else if(width<=128)
						return 2.400;
				}
			}else if(id_ == "Virtex5"){
				if(useDSP){
					if(width<=8)
						return 0.000;
					else if(width<=32)
						return 0.000;
					else if(width<=64)
						return 0.000;
					else if(width<=128)
						return 0.000;
				}else{
					if(width<=8)
						return 1.000;
					else if(width<=32)
						return 1.125;
					else if(width<=64)
						return 1.843;
					else if(width<=128)
						return 2.400;
				}
			}else if(id_ == "Virtex6"){
				if(useDSP){
					if(width<=8)
						return 0.000;
					else if(width<=32)
						return 0.000;
					else if(width<=64)
						return 0.000;
					else if(width<=128)
						return 0.000;
				}else{
					if(width<=8)
						return 1.000;
					else if(width<=32)
						return 1.125;
					else if(width<=64)
						return 1.843;
					else if(width<=128)
						return 2.400;
				}
			}
		}else if(vendor_.compare("Altera") == 0){
			if(id_ == "StratixII"){
				if(width<=8)
					return 4.000;
				else if(width<=64)
					return 3.953;
			}else if(id_ == "StratixIII"){
				if(width<=8)
					return 4.000;
				else if(width<=64)
					return 3.953;
			}else if(id_ == "StratixIV"){
				if(width<=8)
					return 4.000;
				else if(width<=64)
					return 3.953;
			}
		}
		
		return 0;
	}
	
	
	//TODO: get more precise estimations for Virtex4, Virtex6, StratixII
	//		get more specific estimates for Altera (several possibilities)
	double Target::getFFPerAccumulator(int width, bool useDSP){
		
		if(vendor_ == "Xilinx"){
			if(id_ == "Spartan3"){
				if(useDSP){
					if(width<=8)
						return 0.000;
					else if(width<=32)
						return 0.000;
					else if(width<=64)
						return 0.000;
					else if(width<=128)
						return 0.000;
				}else{
					if(width<=8)
						return 1.000;
					else if(width<=32)
						return 2.906;
					else if(width<=64)
						return 3.765;
					else if(width<=128)
						return 4.190;
				}
			}else if(id_ == "Virtex4"){
				if(useDSP){
					if(width<=8)
						return 0.000;
					else if(width<=32)
						return 0.000;
					else if(width<=64)
						return 0.000;
					else if(width<=128)
						return 0.000;
				}else{
					if(width<=8)
						return 1.000;
					else if(width<=32)
						return 2.906;
					else if(width<=64)
						return 3.765;
					else if(width<=128)
						return 4.190;
				}
			}else if(id_ == "Virtex5"){
				if(useDSP){
					if(width<=8)
						return 0.000;
					else if(width<=32)
						return 0.000;
					else if(width<=64)
						return 0.000;
					else if(width<=128)
						return 0.000;
				}else{
					if(width<=8)
						return 1.000;
					else if(width<=32)
						return 2.281;
					else if(width<=64)
						return 3.390;
					else if(width<=128)
						return 3.870;
				}
			}else if(id_ == "Virtex6"){
				if(useDSP){
					if(width<=8)
						return 0.000;
					else if(width<=32)
						return 0.000;
					else if(width<=64)
						return 0.000;
					else if(width<=128)
						return 0.000;
				}else{
					if(width<=8)
						return 1.000;
					else if(width<=32)
						return 2.281;
					else if(width<=64)
						return 3.390;
					else if(width<=128)
						return 3.870;
				}
			}
		}else if(vendor_.compare("Altera") == 0){
			if(id_ == "StratixII"){
				if(width<=8)
					return 6.375;
				else if(width<=64)
					return 3.531;
			}else if(id_ == "StratixIII"){
				if(width<=8)
					return 6.375;
				else if(width<=64)
					return 3.531;
			}else if(id_ == "StratixIV"){
				if(width<=8)
					return 6.375;
				else if(width<=64)
					return 3.593;
			}
		}
		
		return 0;
	}
	
	
	double Target::getDSPPerAccumulator(int width){
		
		if(vendor_ == "Xilinx"){
			if(id_ == "Spartan3"){
				if(width<=48)
					return 1.000;
				else 
					return ceil(width/48.0);
			}else if(id_ == "Virtex4"){
				if(width<=48)
					return 1.000;
				else 
					return ceil(width/48.0);
			}else if(id_ == "Virtex5"){
				if(width<=48)
					return 1.000;
				else 
					return ceil(width/48.0);
			}else if(id_ == "Virtex6"){
				if(width<=48)
					return 1.000;
				else 
					return ceil(width/48.0);
			}
		}else if(vendor_.compare("Altera") == 0){
			if(id_ == "StratixII"){
				if(width<=38)
					return 1.000;
				else 
					return ceil(width/38.0);
			}else if(id_ == "StratixIII"){
				if(width<=38)
					return 1.000;
				else 
					return ceil(width/38.0);
			}else if(id_ == "StratixIV"){
				if(width<=38)
					return 1.000;
				else 
					return ceil(width/38.0);
			}
		}
		
		return 0;
	}
	
	//TODO: get more precise data for Xilinx architectures and for
	//		StratixII
	double Target::getLUTPerDecoder(int width){
		
		if(vendor_ == "Xilinx"){
			return 1.000;
		}else if(vendor_.compare("Altera") == 0){
			if(id_ == "StratixII"){
				if(width<=12)
					return 1.583;
				else if(width<=29)
					return 1.894;
				else if(width<=32)
					return 1.250;
				else if(width<=64)
					return 1.234;
			}else if(id_ == "StratixIII"){
				if(width<=12)
					return 1.583;
				else if(width<=29)
					return 1.894;
				else if(width<=32)
					return 1.250;
				else if(width<=64)
					return 1.234;
			}else if(id_ == "StratixIV"){
				if(width<=12)
					return 1.583;
				else if(width<=29)
					return 1.310;
				else if(width<=32)
					return 1.312;
				else if(width<=64)
					return 1.218;
			}
		}
		
		return 0;
	}
	
	
	//TODO: get more precise data for Xilinx architectures and for
	//		StratixII
	double Target::getFFPerDecoder(int width){
		
		if(vendor_ == "Xilinx"){
			return 1.000;
		}else if(vendor_.compare("Altera") == 0){
			if(id_ == "StratixII"){
				if(width<=12)
					return 2.500;
				else if(width<=29)
					return 2.241;
				else if(width<=32)
					return 2.218;
				else if(width<=64)
					return 2.125;
			}else if(id_ == "StratixIII"){
				if(width<=12)
					return 2.500;
				else if(width<=29)
					return 2.241;
				else if(width<=32)
					return 2.218;
				else if(width<=64)
					return 2.125;
			}else if(id_ == "StratixIV"){
				if(width<=12)
					return 2.500;
				else if(width<=29)
					return 2.241;
				else if(width<=32)
					return 2.218;
				else if(width<=64)
					return 2.125;
			}
		}
		
		return 0;
	}
	
	/*----------------------------------------------------------------*/

}















