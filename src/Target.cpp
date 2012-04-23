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

	/*-------------- Resource Estimation related items	----------*/
	
	
	virtual int Target::suggestLUTType(){
		cerr << "Error: function suggestLUTType() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	 
	virtual int Target::suggestLUTCount(int count, int type){
		cerr << "Error: function suggestLUTCount() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int Target::suggestDSPfromMultiplier(int count, int widthX, int widthY){
		cerr << "Error: function suggestDSPfromMultiplier() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int Target::suggestMultiplierCount(int count, int widthX, int widthY){
		cerr << "Error: function suggestMultiplierCount() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int Target::suggestMemoryCount(int count, int size, int width, int type){
		cerr << "Error: function suggestMemoryCount() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int Target::suggestSRLCount(int count, int width){
		cerr << "Error: function suggestSRLCount() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int Target::suggestLUTfromSRL(int count, int width){
		cerr << "Error: function suggestLUTfromSRL() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int Target::suggestFFfromSRL(int count, int width){
		cerr << "Error: function suggestFFfromSRL() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int Target::suggestMuxCount(int count, int nrInputs, int width){
		cerr << "Error: function suggestMuxCount() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int Target::suggestLUTfromMux(int count, int nrInputs, int width){
		cerr << "Error: function suggestLUTfromMux() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int Target::suggestLUTfromCounter(int count, int width){
		cerr << "Error: function suggestLUTfromCounter() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int Target::suggestFFfromCounter(int count, int width){
		cerr << "Error: function suggestFFfromCounter() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int Target::suggestLUTfromAccumulator(int count, int width){
		cerr << "Error: function suggestLUTfromAccumulator() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int Target::suggestFFfromAccumulator(int count, int width){
		cerr << "Error: function suggestFFfromAccumulator() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int Target::suggestDSPfromAccumulator(int count, int width){
		cerr << "Error: function suggestDSPfromAccumulator() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int Target::suggestLUTfromDecoder(int count, int width){
		cerr << "Error: function suggestLUTfromDecoder() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int Target::suggestFFfromDecoder(int count, int width){
		cerr << "Error: function suggestFFfromDecoder() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int Target::suggestRAMfromDecoder(int count, int width){
		cerr << "Error: function suggestRAMfromDecoder() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int Target::suggestLUTfromArithmeticOperator(int count, int nrInputs, int width){
		cerr << "Error: function suggestLUTfromArithmeticOperator() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int Target::suggestLUTfromFSM(int count, int nrStates, int nrTransitions){
		cerr << "Error: function suggestLUTfromFSM() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int Target::suggestFFfromFSM(int count, int nrStates, int nrTransitions){
		cerr << "Error: function suggestFFfromFSM() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int Target::suggestRAMfromFSM(int count, int nrStates, int nrTransitions){
		cerr << "Error: function suggestRAMfromFSM() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	/*------------------------------------------------------------*/

}
