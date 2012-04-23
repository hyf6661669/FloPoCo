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
	
	
	virtual int suggestLUTType(){
		cerr << "Error: function suggestLUTType() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	 
	virtual int suggestLUTCount(int count, int type){
		cerr << "Error: function suggestLUTCount() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int suggestDSPfromMultiplier(int count, int widthX, int widthY){
		cerr << "Error: function suggestDSPfromMultiplier() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int suggestMultiplierCount(int count, int widthX, int widthY){
		cerr << "Error: function suggestMultiplierCount() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int suggestMemoryCount(int count, int size, int width, int type){
		cerr << "Error: function suggestMemoryCount() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int suggestSRLCount(int count, int width){
		cerr << "Error: function suggestSRLCount() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int suggestLUTfromSRL(int count, int width){
		cerr << "Error: function suggestLUTfromSRL() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int suggestFFfromSRL(int count, int width){
		cerr << "Error: function suggestFFfromSRL() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int suggestMuxCount(int count, int nrInputs, int width){
		cerr << "Error: function suggestMuxCount() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int suggestLUTfromMux(int count, int nrInputs, int width){
		cerr << "Error: function suggestLUTfromMux() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int suggestLUTfromCounter(int count, int width){
		cerr << "Error: function suggestLUTfromCounter() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int suggestFFfromCounter(int count, int width){
		cerr << "Error: function suggestFFfromCounter() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int suggestLUTfromAccumulator(int count, int width){
		cerr << "Error: function suggestLUTfromAccumulator() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int suggestFFfromAccumulator(int count, int width){
		cerr << "Error: function suggestFFfromAccumulator() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int suggestDSPfromAccumulator(int count, int width){
		cerr << "Error: function suggestDSPfromAccumulator() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int suggestLUTfromDecoder(int count, int width){
		cerr << "Error: function suggestLUTfromDecoder() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int suggestFFfromDecoder(int count, int width){
		cerr << "Error: function suggestFFfromDecoder() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int suggestRAMfromDecoder(int count, int width){
		cerr << "Error: function suggestRAMfromDecoder() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int suggestLUTfromArithmeticOperator(int count, int nrInputs, int width){
		cerr << "Error: function suggestLUTfromArithmeticOperator() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int suggestLUTfromFSM(int count, int nrStates, int nrTransitions){
		cerr << "Error: function suggestLUTfromFSM() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int suggestFFfromFSM(int count, int nrStates, int nrTransitions){
		cerr << "Error: function suggestFFfromFSM() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	
	virtual int suggestRAMfromFSM(int count, int nrStates, int nrTransitions){
		cerr << "Error: function suggestRAMfromFSM() not implemented for this target FPGA." << endl;
		exit(1);
	}
	
	/*------------------------------------------------------------*/

}
