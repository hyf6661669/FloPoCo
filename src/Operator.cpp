/*
the base Operator class, every operator should inherit it

Author : Florent de Dinechin, Bogdan Pasca

Initial software.
Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
2008-2010.
  All rights reserved.

*/


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include "Operator.hpp"
#include "utils.hpp"


namespace flopoco{
	

	// global variables used through most of FloPoCo,
	// to be encapsulated in something, someday?
	int Operator::uid = 0; //init of the uid static member of Operator
	int verbose=0;
	
	Operator::Operator(Target* target, map<string, double> inputDelays){
		target_                     = target;
		numberOfInputs_             = 0;
		numberOfOutputs_            = 0;
		hasRegistersWithoutReset_   = false;
		hasRegistersWithAsyncReset_ = false;
		hasRegistersWithSyncReset_  = false;
		pipelineDepth_              = 0;
		currentCycle_               = 0;
		criticalPath_               = 0;
		needRecirculationSignal_    = false;
		inputDelayMap               = inputDelays;
		myuid                       = getNewUId();
		architectureName_			= "arch";
		indirectOperator_           =NULL;
		
		if (target_->isPipelined())
			setSequential();
		else
			setCombinatorial();	
		
		vhdl.disableParsing(!target_->isPipelined());	
	}
	
	
	void Operator::addInput(const std::string name, const int width, const bool isBus) {
		if (signalMap_.find(name) != signalMap_.end()) {
			std::ostringstream o;
			o << srcFileName << " (" << uniqueName_ << "): ERROR in addInput, signal " << name<< " seems to already exist";
			throw o.str();
		}
		Signal *s = new Signal(name, Signal::in, width, isBus) ; // default TTL and cycle OK
		s->setCycle(0);
		ioList_.push_back(s);
		signalMap_[name] = s ;
		numberOfInputs_ ++;
		declareTable[name] = s->getCycle();
	}
	
	void Operator::addOutput(const std::string name, const int width, const int numberOfPossibleOutputValues, const bool isBus) {
		if (signalMap_.find(name) != signalMap_.end()) {
			std::ostringstream o;
			o  << srcFileName << " (" << uniqueName_ << "): ERROR in addOutput, signal " << name << " seems to already exist";
			throw o.str();
		}
		Signal *s = new Signal(name, Signal::out, width, isBus) ;
		s -> setNumberOfPossibleValues(numberOfPossibleOutputValues);
		ioList_.push_back(s);
		for(int i=0; i<numberOfPossibleOutputValues; i++) 
			testCaseSignals_.push_back(s);
		signalMap_[name] = s ;
		numberOfOutputs_ ++;
		//		declareTable[name] = s->getCycle();
	}
	
	void Operator::addFPInput(const std::string name, const int wE, const int wF) {
		if (signalMap_.find(name) != signalMap_.end()) {
			std::ostringstream o;
			o << srcFileName << " (" << uniqueName_ << "): ERROR in addFPInput, signal " << name<< " seems to already exist";
			throw o.str();
		}
		Signal *s = new Signal(name, Signal::in, wE, wF);
		s->setCycle(0);
		ioList_.push_back(s);
		signalMap_[name] = s ;
		numberOfInputs_ ++;
		declareTable[name] = s->getCycle();
	}
	
	void Operator::addFPOutput(const std::string name, const int wE, const int wF, const int numberOfPossibleOutputValues) {
		if (signalMap_.find(name) != signalMap_.end()) {
			std::ostringstream o;
			o << srcFileName << " (" << uniqueName_ << "): ERROR in addFPOutput, signal " << name<< " seems to already exist";
			throw o.str();
		}
		Signal *s = new Signal(name, Signal::out, wE, wF) ;
		s -> setNumberOfPossibleValues(numberOfPossibleOutputValues);
		ioList_.push_back(s);
		for(int i=0; i<numberOfPossibleOutputValues; i++) 
			testCaseSignals_.push_back(s);
		signalMap_[name] = s ;
		numberOfOutputs_ ++;
		//		declareTable[name] = s->getCycle();
	}
	
	
	void Operator::addIEEEInput(const std::string name, const int wE, const int wF) {
		if (signalMap_.find(name) != signalMap_.end()) {
			std::ostringstream o;
			o << srcFileName << " (" << uniqueName_ << "): ERROR in addIEEEInput, signal " << name<< " seems to already exist";
			throw o.str();
		}
		Signal *s = new Signal(name, Signal::in, wE, wF, true);
		s->setCycle(0);
		ioList_.push_back(s);
		signalMap_[name] = s ;
		numberOfInputs_ ++;
		declareTable[name] = s->getCycle();
	}
	
	void Operator::addIEEEOutput(const std::string name, const int wE, const int wF, const int numberOfPossibleOutputValues) {
		if (signalMap_.find(name) != signalMap_.end()) {
			std::ostringstream o;
			o << srcFileName << " (" << uniqueName_ << "): ERROR in addIEEEOutput, signal " << name<< " seems to already exist";
			throw o.str();
		}
		Signal *s = new Signal(name, Signal::out, wE, wF, true) ;
		s -> setNumberOfPossibleValues(numberOfPossibleOutputValues);
		ioList_.push_back(s);
		for(int i=0; i<numberOfPossibleOutputValues; i++) 
			testCaseSignals_.push_back(s);
		signalMap_[name] = s ;
		numberOfOutputs_ ++;
		//		declareTable[name] = s->getCycle();
	}
	
	
	
	Signal * Operator::getSignalByName(string name) {
		ostringstream e;
		if(signalMap_.find(name) ==  signalMap_.end()) {
			e << srcFileName << " (" << uniqueName_ << "): ERROR in getSignalByName, signal " << name<< " not declared";
			throw e.str();
		}
		return signalMap_[name];
	}
	
	bool Operator::isSignalDeclared(string name){
		ostringstream e;
		if(signalMap_.find(name) ==  signalMap_.end()) {
			return false;
		}
		return true;
	}
	
	void Operator::setName(std::string prefix, std::string postfix){
		ostringstream pr, po;
		if (prefix.length()>0)
			pr << prefix << "_"; 
		else 
			pr << "";
		if (postfix.length()>0)
			po << "_"<<postfix;
		else
			po << "";
		uniqueName_ = pr.str() + uniqueName_ + po.str();
	}
	
	void Operator::setName(std::string operatorName){
		uniqueName_ = operatorName;
	}
	
	
	void  Operator::changeName(std::string operatorName){
		commentedName_ = uniqueName_;
		uniqueName_ = operatorName;
	}
	
	
	string Operator::getName() const{
		return uniqueName_;
	}
	
	int Operator::getIOListSize() const{
		return ioList_.size();
	}
	
	vector<Signal*> * Operator::getIOList(){
		return &ioList_; 
	}
	
	Signal * Operator::getIOListSignal(int i){
		return ioList_[i];
	}
	
	
	
	void  Operator::outputVHDLSignalDeclarations(std::ostream& o) {
		for (unsigned int i=0; i < this->signalList_.size(); i++){
			Signal* s = this->signalList_[i];
			o<<tab<<  s->toVHDL() << ";" << endl;
		}
		
	}
	
	void  Operator::outputVHDLRegisters(std::ostream& o) {
		unsigned int i;
		// execute only if the operator is sequential, otherwise output nothing
		if (isSequential()){
			// First registers without a reset
			if (hasRegistersWithoutReset_) {
				o << tab << "process(clk)  begin\n"
				<< tab << tab << "if clk'event and clk = '1' then\n";
				for(i=0; i<signalList_.size(); i++) {
					Signal *s = signalList_[i];
					if(s->type()==Signal::registeredWithoutReset) 
						o << tab <<tab << tab << s->getName() <<"_d <=  " << s->getName() <<";\n";
				}
				o << tab << tab << "end if;\n";
				o << tab << "end process;\n"; 
			}
			
			// then registers with a reset
			if (hasRegistersWithAsyncReset_) {
				o << tab << "process(clk, rst)" << endl;
				o << tab << tab << "begin" << endl;
				o << tab << tab << tab << "if rst = '1' then" << endl;
				for(i=0; i<signalList_.size(); i++) {
					Signal *s = signalList_[i];
					if(s->type()==Signal::registeredWithAsyncReset) {
						if ((s->width()>1)||(s->isBus())) 
							o << tab <<tab << tab << s->getName() <<"_d <=  (" << s->width()-1 <<" downto 0 => '0');\n";
						else
							o << tab <<tab << tab << s->getName() <<"_d <=  '0';\n";
					}
				}
				o << tab << tab << tab << "elsif clk'event and clk = '1' then" << endl;
				for(i=0; i<signalList_.size(); i++) {
					Signal *s = signalList_[i];
					if(s->type()==Signal::registeredWithAsyncReset) 
						o << tab <<tab << tab << s->getName() <<"_d <=  " << s->getName() <<";\n";
				}
				o << tab << tab << tab << "end if;" << endl;
				o << tab << tab << "end process;" << endl;
			}
			
			// then registers with synchronous reset
			if (hasRegistersWithSyncReset_) {
				o << tab << "process(clk, rst)" << endl;
				o << tab << tab << "begin" << endl;
				o<<  "    if clk'event and clk = '1' then" << endl;
				o << tab << tab << tab << "if rst = '1' then" << endl;
				for(i=0; i<signalList_.size(); i++) {
					Signal *s = signalList_[i];
					if(s->type()==Signal::registeredWithSyncReset) {
						if ((s->width()>1)||(s->isBus())) 
							o << tab <<tab << tab << s->getName() <<"_d <=  (" << s->width()-1 <<" downto 0 => '0');\n";
						else
							o << tab <<tab << tab << s->getName() <<"_d <=  '0';\n";
					}
				}
				o << tab << tab << tab << "else" << endl;
				for(i=0; i<signalList_.size(); i++) {
					Signal *s = signalList_[i];
					if(s->type()==Signal::registeredWithSyncReset) 
						o << tab <<tab << tab << s->getName() <<"_d <=  " << s->getName() <<";\n";
				}
				o << tab << tab << tab << "end if;" << endl;
				o << tab << tab << "end if;" << endl;
				o << tab << tab << "end process;" << endl;
			}
		}
	}
	
	
	void Operator::outputVHDLComponent(std::ostream& o, std::string name) {
		unsigned int i;
		o << tab << "component " << name << " is" << endl;
		if (ioList_.size() > 0)
		{
			o << tab << tab << "port ( ";
			if(isSequential()) {
				// add clk and rst signals which are no longer member of iolist
				o << "clk, rst : in std_logic;" <<endl;
			}
			if(isRecirculatory()) {
				// add clk and rst signals which are no longer member of iolist
				o << "stall_s : in std_logic;" <<endl;
			}
			for (i=0; i<this->ioList_.size(); i++){
				Signal* s = this->ioList_[i];
				if (i>0 || isSequential()) // align signal names 
					o<<tab<<"          ";
				o<<  s->toVHDL();
				if(i < this->ioList_.size()-1)  o<<";" << endl;
			}
			o << tab << ");"<<endl;
		}
		o << tab << "end component;" << endl;
	}
	
	void Operator::outputVHDLComponent(std::ostream& o) {
		this->outputVHDLComponent(o,  this->uniqueName_); 
	}
	
	
	void Operator::outputVHDLEntity(std::ostream& o) {
		unsigned int i;
		o << "entity " << uniqueName_ << " is" << endl;
		if (ioList_.size() > 0)
		{
			o << tab << "port ( ";
			
			if(isSequential()) {
				// add clk and rst signals which are no longer member of iolist
				o << "clk, rst : in std_logic;" <<endl;
			}
			if(isRecirculatory()) {
				// add stall signals to stop pipeline work 
				o << "stall_s : in std_logic;" <<endl;
			}
			
			for (i=0; i<this->ioList_.size(); i++){
				Signal* s = this->ioList_[i];
				if (i>0 || isSequential()) // align signal names 
					o<<"          ";
				o<<  s->toVHDL();
				if(i < this->ioList_.size()-1)  o<<";" << endl;
			}
			
			o << tab << ");"<<endl;
		}
		o << "end entity;" << endl << endl;
	}
	
	
	void Operator::setCopyrightString(std::string authorsYears){
		copyrightString_ = authorsYears;
	}
	
	
	void Operator::licence(std::ostream& o){
		licence(o, copyrightString_);
	}
	
	
	void Operator::licence(std::ostream& o, std::string authorsyears){
		o<<"--------------------------------------------------------------------------------"<<endl;
		// centering the unique name
		int s, i;
		if(uniqueName_.size()<76) s = (76-uniqueName_.size())/2; else s=0;
		o<<"--"; for(i=0; i<s; i++) o<<" "; o  << uniqueName_ << endl; 
		
		// if this operator was renamed from the command line, show the original name
		if(commentedName_!="") {
			if(commentedName_.size()<74) s = (74-commentedName_.size())/2; else s=0;
			o<<"--"; for(i=0; i<s; i++) o<<" "; o << "(" << commentedName_ << ")" << endl; 
		}
		
		o<<"-- This operator is part of the Infinite Virtual Library FloPoCoLib"<<endl;
		o<<"-- All rights reserved "<<endl;
		o<<"-- Authors: " << authorsyears <<endl;
		o<<"--------------------------------------------------------------------------------"<<endl;
	}
	
	
	
	void Operator::pipelineInfo(std::ostream& o){
		if(isSequential())
			o<<"-- Pipeline depth: " <<getPipelineDepth() << " cycles"  <<endl <<endl;
		else 
			o<<"-- combinatorial"  <<endl <<endl;
	}
	
	void Operator::outputVHDL(std::ostream& o) {
		this->outputVHDL(o, this->uniqueName_); 
	}
	
	bool Operator::isSequential() {
		return isSequential_; 
	}
	
	bool Operator::isRecirculatory() {
		return needRecirculationSignal_; 
	}
	
	void Operator::setSequential() {
		isSequential_=true; 
		vhdl.disableParsing(false); 
	}
	
	void Operator::setCombinatorial() {
		isSequential_=false;
		vhdl.disableParsing(true); 
	}
	
	void Operator::setRecirculationSignal() {
		needRecirculationSignal_ = true;
	}
	
	
	int Operator::getPipelineDepth() {
		return pipelineDepth_; 
	}
	
	void Operator::outputFinalReport(int level) {

		if (getIndirectOperator()!=NULL){ // interface operator
			if(getOpList().size()!=1){
				ostringstream o;
				o << "!?! Operator " << getUniqueName() << " is an interface operator with " << getOpList().size() << "children";
				throw o.str();
			}
			getOpListR()[0]->outputFinalReport(level);
		}

		else{ // Hard operator
			for (unsigned i=0; i< getOpList().size(); i++)
				if (! getOpListR().empty())
					getOpListR()[i]->outputFinalReport(level+1);	
			
			ostringstream tabs, ctabs;
			for (int i=0;i<level-1;i++){
				tabs << "|" << tab;
				ctabs << "|" << tab;
			}
			
			if (level>0){
				tabs << "|" << "---";
				ctabs << "|" << tab;
			}
			
			cerr << tabs.str() << "Entity " << uniqueName_ <<":"<< endl;
			if(this->getPipelineDepth()!=0)
				cerr << ctabs.str() << tab << "Pipeline depth = " << getPipelineDepth() << endl;
			else
				cerr << ctabs.str() << tab << "Not pipelined"<< endl;
		}
	}


	void Operator::setCycle(int cycle, bool report) {
		// lexing part
		vhdl.flush(currentCycle_);
		if(isSequential()) {
			currentCycle_=cycle;
			vhdl.setCycle(currentCycle_);
			if(report){
				vhdl << tab << "----------------Synchro barrier, entering cycle " << currentCycle_ << "----------------" << endl ;
			}
			// automatically update pipeline depth of the operator 
			if (currentCycle_ > pipelineDepth_) 
				pipelineDepth_ = currentCycle_;
		}
	}
	
	int Operator::getCurrentCycle(){
		return currentCycle_;
	} 
	
	void Operator::nextCycle(bool report) {
		// lexing part
		vhdl.flush(currentCycle_);
		
		if(isSequential()) {
			
			currentCycle_ ++; 
			vhdl.setCycle(currentCycle_);
			if(report)
				vhdl << tab << "----------------Synchro barrier, entering cycle " << currentCycle_ << "----------------" << endl;
			
			criticalPath_ = 0;
			// automatically update pipeline depth of the operator 
			if (currentCycle_ > pipelineDepth_) 
				pipelineDepth_ = currentCycle_;
		}
	}

	void Operator::previousCycle(bool report) {
		// lexing part
		vhdl.flush(currentCycle_);
		
		if(isSequential()) {
			
			currentCycle_ --; 
			vhdl.setCycle(currentCycle_);
			if(report)
				vhdl << tab << "----------------Synchro barrier, entering cycle " << currentCycle_ << "----------------" << endl;
			
		}
	}
	
	
	void Operator::setCycleFromSignal(string name, bool report) {
		setCycleFromSignal(name, 0.0, report);
	}
	
	
	void Operator::setCycleFromSignal(string name, double criticalPath, bool report) {
		// lexing part
		vhdl.flush(currentCycle_);
		
		ostringstream e;
		e << srcFileName << " (" << uniqueName_ << "): ERROR in setCycleFromSignal, "; // just in case
		
		if(isSequential()) {
			Signal* s;
			try {
				s=getSignalByName(name);
			}
			catch (string e2) {
				e << endl << tab << e2;
				throw e.str();
			}
			
			if( s->getCycle() < 0 ) {
				ostringstream o;
				o << "signal " << name<< " doesn't have (yet?) a valid cycle";
				throw o.str();
			} 
			
			currentCycle_ = s->getCycle();
			criticalPath_ = criticalPath;
			vhdl.setCycle(currentCycle_);
			
			if(report)
				vhdl << tab << "---------------- cycle " << currentCycle_ << "----------------" << endl ;
			// automatically update pipeline depth of the operator 
			if (currentCycle_ > pipelineDepth_) 
				pipelineDepth_ = currentCycle_;
		}
	}
	
	
	int Operator::getCycleFromSignal(string name, bool report) {
		// lexing part
		vhdl.flush(currentCycle_);
		
		ostringstream e;
		e << srcFileName << " (" << uniqueName_ << "): ERROR in getCycleFromSignal, "; // just in case
		
		if(isSequential()) {
			Signal* s;
			try {
				s=getSignalByName(name);
			}
			catch (string e2) {
				e << endl << tab << e2;
				throw e.str();
			}
			
			if( s->getCycle() < 0 ) {
				ostringstream o;
				o << "signal " << name<< " doesn't have (yet?) a valid cycle";
				throw o.str();
			} 
			
			return s->getCycle();
		}else{
			return 0; //if combinatorial everything happens at cycle 0
		}
	}
	
	
	bool Operator::syncCycleFromSignal(string name, bool report) {
		return(syncCycleFromSignal(name, 0.0, report));
	}



	bool Operator::syncCycleFromSignal(string name, double criticalPath, bool report) {

		bool advanceCycle = false;

		// lexing part
		vhdl.flush(currentCycle_);
		ostringstream e;
		e << srcFileName << " (" << uniqueName_ << "): ERROR in syncCycleFromSignal, "; // just in case
		
		if(isSequential()) {
			Signal* s;
			try {
				s=getSignalByName(name);
			}
			catch (string e2) {
				e << endl << tab << e2;
				throw e.str();
			}
			
			if( s->getCycle() < 0 ) {
				ostringstream o;
				o << "signal " << name << " doesn't have (yet?) a valid cycle";
				throw o.str();
			} 

			if (s->getCycle() == currentCycle_){
				advanceCycle = false;
				if (criticalPath>criticalPath_)
					criticalPath_=criticalPath ;
			}
						
			if (s->getCycle() > currentCycle_){
				advanceCycle = true;
				currentCycle_ = s->getCycle();
				criticalPath_= criticalPath;
				vhdl.setCycle(currentCycle_);
			}
			
			// if (s->getCycle() < currentCycle_) do nothing: 
			//   the argument signal will be delayed, so its critical path will be 0

			// cout << tab << "----------------Synchro barrier on " << s->getName() << ",  entering cycle " << currentCycle_ << "----------------"  ;

			if(report && advanceCycle)
				vhdl << tab << "----------------Synchro barrier, entering cycle " << currentCycle_ << "----------------" << endl ;

			// automatically update pipeline depth of the operator 
			if (currentCycle_ > pipelineDepth_) 
				pipelineDepth_ = currentCycle_;
		}
		
		return advanceCycle;
	}

	void Operator::setSignalDelay(string name, double delay){
		Signal* s;
		try {
			s=getSignalByName(name);
		}
		catch (string e2) {
			cout << "WARNING: signal " << name << " was not found in file " << srcFileName << " when called using setSignalDelay" << endl;
			return;
		}

		s->setDelay(delay);		
	}	

	double Operator::getSignalDelay(string name){
		Signal* s;
		try {
			s=getSignalByName(name);
		}
		catch (string e2) {
			cout << "WARNING: signal " << name << " was not found in file " << srcFileName << " when called using getSignalDelay" << endl;
			return 0.0;
		}

		return s->getDelay();		
	}

	double Operator::getCriticalPath() {return criticalPath_;}
	
	void Operator::setCriticalPath(double delay) {criticalPath_=delay;}
	
	void Operator::addToCriticalPath(double delay){
		criticalPath_ += delay;
	}
	
//	bool Operator::manageCriticalPath(double delay, bool report){
//		//		criticalPath_ += delay;
//		if ( target_->ffDelay() + (criticalPath_ + delay) + target_->localWireDelay() > (1.0/target_->frequency())){
//			nextCycle(report); //TODO Warning
//			criticalPath_ = min(delay, 1.0/target_->frequency());
//			return true;
//		}
//		else{
//			criticalPath_ += delay;
//			return false;
//		}
//	}

	bool Operator::manageCriticalPath(double delay, bool report){
		//		criticalPath_ += delay;
		if ( target_->ffDelay() + (criticalPath_ + delay) > (1.0/target_->frequency())){
			nextCycle(report); //TODO Warning
			criticalPath_ = min(delay, 1.0/target_->frequency());
			return true;
		}
		else{
			criticalPath_ += delay;
			return false;
		}
	}

	
	double Operator::getOutputDelay(string s) {return outDelayMap[s];}  // TODO add checks
	
	string Operator::declare(string name, const int width, bool isbus, Signal::SignalType regType) {
		Signal* s;
		ostringstream e;
		// check the signals doesn't already exist
		if(signalMap_.find(name) !=  signalMap_.end()) {
			e << srcFileName << " (" << uniqueName_ << "): ERROR in declare(), signal " << name<< " already exists";
			throw e.str();
		}
		// construct the signal (lifeSpan and cycle are reset to 0 by the constructor)
		s = new Signal(name, regType, width, isbus);
		if(regType==Signal::registeredWithoutReset)
			hasRegistersWithoutReset_ = true;
		if(regType==Signal::registeredWithSyncReset)
			hasRegistersWithSyncReset_ = true;
		if(regType==Signal::registeredWithAsyncReset)
			hasRegistersWithAsyncReset_ = true;
		
		// define its cycle 
		if(isSequential())
			s->setCycle(this->currentCycle_);
		
		// add this signal to the declare table
		declareTable[name] = s->getCycle();
		
		// add the signal to signalMap and signalList
		signalList_.push_back(s);    
		signalMap_[name] = s ;
		return name;
	}
	
	
	#if 1
	string Operator::use(string name) {
		ostringstream e;
		e << "ERROR in use(), "; // just in case
		
		if(isSequential()) {
			Signal *s;
			try {
				s=getSignalByName(name);
			}
			catch (string e2) {
				e << endl << tab << e2;
				throw e.str();
			}
			if(s->getCycle() < 0) {
				e << "signal " << name<< " doesn't have (yet?) a valid cycle";
				throw e.str();
			} 
			if(s->getCycle() > currentCycle_) {
				ostringstream e;
				e << "active cycle of signal " << name<< " is later than current cycle, cannot delay it";
				throw e.str();
			} 
			// update the lifeSpan of s
			s->updateLifeSpan( currentCycle_ - s->getCycle() );
			//return s->delayedName( currentCycle_ - s->getCycle() );
			return s->delayedName( 0 );
		}
		else
			return name;
	}
	
	string Operator::use(string name, int delay) {
		
		ostringstream e;
		e << "ERROR in use(), "; // just in case
		
		if(isSequential()) {
			Signal *s;
			try {
				s=getSignalByName(name);
			}
			catch (string e2) {
				e << endl << tab << e2;
				throw e.str();
			}
			// update the lifeSpan of s
			
			s->updateLifeSpan( delay );
			//return s->delayedName( currentCycle_ - s->getCycle() );
			return s->delayedName( delay );
		}else
			return name;
	}
	
	#endif
	
	void Operator::outPortMap(Operator* op, string componentPortName, string actualSignalName, bool newSignal){
		Signal* formal;
		Signal* s;
		ostringstream e;
		e << srcFileName << " (" << uniqueName_ << "): ERROR in outPortMap() for entity " << op->getName()  << ", "; // just in case
		// check the signals doesn't already exist
		if(signalMap_.find(actualSignalName) !=  signalMap_.end() && newSignal) {
			e << "signal " << actualSignalName << " already exists";
			throw e.str();
		}
		try {
			formal=op->getSignalByName(componentPortName);
		}
		catch (string e2) {
			e << endl << tab << e2;
			throw e.str();
		}
		if (formal->type()!=Signal::out){
			e << "signal " << componentPortName << " of component " << op->getName() 
			<< " doesn't seem to be an output port";
			throw e.str();
		}
		if (newSignal) {
			int width = formal -> width();
			bool isbus = formal -> isBus();
			// construct the signal (lifeSpan and cycle are reset to 0 by the constructor)
			s = new Signal(actualSignalName, Signal::wire, width, isbus);
			// define its cycle 
			if(isSequential())
				s->setCycle( this->currentCycle_ + op->getPipelineDepth() );
		
			// add this signal to the declare table
			declareTable[actualSignalName] = s->getCycle();
			
			// add the signal to signalMap and signalList
			signalList_.push_back(s);    
			signalMap_[actualSignalName] = s ;
		};
		// add the mapping to the mapping list of Op
		op->portMap_[componentPortName] = actualSignalName;
	}
	
	
	void Operator::inPortMap(Operator* op, string componentPortName, string actualSignalName){
		Signal* formal;
		ostringstream e;
		string name;
		e  << srcFileName << " (" << uniqueName_ << "): ERROR in inPortMap() for entity " << op->getName() << ","; // just in case
		
		if(isSequential()) {
			Signal *s;
			try {
				s=getSignalByName(actualSignalName);
			}
			catch (string e2) {
				e << endl << tab << e2;
				throw e.str();
			}
			if(s->getCycle() < 0) {
				ostringstream e;
				e << "signal " << actualSignalName<< " doesn't have (yet?) a valid cycle";
				throw e.str();
			} 
			if(s->getCycle() > currentCycle_) {
				ostringstream e;
				e << "active cycle of signal " << actualSignalName<< " is later than current cycle, cannot delay it";
				throw e.str();
			} 
			// update the lifeSpan of s
			s->updateLifeSpan( currentCycle_ - s->getCycle() );
			name = s->delayedName( currentCycle_ - s->getCycle() );
		}
		else
			name = actualSignalName;
		
		try {
			formal=op->getSignalByName(componentPortName);
		}
		catch (string e2) {
			e << endl << tab << e2;
			throw e.str();
		}
		if (formal->type()!=Signal::in){
			e << "signal " << componentPortName << " of component " << op->getName() 
			<< " doesn't seem to be an input port";
			throw e.str();
		}
		
		// add the mapping to the mapping list of Op
		op->portMap_[componentPortName] = name;
	}
	
	
	
	void Operator::inPortMapCst(Operator* op, string componentPortName, string actualSignal){
		Signal* formal;
		ostringstream e;
		string name;
		e << srcFileName << " (" << uniqueName_ << "): ERROR in inPortMapCst() for entity " << op->getName()  << ", "; // just in case
		
		try {
			formal=op->getSignalByName(componentPortName);
		}
		catch (string e2) {
			e << endl << tab << e2;
			throw e.str();
		}
		if (formal->type()!=Signal::in){
			e << "signal " << componentPortName << " of component " << op->getName() 
			<< " doesn't seem to be an input port";
			throw e.str();
		}
		
		// add the mapping to the mapping list of Op
		op->portMap_[componentPortName] = actualSignal;
	}
	
	
	string Operator::instance(Operator* op, string instanceName){
		ostringstream o;
		// TODO add checks here? Check that all the signals are covered for instance
		
		o << tab << instanceName << ": " << op->getName();
		if (op->isSequential()) 
			o << "  -- pipelineDepth="<< op->getPipelineDepth() << " maxInDelay=" << getMaxInputDelays(op->inputDelayMap);
		o << endl;
		o << tab << tab << "port map ( ";
		// build vhdl and erase portMap_
		map<string,string>::iterator it;
		if(op->isSequential()) {
			o << "clk  => clk";
			o << "," << endl << tab << tab << "           rst  => rst";
		}
		if (op->isRecirculatory()) {
			o << "," << endl << tab << tab << "           stall_s => stall_s";
		};
		
		for (it=op->portMap_.begin()  ; it != op->portMap_.end(); it++ ) {
			bool outputSignal = false;
			for ( int k = 0; k < int(op->ioList_.size()); k++){
				if ((op->ioList_[k]->type() == Signal::out) && ( op->ioList_[k]->getName() == (*it).first )){ 
					outputSignal = true;
				}
			}
			
			bool parsing = vhdl.isParsing();
			
			if ( outputSignal && parsing){
				vhdl.flush(currentCycle_);
				vhdl.disableParsing(true);
			}
			
			if (it!=op->portMap_.begin() || op->isSequential())				
				o << "," << endl <<  tab << tab << "           ";
				
				o << (*it).first << " => "  << (*it).second;
			
			if ( outputSignal && parsing ){
				vhdl << o.str();
				vhdl.flush(currentCycle_);
				o.str("");
				vhdl.disableParsing(!parsing);
			}
			//op->portMap_.erase(it);
		}
		o << ");" << endl;
		
		
		
		//Floorplanning related-----------------------------------------
		/*
		floorplan << manageFloorplan();
		flComponentList.push_back(op->getName());
		flInstanceNames[op->getName()] = instanceName;
		*/
		//--------------------------------------------------------------
		
		
		
		// add the operator to the subcomponent list 
		subComponents_[op->getName()]  = op;
		return o.str();
	}


	
	string Operator::buildVHDLSignalDeclarations() {
		ostringstream o;
		for(unsigned int i=0; i<signalList_.size(); i++) {
			Signal *s = signalList_[i];
			o << s->toVHDLDeclaration() << endl;
		}
		//now the signals from the I/O List which have the cycle>0
		for (unsigned int i=0; i<ioList_.size(); i++) {
			Signal *s = ioList_[i];
			if (s->getLifeSpan()>0){
				o << s->toVHDLDeclaration() << endl;	
			}
			
		}
		
		return o.str();	
	}
	
	
	void Operator::useHardRAM(Operator* t) {
		if (target_->getVendor() == "Xilinx") 
		{
			addAttribute("rom_extract", "string", t->getName()+": component", "yes");
			addAttribute("rom_style", "string", t->getName()+": component", "block");
		}
		if (target_->getVendor() == "Altera") 
			addAttribute("altera_attribute", "string", t->getName()+": component", "-name ALLOW_ANY_ROM_SIZE_FOR_RECOGNITION ON");
	}
	
	void Operator::useSoftRAM(Operator* t) {
		if (target_->getVendor() == "Xilinx") 
		{
			addAttribute("rom_extract", "string", t->getName()+": component", "yes");
			addAttribute("rom_style", "string", t->getName()+": component", "distributed");
		}
		if (target_->getVendor() == "Altera") 
			addAttribute("altera_attribute", "string", t->getName()+": component", "-name ALLOW_ANY_ROM_SIZE_FOR_RECOGNITION OFF");
	}
	
	
	
	string Operator::buildVHDLComponentDeclarations() {
		ostringstream o;
		for(map<string, Operator*>::iterator it = subComponents_.begin(); it !=subComponents_.end(); it++) {
			Operator *op = it->second;
			op->outputVHDLComponent(o);
			o<< endl;
		}
		return o.str();	
	}
	
	
	void Operator::addConstant(std::string name, std::string t, mpz_class v) {
		ostringstream tmp; 
		tmp << v;
		constants_[name] =  make_pair(t, tmp.str());
	}
	
	void Operator::addType(std::string name, std::string value) {
		types_ [name] =  value;
	}
	
	void Operator::addConstant(std::string name, std::string t, int v) {
		ostringstream tmp; 
		tmp << v;
		constants_[name] =  make_pair(t, tmp.str());
	}
	
	void Operator::addConstant(std::string name, std::string t, string v) {
		constants_[name] =  make_pair(t, v);
	}
	
	
	void Operator::addAttribute(std::string attributeName,  std::string attributeType,  std::string object, std::string value ) {
		// TODO add some checks ?
		attributes_[attributeName] = attributeType;
		pair<string,string> p = make_pair(attributeName,object);
		attributesValues_[p] = value;
	}
	
	
	string Operator::buildVHDLTypeDeclarations() {
		ostringstream o;
		for(map<string, string >::iterator it = types_.begin(); it !=types_.end(); it++) {
			string name  = it->first;
			string value = it->second;
			o <<  "type " << name << " is "  << value << ";" << endl;
		}
		return o.str();	
	}
	
	
	string Operator::buildVHDLConstantDeclarations() {
		ostringstream o;
		for(map<string, pair<string, string> >::iterator it = constants_.begin(); it !=constants_.end(); it++) {
			string name  = it->first;
			string type = it->second.first;
			string value = it->second.second;
			o <<  "constant " << name << ": " << type << " := " << value << ";" << endl;
		}
		return o.str();	
	}
	
	
	
	string Operator::buildVHDLAttributes() {
		ostringstream o;
		// First add the declarations of attribute names
		for(map<string, string>::iterator it = attributes_.begin(); it !=attributes_.end(); it++) {
			string name  = it->first;
			string type = it->second;
			o <<  "attribute " << name << ": " << type << ";" << endl;
		}
		// Then add the declarations of attribute values
		for(map<pair<string, string>, string>::iterator it = attributesValues_.begin(); it !=attributesValues_.end(); it++) {
			string name  = it->first.first;
			string object = it->first.second;
			string value = it->second;
			if(attributes_[name]=="string")
				value = '"' + value + '"';
			o <<  "attribute " << name << " of " << object << " is " << value << ";" << endl;
		}
		return o.str();	
	}
	
	string  Operator::buildVHDLRegisters() {
		ostringstream o;
		
		// execute only if the operator is sequential, otherwise output nothing
		string recTab = "";
		if (isRecirculatory()) recTab = tab;
		if (isSequential()){
			o << tab << "process(clk)" << endl;
			o << tab << tab << "begin" << endl;
			o << tab << tab << tab << "if clk'event and clk = '1' then" << endl;
			if (isRecirculatory()) o << tab << tab << tab << tab << "if stall_s = '0' then" << endl;
			for(unsigned int i=0; i<signalList_.size(); i++) {
				Signal *s = signalList_[i];
				if ((s->type() == Signal::registeredWithoutReset) || (s->type() == Signal::wire)) 
					if(s->getLifeSpan() >0) {
						for(int j=1; j <= s->getLifeSpan(); j++)
							
							o << recTab << tab << tab <<tab << tab << s->delayedName(j) << " <=  " << s->delayedName(j-1) <<";" << endl;
					}
			}
			for(unsigned int i=0; i<ioList_.size(); i++) {
				Signal *s = ioList_[i];
				if(s->getLifeSpan() >0) {
					for(int j=1; j <= s->getLifeSpan(); j++)
						o << recTab << tab << tab <<tab << tab << s->delayedName(j) << " <=  " << s->delayedName(j-1) <<";" << endl;
				}
			}
			if (isRecirculatory()) o << tab << tab << tab << tab << "end if;" << endl;
			o << tab << tab << tab << "end if;\n";
			o << tab << tab << "end process;\n"; 
			
			// then registers with a reset
			if (hasRegistersWithAsyncReset_) {
				o << tab << "process(clk, rst)" << endl;
				o << tab << tab << "begin" << endl;
				o << tab << tab << tab << "if rst = '1' then" << endl;
				for(unsigned int i=0; i<signalList_.size(); i++) {
					Signal *s = signalList_[i];
					if (s->type() == Signal::registeredWithAsyncReset)  
						if(s->getLifeSpan() >0) {
							for(int j=1; j <= s->getLifeSpan(); j++){
								if ( (s->width()>1) || (s->isBus()))
									o << tab << tab <<tab << tab << s->delayedName(j) << " <=  (others => '0');" << endl;
								else
									o << tab <<tab << tab << tab << s->delayedName(j) << " <=  '0';" << endl;
							}
						}
				}			
				o << tab << tab << tab << "elsif clk'event and clk = '1' then" << endl;
				for(unsigned int i=0; i<signalList_.size(); i++) {
					Signal *s = signalList_[i];
					if (s->type() == Signal::registeredWithAsyncReset)  
						if(s->getLifeSpan() >0) {
							for(int j=1; j <= s->getLifeSpan(); j++)
								o << tab <<tab << tab << tab << s->delayedName(j) << " <=  " << s->delayedName(j-1) <<";" << endl;
						}
				}			o << tab << tab << tab << "end if;" << endl;
				o << tab << tab <<"end process;" << endl;
			}
			
			// then registers with synchronous reset
			if (hasRegistersWithSyncReset_) {
				o << tab << "process(clk, rst)" << endl;
				o << tab << tab << "begin" << endl;
				o << tab << tab << tab << "if clk'event and clk = '1' then" << endl;
				o << tab << tab << tab << tab << "if rst = '1' then" << endl;
				for(unsigned int i=0; i<signalList_.size(); i++) {
					Signal *s = signalList_[i];
					if (s->type() == Signal::registeredWithSyncReset)  
						if(s->getLifeSpan() >0) {
							for(int j=1; j <= s->getLifeSpan(); j++){
								if ( (s->width()>1) || (s->isBus()))
									o << tab <<tab << tab <<tab << tab << s->delayedName(j) << " <=  (others => '0');" << endl;
								else
									o << tab <<tab << tab <<tab << tab << s->delayedName(j) << " <=  '0';" << endl;
							}
						}
				}			
				o << tab << tab << tab << tab << "else" << endl;
				for(unsigned int i=0; i<signalList_.size(); i++) {
					Signal *s = signalList_[i];
					if (s->type() == Signal::registeredWithSyncReset)  
						if(s->getLifeSpan() >0) {
							for(int j=1; j <= s->getLifeSpan(); j++)
								o << tab <<tab << tab <<tab << tab << s->delayedName(j) << " <=  " << s->delayedName(j-1) <<";" << endl;
						}
				}			
				o << tab << tab << tab << tab << "end if;" << endl;
				o << tab << tab << tab << "end if;" << endl;
				o << tab << tab << "end process;" << endl;
			}
		}
		return o.str();
	}
	
	
	void Operator::buildStandardTestCases(TestCaseList* tcl) {
		// Each operator should overload this method. If not, it is mostly harmless but deserves a warning.
		cerr << "WARNING: No standard test cases implemented for this operator" << endl;
	}
	
	
	
	
	void Operator::buildRandomTestCaseList(TestCaseList* tcl, int n){
		
		TestCase *tc;
		/* Generate test cases using random input numbers */
		for (int i = 0; i < n; i++) {
			// TODO free all this memory when exiting TestBench
			tc = buildRandomTestCase(i); 
			tcl->add(tc);
		}
	}
	
	TestCase* Operator::buildRandomTestCase(int i){
		TestCase *tc = new TestCase(this);
		/* Generate test cases using random input numbers */
		// TODO free all this memory when exiting TestBench
		/* Fill inputs */
		for (unsigned int j = 0; j < ioList_.size(); j++) {
			Signal* s = ioList_[j]; 
			if (s->type() == Signal::in) {
				mpz_class a = getLargeRandom(s->width());
				tc->addInput(s->getName(), a);
			}
		}
		/* Get correct outputs */
		emulate(tc);
		
		//		cout << tc->getInputVHDL();
		//    cout << tc->getExpectedOutputVHDL();
		
		
		// add to the test case list
		return tc;
	}
	
	map<string, double> Operator::getOutDelayMap(){
		return outDelayMap;
	}
	
	map<string, int> Operator::getDeclareTable(){
		return declareTable;
	}
	
	void Operator::outputVHDL(std::ostream& o, std::string name) {
		if (! vhdl.isEmpty() ){
			licence(o);
			pipelineInfo(o);
			stdLibs(o);
			outputVHDLEntity(o);
			newArchitecture(o,name);
			o << buildVHDLComponentDeclarations();	
			o << buildVHDLSignalDeclarations();
			o << buildVHDLTypeDeclarations();
			o << buildVHDLConstantDeclarations();
			o << buildVHDLAttributes();
			beginArchitecture(o);		
			o<<buildVHDLRegisters();
			if(getIndirectOperator())
				o << getIndirectOperator()->vhdl.str();
			else
				o << vhdl.str();
			endArchitecture(o);
		}
	}
	
	void Operator::parse2(){
		REPORT(DEBUG, "Starting second-level parsing for operator "<<srcFileName);
		vector<pair<string,int> >:: iterator iterUse;
		map<string, int>::iterator iterDeclare;
		
		string name;
		int declareCycle, useCycle;
		
		string str (vhdl.str());
		
		/* parse the useTable and check that the declarations are ok */
		for (iterUse = (vhdl.useTable).begin(); iterUse!=(vhdl.useTable).end();++iterUse){
			name     = (*iterUse).first;
			useCycle = (*iterUse).second;
			
			ostringstream tSearch;
			ostringstream tReplace;
			string replaceString;
			
			tSearch << "__"<<name<<"__"<<useCycle<<"__"; 
			string searchString (tSearch.str());
			
			iterDeclare = declareTable.find(name);
			declareCycle = iterDeclare->second;
			
			if (iterDeclare != declareTable.end()){
				tReplace << use(name, useCycle - declareCycle); 
				replaceString = tReplace.str();
				if (useCycle<declareCycle){
					cerr << srcFileName << " (" << uniqueName_ << "): WARNING: Signal " << name <<" defined @ cycle "<<declareCycle<<" and used @ cycle " << useCycle <<endl;
					cerr << srcFileName << " (" << uniqueName_ << "): If this is a feedback signal you may ignore this warning"<<endl;
				}
			}else{
				/* parse the declare by hand and check lower/upper case */
				bool found = false;
				string tmp;
				for (iterDeclare = declareTable.begin(); iterDeclare!=declareTable.end();++iterDeclare){
					tmp = iterDeclare->first;
					if ( (to_lowercase(tmp)).compare(to_lowercase(name))==0){
						found = true;
						break;
					}
				}
				
				if (found == true){
					cerr  << srcFileName << " (" << uniqueName_ << "): ERROR: Clash on signal:"<<name<<". Definition used signal name "<<tmp<<". Check signal case!"<<endl;
					exit(-1);
				}
				
				tReplace << name;
				replaceString = tReplace.str();
			}
			
			if ( searchString != replaceString ){
				string::size_type pos = 0;
				while ( (pos = str.find(searchString, pos)) != string::npos ) {
					str.replace( pos, searchString.size(), replaceString );
					pos++;
				}
			}				
		}
		for (iterDeclare = declareTable.begin(); iterDeclare!=declareTable.end();++iterDeclare){
			name = iterDeclare->first;
			useCycle = iterDeclare->second;
			
			ostringstream tSearch;
			tSearch << "__"<<name<<"__"<<useCycle<<"__"; 
			//			cout << "searching for: " << tSearch.str() << endl;
			string searchString (tSearch.str());
			
			ostringstream tReplace;
			tReplace << name;
			string replaceString(tReplace.str());
			
			if ( searchString != replaceString ){
				
				string::size_type pos = 0;
				while ( (pos = str.find(searchString, pos)) != string::npos ) {
					str.replace( pos, searchString.size(), replaceString );
					pos++;
				}
			}				
		}
		vhdl.setSecondLevelCode(str);
	}
	
	
	void Operator::emulate(TestCase * tc) {
		throw std::string("emulate() not implemented for ") + uniqueName_;
	}
	
	bool Operator::hasComponent(string s){
		map<string, Operator*>::iterator theIterator;
		
		theIterator = subComponents_.find(s);
		if (theIterator != subComponents_.end() )
			return true;
		else
			return false;
	}

	void Operator::addComment(string comment, string align){
		vhdl << align << "-- " << comment << endl;
	}

	void Operator::addFullComment(string comment, int lineLength) {
		string align = "--";
		// - 2 for the two spaces
		for (unsigned i = 2; i < (lineLength - 2- comment.size()) / 2; i++) align += "-";
		vhdl << align << " " << comment << " " << align << endl; 
	}
	

	void Operator::outputVHDLToFile(vector<Operator*> &oplist, ofstream& file){
		string srcFileName = "Operator.cpp"; // for REPORT
		for(unsigned i=0; i<oplist.size(); i++) {
			try {
				REPORT(FULL, "OPERATOR:"<<oplist[i]->getName());
				REPORT(FULL, "--DECLARE LIST---------------------------------------------------");
				REPORT(FULL, printMapContent(oplist[i]->getDeclareTable()) );
				REPORT(FULL, "--USE LIST-------------------------------------------------------");
				REPORT(FULL, printVectorContent(  (oplist[i]->getFlopocoVHDLStream())->getUseTable()) );
				
				// check for subcomponents 
				if (! oplist[i]->getOpListR().empty() ){
					//recursively call to print subcomponent
					outputVHDLToFile(oplist[i]->getOpListR(), file);	
				}
				oplist[i]->getFlopocoVHDLStream()->flush();

				/* second parse is only for sequential operators */
				if (oplist[i]->isSequential()){
					REPORT (FULL, "--2nd PASS-------------------------------------------------------");
					oplist[i]->parse2();
				}
				oplist[i]->outputVHDL(file);			
			
			} catch (std::string s) {
					cerr << "Exception while generating '" << oplist[i]->getName() << "': " << s <<endl;
			}
		}
	}
	

	void Operator::outputVHDLToFile(ofstream& file){
		vector<Operator*> oplist;
		oplist.push_back(this);
		Operator::outputVHDLToFile(oplist, file);
	}
	
	
	/////////////////////////////////////////////////////////////////////////////////////////////////
	////////////Functions used for resource estimations
	
	//--Logging functions
	//---General resources
	void Operator::initResourceEstimation(){
		
		resourceEstimate << "Starting Resource estimation report for entity: " << uniqueName_ << " --------------- " << endl;
		resourceEstimateReport << "";
		
		reHelper = new ResourceEstimationHelper(target_, this);
		reHelper->initResourceEstimation();
	}
	
	std::string Operator::addFF(int count){
		
		return reHelper->addFF(count);
	}
	
	std::string Operator::addLUT(int nrInputs, int count){
		
		return reHelper->addLUT(nrInputs, count);
	}
	
	std::string Operator::addReg(int width, int count){
				
		return reHelper->addReg(width, count);
	}
	
	//TODO: verify increase in the DSP count
	std::string Operator::addMultiplier(int count){
		
		return reHelper->addMultiplier(count);
	}
	
	//TODO: verify increase in the DSP count 
	std::string Operator::addMultiplier(int widthX, int widthY, double ratio, int count){
		
		return reHelper->addMultiplier(widthX, widthY, ratio, count);
	}
	
	//TODO: verify increase in the element count
	std::string Operator::addAdderSubtracter(int widthX, int widthY, double ratio, int count){
		
		return reHelper->addAdderSubtracter(widthX, widthY, ratio, count);
	}
	
	//TODO: take into account the memory type (RAM or ROM); depending on 
	//		the type, might be implemented through distributed memory or
	//		dedicated memory blocks
	std::string Operator::addMemory(int size, int width, int type, int count){
		
		return reHelper->addMemory(size, width, type, count);
	}
	
	//---More particular resource logging
	std::string Operator::addDSP(int count){
		
		return reHelper->addDSP(count);
	}
	
	std::string Operator::addRAM(int count){
		
		return reHelper->addRAM(count);
	}
	
	std::string Operator::addROM(int count){
		
		return reHelper->addROM(count);
	}
	
	//TODO: should count the shift registers according to their bitwidths
	std::string Operator::addSRL(int width, int depth, int count){
				
		return reHelper->addSRL(width, depth, count);
	}
	
	std::string Operator::addWire(int count, std::string signalName){
		
		return reHelper->addWire(count, signalName);
	}
	
	std::string Operator::addIOB(int count, std::string portName){
		
		return reHelper->addIOB(count, portName);
	}
	
	//---Even more particular resource logging-------------------------
	
	//TODO: get a more accurate count of the number of multiplexers 
	//		needed; currently specific resources are not taken into account
	std::string Operator::addMux(int width, int nrInputs, int count){
				
		return reHelper->addMux(width, nrInputs, count);
	}
	
	//TODO: count the counters according to their bitwidth
	//TODO: get estimations when using specific resources (like DSPs)
	//		involves also changes to getLUTPerCounter() getFFPerCounter()
	std::string Operator::addCounter(int width, int count){
				
		return reHelper->addCounter(width, count);
	}
	
	//TODO: count the accumulators according to their bitwidth
	std::string Operator::addAccumulator(int width, bool useDSP, int count){
		
		return reHelper->addAccumulator(width, useDSP, count);
	}
	
	//TODO: count the decoders according to their input and output 
	//		bitwidths
	std::string Operator::addDecoder(int wIn, int wOut, int count){
				
		return reHelper->addDecoder(wIn, wOut, count);
	}
	
	std::string Operator::addArithOp(int width, int nrInputs, int count){
				
		return reHelper->addArithOp(width, nrInputs, count);
	}
	
	//TODO: find a better approximation for the resources
	//		currently just logic corresponding to the multiplexers
	//TODO: find a better approximation for the resources
	//		currently just logic corresponding to the state register
	//TODO: find a better approximation for the resources
	//		for now, RAM blocks are not used
	std::string Operator::addFSM(int nrStates, int nrTransitions, int count){
				
		return reHelper->addFSM(nrStates, nrTransitions, count);
	}
	
	//--Resource usage statistics---------------------------------------
	std::string Operator::generateStatistics(int detailLevel){
				
		return reHelper->generateStatistics(detailLevel);
	}
	
	//--Utility functions related to the generation of resource usage statistics
	
	//TODO: find a more precise way to determine the required number of
	//		registers due to pipeline
	std::string Operator::addPipelineFF(){
				
		return reHelper->addPipelineFF();
	}
	
	std::string Operator::addWireCount(){
				
		return reHelper->addWireCount();
	}
	
	std::string Operator::addPortCount(){
				
		return reHelper->addPortCount();
	}
	
	//TODO: add function to add resource count from specified component
	std::string Operator::addComponentResourceCount(){
				
		return reHelper->addComponentResourceCount();
	}
	
	void Operator::addAutomaticResourceEstimations(){
				
		resourceEstimate << reHelper->addAutomaticResourceEstimations();
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////
	
	
	/////////////////////////////////////////////////////////////////////////////////////////////////
	////////////Functions used for floorplanning
	/*void Operator::initFloorplanning(double ratio){
		
		floorplan << "";
		
		floorplanningRatio = ratio;
		
		virtualModuleId = 0;
		
		prevEstimatedCountFF 		 = 0;
		prevEstimatedCountLUT 		 = 0;
		prevEstimatedCountMemory 	 = 0;
		prevEstimatedCountMultiplier = 0;
	}
	
	std::string Operator::manageFloorplan(){
		//create a new virtual module for the resources that are between
		//	two modules, and add it to the corresponding lists
		std::ostringstream result;
		Operator* virtualComponent = new Operator(target_); 
		std::string moduleName = join("virtual_module_", virtualModuleId);
		
		result << "";
		
		virtualModuleId++;
		
		virtualComponent->setName(moduleName);
		virtualComponent->initResourceEstimation();
		virtualComponent->estimatedCountLUT 		= estimatedCountLUT-prevEstimatedCountLUT;
		virtualComponent->estimatedCountFF  		= estimatedCountFF-prevEstimatedCountFF;
		virtualComponent->estimatedCountMultiplier 	= estimatedCountMultiplier-prevEstimatedCountMultiplier;
		virtualComponent->estimatedCountMemory 		= estimatedCountMemory-prevEstimatedCountMemory;
		
		if((virtualComponent->estimatedCountLUT!=0) || (virtualComponent->estimatedCountFF!=0) 
				|| (virtualComponent->estimatedCountMultiplier!=0) || (virtualComponent->estimatedCountMemory!=0)){
			flComponentList.push_back(moduleName);
			flVirtualComponentList[moduleName] = virtualComponent;
			
			result << "Created virtual module - " << moduleName << endl;
			result << tab << "Added " << estimatedCountLUT-prevEstimatedCountLUT << " function generators" << endl;
			result << tab << "Added " << estimatedCountFF-prevEstimatedCountFF << " registers" << endl;
			result << tab << "Added " << estimatedCountMultiplier-prevEstimatedCountMultiplier << " multipliers" << endl;
			result << tab << "Added " << estimatedCountMemory-prevEstimatedCountMemory << " memories" << endl;
			
			prevEstimatedCountLUT = estimatedCountLUT;
			prevEstimatedCountFF = estimatedCountFF;
			prevEstimatedCountMultiplier = estimatedCountMultiplier;
			prevEstimatedCountMemory = estimatedCountMemory;
		}
		
		return result.str();
	}
	
	std::string Operator::addPlacementConstraint(std::string source, std::string sink, int type){
		std::ostringstream result;
		constraintType newConstraint;
		
		//check to see if the type of the constraint is valid
		if(!((type==TO_LEFT_OF) || (type==TO_RIGHT_OF) || (type==ABOVE) || (type==UNDER) || 
				(type==TO_LEFT_OF_WITH_EXTRA) || (type==TO_RIGHT_OF_WITH_EXTRA) || (type==ABOVE_WITH_EXTRA) || (type==UNDER_WITH_EXTRA))){
			cerr << "Error: Trying to add a placement constraint of undefined type." << endl;
			exit(1);
		}
		//check if the source is a valid sub-component name
		if(subComponents_.find(source)==subComponents_.end()){
			cerr << "Error: source sub-component " << source << " was not found" << endl;
			exit(1);
		}
		//check if the sink is a valid sub-component name
		if(subComponents_.find(sink)==subComponents_.end()){
			cerr << "Error: sink sub-component " << sink << " was not found" << endl;
			exit(1);
		}
		
		newConstraint.type = PLACEMENT;
		newConstraint.source = source;
		newConstraint.sink = sink;
		newConstraint.value = type;
		
		flConstraintList.push_back(newConstraint);
		
		result << "Created new placement constraint" << endl;
		result << tab << "Sink module " << sink
				<< " is " << ((type==TO_LEFT_OF) ? "to the left of" : (type==TO_RIGHT_OF) ? "to the right of" : (type==ABOVE) ? "above" : (type==UNDER) ? "under" : (type==TO_LEFT_OF_WITH_EXTRA) ? "to the left of" : (type==TO_RIGHT_OF_WITH_EXTRA) ? "to the right of" : (type==ABOVE_WITH_EXTRA) ? "above" : "under")
				<< " source module " << source << endl;
		
		return result.str();
	}
	
	std::string Operator::addConnectivityConstraint(std::string source, std::string sink, int nrWires){
		std::ostringstream result;
		constraintType newConstraint;
		
		//no non-positive values allowed for the number of wires
		if(nrWires<1){
			cerr << "Error: trying to add an invalid number of wires:" << nrWires << endl;
			exit(1);
		}
		//check if the source is a valid sub-component name
		if(subComponents_.find(source)==subComponents_.end()){
			cerr << "Error: source sub-component " << source << " was not found" << endl;
			exit(1);
		}
		//check if the sink is a valid sub-component name
		if(subComponents_.find(sink)==subComponents_.end()){
			cerr << "Error: sink sub-component " << sink << " was not found" << endl;
			exit(1);
		}
		
		newConstraint.type = CONNECTIVITY;
		newConstraint.source = source;
		newConstraint.sink = sink;
		newConstraint.value = nrWires;
		
		flConstraintList.push_back(newConstraint);
		
		result << "Created new connectivity constraint" << endl;
		result << tab << "Source module " << source 
				<< "is connected to sink module " << sink << " by " << nrWires << " wires" << endl;
		
		return result.str();
	}
	
	std::string Operator::processPlacementConstraints(){
		std::ostringstream result;
		
		//process each constraint, one by one, and update the virtual grid
		for(unsigned int i=0; i<flConstraintList.size(); i++){
			constraintType newConstraint = flConstraintList[i];
			coordinateType sourceCoordinate, sinkCoordinate, newCoordinate;
			bool positionOccupied = false;
			
			result << "Processed placement constraint" << endl;
			
			//get coordinates of source and sink modules
			sourceCoordinate = flComponentCordVirtual[newConstraint.source];
			sinkCoordinate = flComponentCordVirtual[newConstraint.sink];
			if(newConstraint.type==PLACEMENT){
				if(newConstraint.value==TO_LEFT_OF){//place sink to the left of source
					//new coordinates of sink module
					newCoordinate.x = sourceCoordinate.x - 1;
					newCoordinate.y = sourceCoordinate.y;
					
					if(sinkCoordinate.x<sourceCoordinate.x && sinkCoordinate.y==sourceCoordinate.y)
						continue;
					
					//as long as they are taken, go to the left of source
					do{
						int iterationsPerformed = 0;
						
						for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
							coordinateType tempCoordinate = it->second;
							
							if((tempCoordinate.x==newCoordinate.x) && (tempCoordinate.y==newCoordinate.y)){
								positionOccupied = true;
								break;
							}
							iterationsPerformed++;
						}
						if(iterationsPerformed == flComponentCordVirtual.size()){
							positionOccupied = false;
						}
						if(positionOccupied == true){
							newCoordinate.x -= 1;
						}
					}while(positionOccupied == true);
					
					//update the coordinates in the virtual grid
					flComponentCordVirtual[newConstraint.sink] = newCoordinate;
					
					result << tab << "Updated coordinates of module " << newConstraint.sink << endl;
					result << tab << tab << " from x=" << sinkCoordinate.x << " and y=" << sinkCoordinate.y << endl;
					result << tab << tab << " moving to the left of module " << newConstraint.source << endl;
					result << tab << tab << " to x=" << newCoordinate.x << " and y=" << newCoordinate.y << endl;
				}else if(newConstraint.value==TO_RIGHT_OF){//place sink to the left of source
					//new coordinates of sink module
					newCoordinate.x = sourceCoordinate.x + 1;
					newCoordinate.y = sourceCoordinate.y;
					
					if(sinkCoordinate.x>sourceCoordinate.x && sinkCoordinate.y==sourceCoordinate.y)
						continue;
					
					//as long as they are taken, go to the right of source
					do{
						int iterationsPerformed = 0;
						
						for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
							coordinateType tempCoordinate = it->second;
							
							if((tempCoordinate.x==newCoordinate.x) && (tempCoordinate.y==newCoordinate.y)){
								positionOccupied = true;
								break;
							}
							iterationsPerformed++;
						}
						if(iterationsPerformed == flComponentCordVirtual.size()){
							positionOccupied = false;
						}
						if(positionOccupied == true){
							newCoordinate.x += 1;
						}
					}while(positionOccupied == true);
					
					//update the coordinates in the virtual grid
					flComponentCordVirtual[newConstraint.sink] = newCoordinate;
					
					result << tab << "Updated coordinates of module " << newConstraint.sink << endl;
					result << tab << tab << " from x=" << sinkCoordinate.x << " and y=" << sinkCoordinate.y << endl;
					result << tab << tab << " moving to the right of module " << newConstraint.source << endl;
					result << tab << tab << " to x=" << newCoordinate.x << " and y=" << newCoordinate.y << endl;
				}else if(newConstraint.value==ABOVE){//place sink above source
					//new coordinates of sink module
					newCoordinate.x = sourceCoordinate.x;
					newCoordinate.y = sourceCoordinate.y - 1;
					
					if(sinkCoordinate.y<sourceCoordinate.y && sinkCoordinate.x==sourceCoordinate.x)
						continue;
					
					//as long as they are taken, go upwards of source
					do{
						int iterationsPerformed = 0;
						
						for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
							coordinateType tempCoordinate = it->second;
							
							if((tempCoordinate.x==newCoordinate.x) && (tempCoordinate.y==newCoordinate.y)){
								positionOccupied = true;
								break;
							}
							iterationsPerformed++;
						}
						if(iterationsPerformed == flComponentCordVirtual.size()){
							positionOccupied = false;
						}
						if(positionOccupied == true){
							newCoordinate.y -= 1;
						}
					}while(positionOccupied == true);
					
					//update the coordinates in the virtual grid
					flComponentCordVirtual[newConstraint.sink] = newCoordinate;
					
					result << tab << "Updated coordinates of module " << newConstraint.sink << endl;
					result << tab << tab << " from x=" << sinkCoordinate.x << " and y=" << sinkCoordinate.y << endl;
					result << tab << tab << " moving above module " << newConstraint.source << endl;
					result << tab << tab << " to x=" << newCoordinate.x << " and y=" << newCoordinate.y << endl;
				}else if(newConstraint.value==UNDER){//place sink under source
					//new coordinates of sink module
					newCoordinate.x = sourceCoordinate.x;
					newCoordinate.y = sourceCoordinate.y + 1;
					
					if(sinkCoordinate.y>sourceCoordinate.y && sinkCoordinate.x==sourceCoordinate.x)
						continue;
					
					//as long as they are taken, go downwards of source
					do{
						int iterationsPerformed = 0;
						
						for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
							coordinateType tempCoordinate = it->second;
							
							if((tempCoordinate.x==newCoordinate.x) && (tempCoordinate.y==newCoordinate.y)){
								positionOccupied = true;
								break;
							}
							iterationsPerformed++;
						}
						if(iterationsPerformed == flComponentCordVirtual.size()){
							positionOccupied = false;
						}
						if(positionOccupied == true){
							newCoordinate.y += 1;
						}
					}while(positionOccupied == true);
					
					//update the coordinates in the virtual grid
					flComponentCordVirtual[newConstraint.sink] = newCoordinate;
					
					result << tab << "Updated coordinates of module " << newConstraint.sink << endl;
					result << tab << tab << " from x=" << sinkCoordinate.x << " and y=" << sinkCoordinate.y << endl;
					result << tab << tab << " moving under module " << newConstraint.source << endl;
					result << tab << tab << " to x=" << newCoordinate.x << " and y=" << newCoordinate.y << endl;
				}else if(newConstraint.value==TO_LEFT_OF_WITH_EXTRA){
					//new coordinates of sink module
					newCoordinate.x = sourceCoordinate.x - 1;
					newCoordinate.y = sourceCoordinate.y;
					
					if(sinkCoordinate.x<sourceCoordinate.x && sinkCoordinate.y==sourceCoordinate.y)
						continue;
					
					//as long as they are taken, go downwards of source
					do{
						int iterationsPerformed = 0;
						
						for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
							coordinateType tempCoordinate = it->second;
							
							if((tempCoordinate.x==newCoordinate.x) && (tempCoordinate.y==newCoordinate.y)){
								positionOccupied = true;
								break;
							}
							iterationsPerformed++;
						}
						if(iterationsPerformed == flComponentCordVirtual.size()){
							positionOccupied = false;
						}
						if(positionOccupied == true){
							newCoordinate.x -= 1;
						}
					}while(positionOccupied == true);
					
					//look for the virtual module just before the sink module, if it exists
					// if it does, then place it between the source and the sink
					std::string prevModuleName = "";
					
					for(unsigned int i=0; i<flComponentList.size(); i++){
						if(newConstraint.sink == flComponentList[i])
							break;
						else
							prevModuleName = flComponentList[i];
					}
					if(prevModuleName.find("virtual_module_") != string::npos){
						flComponentCordVirtual[prevModuleName] = newCoordinate;
						newCoordinate.x -= 1;
					}
					
					//change the coordinates of the sink module
					flComponentCordVirtual[newConstraint.sink] = newCoordinate;
					
					result << tab << "Updated coordinates of module " << newConstraint.sink << endl;
					result << tab << tab << " from x=" << sinkCoordinate.x << " and y=" << sinkCoordinate.y << endl;
					result << tab << tab << " moving to the left of module " << newConstraint.source << endl;
					result << tab << tab << " to x=" << newCoordinate.x << " and y=" << newCoordinate.y << endl;
				}else if(newConstraint.value==TO_RIGHT_OF_WITH_EXTRA){
					//new coordinates of sink module
					newCoordinate.x = sourceCoordinate.x + 1;
					newCoordinate.y = sourceCoordinate.y;
					
					if(sinkCoordinate.x>sourceCoordinate.x && sinkCoordinate.y==sourceCoordinate.y)
						continue;
					
					//as long as they are taken, go downwards of source
					do{
						int iterationsPerformed = 0;
						
						for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
							coordinateType tempCoordinate = it->second;
							
							if((tempCoordinate.x==newCoordinate.x) && (tempCoordinate.y==newCoordinate.y)){
								positionOccupied = true;
								break;
							}
							iterationsPerformed++;
						}
						if(iterationsPerformed == flComponentCordVirtual.size()){
							positionOccupied = false;
						}
						if(positionOccupied == true){
							newCoordinate.x += 1;
						}
					}while(positionOccupied == true);
					
					//look for the virtual module just before the sink module, if it exists
					// if it does, then place it between the source and the sink
					std::string prevModuleName = "";
					
					for(unsigned int i=0; i<flComponentList.size(); i++){
						if(newConstraint.sink == flComponentList[i])
							break;
						else
							prevModuleName = flComponentList[i];
					}
					if(prevModuleName.find("virtual_module_") != string::npos){
						flComponentCordVirtual[prevModuleName] = newCoordinate;
						newCoordinate.x += 1;
					}
					
					//change the coordinates of the sink module
					flComponentCordVirtual[newConstraint.sink] = newCoordinate;
					
					result << tab << "Updated coordinates of module " << newConstraint.sink << endl;
					result << tab << tab << " from x=" << sinkCoordinate.x << " and y=" << sinkCoordinate.y << endl;
					result << tab << tab << " moving to the right of module " << newConstraint.source << endl;
					result << tab << tab << " to x=" << newCoordinate.x << " and y=" << newCoordinate.y << endl;
				}else if(newConstraint.value==ABOVE_WITH_EXTRA){
					//new coordinates of sink module
					newCoordinate.x = sourceCoordinate.x;
					newCoordinate.y = sourceCoordinate.y - 1;
					
					if(sinkCoordinate.y<sourceCoordinate.y && sinkCoordinate.x==sourceCoordinate.x)
						continue;
					
					//as long as they are taken, go downwards of source
					do{
						int iterationsPerformed = 0;
						
						for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
							coordinateType tempCoordinate = it->second;
							
							if((tempCoordinate.x==newCoordinate.x) && (tempCoordinate.y==newCoordinate.y)){
								positionOccupied = true;
								break;
							}
							iterationsPerformed++;
						}
						if(iterationsPerformed == flComponentCordVirtual.size()){
							positionOccupied = false;
						}
						if(positionOccupied == true){
							newCoordinate.y -= 1;
						}
					}while(positionOccupied == true);
					
					//look for the virtual module just before the sink module, if it exists
					// if it does, then place it between the source and the sink
					std::string prevModuleName = "";
					
					for(unsigned int i=0; i<flComponentList.size(); i++){
						if(newConstraint.sink == flComponentList[i])
							break;
						else
							prevModuleName = flComponentList[i];
					}
					if(prevModuleName.find("virtual_module_") != string::npos){
						flComponentCordVirtual[prevModuleName] = newCoordinate;
						newCoordinate.y -= 1;
					}
					
					//change the coordinates of the sink module
					flComponentCordVirtual[newConstraint.sink] = newCoordinate;
					
					result << tab << "Updated coordinates of module " << newConstraint.sink << endl;
					result << tab << tab << " from x=" << sinkCoordinate.x << " and y=" << sinkCoordinate.y << endl;
					result << tab << tab << " moving above module " << newConstraint.source << endl;
					result << tab << tab << " to x=" << newCoordinate.x << " and y=" << newCoordinate.y << endl;
				}else if(newConstraint.value==UNDER_WITH_EXTRA){
					//new coordinates of sink module
					newCoordinate.x = sourceCoordinate.x;
					newCoordinate.y = sourceCoordinate.y + 1;
					
					if(sinkCoordinate.y>sourceCoordinate.y && sinkCoordinate.x==sourceCoordinate.x)
						continue;
					
					//as long as they are taken, go downwards of source
					do{
						int iterationsPerformed = 0;
						
						for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
							coordinateType tempCoordinate = it->second;
							
							if((tempCoordinate.x==newCoordinate.x) && (tempCoordinate.y==newCoordinate.y)){
								positionOccupied = true;
								break;
							}
							iterationsPerformed++;
						}
						if(iterationsPerformed == flComponentCordVirtual.size()){
							positionOccupied = false;
						}
						if(positionOccupied == true){
							newCoordinate.y += 1;
						}
					}while(positionOccupied == true);
					
					//look for the virtual module just before the sink module, if it exists
					// if it does, then place it between the source and the sink
					std::string prevModuleName = "";
					
					for(unsigned int i=0; i<flComponentList.size(); i++){
						if(newConstraint.sink == flComponentList[i])
							break;
						else
							prevModuleName = flComponentList[i];
					}
					if(prevModuleName.find("virtual_module_") != string::npos){
						flComponentCordVirtual[prevModuleName] = newCoordinate;
						newCoordinate.y += 1;
					}
					
					//change the coordinates of the sink module
					flComponentCordVirtual[newConstraint.sink] = newCoordinate;
					
					result << tab << "Updated coordinates of module " << newConstraint.sink << endl;
					result << tab << tab << " from x=" << sinkCoordinate.x << " and y=" << sinkCoordinate.y << endl;
					result << tab << tab << " moving under module " << newConstraint.source << endl;
					result << tab << tab << " to x=" << newCoordinate.x << " and y=" << newCoordinate.y << endl;
				}
			}
		}
		
		//re-normalize the virtual grid
		//the constraints might have created wholes in the grid, 
		//	leading to inefficiency and unused space
		
		//scan the components' y coordinates; if there are negative 
		//	coordinates, shift all coordinates downwards, so as to
		//	have all coordinates non-negative
		int minLine = 0;
		
		for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
			if((it->second).y < minLine)
				minLine = (it->second).y;
		}
		if(minLine < 0){
			for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
				(it->second).y += abs(minLine);
			}
		}
		
		//scan the components' y coordinates; if there are empty lines,
		//	move the rest of the lines up
		bool hasEmptyLines = false;
		int maxLine = 0;
		
		//compute the lowest line
		for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
			coordinateType tempCoordinate = it->second;
			if(tempCoordinate.y>maxLine)
				maxLine = tempCoordinate.y;
		}
		
		//as long as the closest line to the current is more than a line 
		//	away, move the line upwards so as to fill the gaps
		do{
			int minGapSize = flComponentCordVirtual.size(), minGapLine = 0;
			
			for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
				coordinateType tempCoordinate = it->second;
				coordinateType closestCoordinate = it->second;
				
				if(tempCoordinate.y == maxLine)
					continue;
				
				closestCoordinate.y +=flComponentCordVirtual.size();
				
				for(map<string, coordinateType>::iterator it2 = flComponentCordVirtual.begin(); it2 != flComponentCordVirtual.end(); it2++){
					coordinateType tempCoordinate2 = it2->second;
					
					if((tempCoordinate2.y > tempCoordinate.y) && (tempCoordinate2.y < closestCoordinate.y) && (it->first != it2->first))
						closestCoordinate = tempCoordinate2;
				}
				
				if((closestCoordinate.y-tempCoordinate.y < minGapSize) && (closestCoordinate.y-tempCoordinate.y > 1)){
					minGapSize = closestCoordinate.y-tempCoordinate.y;
					minGapLine = tempCoordinate.y;
				}
			}
			
			if(minGapSize>1 && minGapSize!=flComponentCordVirtual.size()){
				for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
					if((it->second).y >= minGapLine)
						flComponentCordVirtual[(it->first)].y -= (minGapSize-1);
				}
			}else{
				hasEmptyLines = false;
			}
		}while(hasEmptyLines == true);
		
		//scan the components' x coordinates; if there are negative 
		//	coordinates, shift all coordinates to the right, so as to
		//	have all coordinates non-negative
		int minColumn = 0;
		
		for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
			if((it->second).x < minColumn)
				minColumn = (it->second).x;
		}
		if(minColumn < 0){
			for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
				(it->second).x += abs(minColumn);
			}
		}
		
		return result.str();
	}
	
	std::string Operator::processConnectivityConstraints(){
		//
		 // currently, as the user decides where each module goes, without
		 // the automation of the process, the relevance of the number of
		 // links between two modules is questionable.
		 // futrher modifications to follow.
		 //
		 return "";
	}
	
	std::string Operator::createVirtualGrid(){
		std::ostringstream result;
		
		result << "Created virtual grid of sub-components" << endl;
		
		//create the virtual component grid, and place the modules one
		//	under the other, in the order that they are instantiated
		//virtual modules, for the glue logic, are also considered as
		//	modules
		for(unsigned int i=0; i<flComponentList.size(); i++){
			std::string componentName = flComponentList[i];
			coordinateType newCoordinate;
			
			newCoordinate.x = 0;
			newCoordinate.y = i;
			
			flComponentCordVirtual[componentName] = newCoordinate;
			
			result << tab << "Sub-component " << componentName << " placed by default at x=0 and y=" << i << endl;
		}
		
		return result.str();
	}
	
	std::string Operator::createPlacementGrid(){
		ostringstream result;
		
		result << "Created real grid of components" << endl;
		
		//different processing techniques and efforst for modules that
		//	contain and that don't contain DSP and RAM blocks
		if(estimatedCountMultiplier!=0 || estimatedCountMemory!=0){
			
			result << tab << "Creating the grid of real coordintes with DSP or/and RAM requirements" << endl;
			
			int nrLines = 0, nrColumns = 0;
			int maxDSPperLine = 0, maxDSPperColumn = 0, maxRAMperLine = 0, maxRAMperColumn = 0;
			vector< vector<string> > subcomponentMatrixLines, subcomponentMatrixColumns;
			bool invertAxes = false, breakLines = false;
			int virtualDSPperColumn, virtualDSPperRow, virtualRAMperColumn, virtualRAMperRow;
			int maxComponentHeight = 0;
			
			//determine the maximum number of DSPs in one line of the floorplan
			//	if the nb. is larger the number of lines in the chip, throw warning to the user,
			//	then check to see the maximum number of DSPs in one column
			//-> if viable, invert axes and continue with the floorplan
			//		if not, break the lines longer than the capacity
			
			result << tab << "Created a (almost) 2D version of the virtual component grid, organized by lines" << endl;
			
			//construct a (almost) 2D version of the virtual grid, organized by lines
			for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
				if((it->second).y > nrLines)
					nrLines = (it->second).y;
			}
			for(int i=0; i<=nrLines; i++){
				vector<string> tempLevelList;
				int nrDSP = 0, nrRAM = 0;
				
				for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
					if((it->second).y == i){
						tempLevelList.push_back(it->first);
						
						Operator* tempOperator = subComponents_[it->first];
						nrDSP += (tempOperator->estimatedCountMultiplier > 0) 	? 1 : 0;
						nrRAM += (tempOperator->estimatedCountMemory > 0) 		? 1 : 0;
					}
				}
				subcomponentMatrixLines.push_back(tempLevelList);
				
				if(nrDSP > maxDSPperLine)
					maxDSPperLine = nrDSP;
				if(nrRAM > maxRAMperLine)
					maxRAMperLine = nrRAM;
			}
			
			result << tab << "Created a (almost) 2D version of the virtual component grid, organized by columns" << endl;
			
			//construct a (almost) 2D version of the virtual grid, organized by columns
			for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
				if((it->second).y > nrLines)
					nrColumns = (it->second).x;
			}
			for(int i=0; i<=nrColumns; i++){
				vector<string> tempLevelList;
				int nrDSP = 0, nrRAM = 0;
				
				for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
					if((it->second).x == i){
						tempLevelList.push_back(it->first);
						
						Operator* tempOperator = subComponents_[it->first];
						nrDSP += tempOperator->estimatedCountMultiplier;
						nrRAM += tempOperator->estimatedCountMemory;
					}
				}
				subcomponentMatrixColumns.push_back(tempLevelList);
				
				if(nrDSP > maxDSPperColumn)
					maxDSPperColumn = nrDSP;
				if(nrRAM > maxRAMperColumn)
					maxRAMperColumn = nrRAM;
			}
			
			result << tab << "Testing if columns and lines should/shouldn't be inverted" << endl;
			
			//determine if the axes should be inverted and if the lines
			//	should be broken in two
			if(maxDSPperLine > (target_->multiplierPosition).size()){
				if(maxDSPperColumn < (target_->multiplierPosition).size())
					invertAxes = true;
				else
					breakLines = true;
			}
			if(maxRAMperLine > (target_->memoryPosition).size()){
				if(maxRAMperColumn < (target_->memoryPosition).size())
					invertAxes = true;
				else
					breakLines = true;
			}
			
			//feature not yet finifshed or used
			//if(invertAxes == true){
				//for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
					//int temp = (it->second).x;
					//(it->second).x = (it->second).y;
					//(it->second).y = temp;
				//}
				
				//vector< vector<string> > tempMatrix;
				//tempMatrix = subcomponentMatrixColumns;
				//subcomponentMatrixColumns = subcomponentMatrixLines;
				//subcomponentMatrixLines = tempMatrix;
				
				//virtualDSPperColumn = (target_->multiplierPosition).size();
				//virtualDSPperRow = target_->dspPerColumn;
				
				//virtualRAMperColumn = (target_->memoryPosition).size();
				//virtualRAMperRow = target_->ramPerColumn;
				
				//int temp = maxDSPperColumn;
				//maxDSPperColumn = maxDSPperLine;
				//maxDSPperLine = temp;
				//temp = maxRAMperColumn;
				//maxRAMperColumn = maxRAMperLine;
				//maxRAMperLine = temp;
			//}else{
				//virtualDSPperColumn = target_->dspPerColumn;
				//virtualDSPperRow = (target_->multiplierPosition).size();
				
				//virtualRAMperColumn = target_->ramPerColumn;
				//virtualRAMperRow = (target_->memoryPosition).size();
			//}
			//if(breakLines == true){
				////should break the lines, but this may lead to considerable
				////	differences from the original floorplan
				////  so an Error is signaled to the user to reconsider the strategy
				//cerr << "Warning: the current floorplanning strategy, on the target FPGA, cannot be performed." 
						//<< " Please reconsider. No floorplan has been generated." << endl;
			//}
			//
			
			//determine the maximum height of a bounding box for a 
			//	sub-component based on the largets chain of DSPs or RAMs 
			//	in any of the sub-components
			for(map<string, Operator*>::iterator it = subComponents_.begin(); it != subComponents_.end(); it++){
				Operator* tempOperator = it->second;
				
				if((tempOperator->estimatedCountMultiplier/target_->dspHeightInLUT) > maxComponentHeight)
					maxComponentHeight = tempOperator->estimatedCountMultiplier/target_->dspHeightInLUT;
				if((tempOperator->estimatedCountMemory/target_->ramHeightInLUT) > maxComponentHeight)
					maxComponentHeight = tempOperator->estimatedCountMemory/target_->ramHeightInLUT;
			}
			
			result << tab << "Ordering the elements of the lines in the order that they appear in the virtual grid" << endl;
			
			//order the lines of the 2D virtual placement matrix, in
			//	increasing order of their x coordinate
			for(unsigned int i=0; i<subcomponentMatrixLines.size(); i++){
				vector<string> tempList = subcomponentMatrixLines[i];
				
				for(unsigned int j=0; j<tempList.size(); j++)
					for(unsigned int k=0; k<tempList.size(); k++){
						if((flComponentCordVirtual[tempList[j]]).x > (flComponentCordVirtual[tempList[k]]).x){
							string tempString = tempList[j];
							tempList[j] = tempList[k];
							tempList[k] = tempString;
						}
					}
					
				subcomponentMatrixLines[i] = tempList;
			}
			
			result << tab << "Creating the grid of real coordintes" << endl;
			
			//go through the 2D virtual placement list and generate the
			//	real coordinate list
			int prevCoordX, prevCoordY;
			
			result << tab << "Creating the initial placement in the real coordinate grid" << endl;
			
			//create an initial placement, then move the components around
			//	so as to satisfy their resource requirements
			for(unsigned int i=0; i<subcomponentMatrixLines.size(); i++){
				vector<string> matrixLine = subcomponentMatrixLines[i];
				
				prevCoordX = 0;
				prevCoordY = i * maxComponentHeight * (1.0/floorplanningRatio);
				for(unsigned int j=0; j<matrixLine.size(); j++){
					int lutWidth, ffWidth, componentWidth;
					coordinateType componentLocation, componentDimension;
					Operator* currentComponent = subComponents_[matrixLine[j]];
					
					lutWidth = (currentComponent->estimatedCountLUT/target_->lutPerSlice)/maxComponentHeight;		//width in slices
					ffWidth = (currentComponent->estimatedCountFF/target_->ffPerSlice)/maxComponentHeight;		//width in slices
					componentWidth = (lutWidth>ffWidth) ? lutWidth : ffWidth;
					
					componentLocation.x = prevCoordX + 1;
					componentLocation.y = prevCoordY;
					componentDimension.x = componentWidth * (1.0/floorplanningRatio);
					componentDimension.y = maxComponentHeight * (1.0/floorplanningRatio);
					
					flComponentCordReal[matrixLine[j]] = componentLocation;
					flComponentDimension[matrixLine[j]] = componentDimension;
					
					prevCoordX = prevCoordX + 1 + componentDimension.x;
				}
			}
			
			result << tab << "Adjusting the component grid so that components with RAM and DSP requirement are properly placed" << endl;
			
			//shift the components to the left so that they meet the
			//	DSP and RAM requirements
			for(unsigned int i=0; i<subcomponentMatrixLines.size(); i++){
				vector<string> matrixLine = subcomponentMatrixLines[i];
				
				for(unsigned int j=0; j<matrixLine.size(); j++){
					bool dspSatisfied = false, ramSatisfied = false;
						
					for(unsigned int k=0; k<target_->multiplierPosition.size(); k++){
						if(target_->multiplierPosition[k] == flComponentCordReal[matrixLine[j]].x+1){
							dspSatisfied = true;
							break;
						}
					}
					for(unsigned int k=0; k<target_->memoryPosition.size(); k++){
						if(target_->memoryPosition[k] == flComponentCordReal[matrixLine[j]].x+1){
							dspSatisfied = true;
							break;
						}
					}
					
					if(dspSatisfied && ramSatisfied)
						continue;
					else if(!dspSatisfied && !ramSatisfied){
						int closestDSPColumn = flComponentCordReal[matrixLine[j]].x, closestRAMColumn = flComponentCordReal[matrixLine[j]].x;
						int closestMove, dimensionIncrease;
						
						//decide which is closest, DSP or RAM, then shift to that column, 
						//	and then extend the width to accomodate the other resource
						for(unsigned int k=0; k<target_->multiplierPosition.size(); k++){
							if(target_->multiplierPosition[k] >= flComponentCordReal[matrixLine[j]].x){
								closestDSPColumn = target_->multiplierPosition[k]+1;
								break;
							}
						}
						for(unsigned int k=0; k<target_->memoryPosition.size(); k++){
							if(target_->memoryPosition[k] >= flComponentCordReal[matrixLine[j]].x){
								closestRAMColumn = target_->memoryPosition[k]+1;
								break;
							}
						}
						
						closestMove = (closestDSPColumn>closestRAMColumn) ? closestDSPColumn : closestRAMColumn;
						dimensionIncrease += abs(closestDSPColumn - closestRAMColumn);
						//if needed, shift all the components to the left
						if((closestMove - flComponentCordReal[matrixLine[j]].x) > 0){
							for(unsigned int k=j; k<matrixLine.size(); k++){
								flComponentCordReal[matrixLine[k]].x += dimensionIncrease + 
																			(closestMove - flComponentCordReal[matrixLine[j]].x);
							}
						}
					}else if(!dspSatisfied){
						int closestDSPColumn = flComponentCordReal[matrixLine[j]].x;
						
						for(unsigned int k=0; k<target_->multiplierPosition.size(); k++){
							if(target_->multiplierPosition[k] >= flComponentCordReal[matrixLine[j]].x){
								closestDSPColumn = target_->multiplierPosition[k]+1;
								break;
							}
						}
						
						//if needed, shift all the components to the left
						int shiftSize = closestDSPColumn - flComponentCordReal[matrixLine[j]].x;
						if(shiftSize > 0){
							for(unsigned int k=j; k<matrixLine.size(); k++){
								flComponentCordReal[matrixLine[k]].x += shiftSize;
							}
						}
					}else if(!ramSatisfied){
						int closestRAMColumn = flComponentCordReal[matrixLine[j]].x;
						
						for(unsigned int k=0; k<target_->memoryPosition.size(); k++){
							if(target_->memoryPosition[k] >= flComponentCordReal[matrixLine[j]].x){
								closestRAMColumn = target_->memoryPosition[k]+1;
								break;
							}
						}
						
						//if needed, shift all the components to the left
						int shiftSize = closestRAMColumn - flComponentCordReal[matrixLine[j]].x;
						if(shiftSize > 0){
							for(unsigned int k=j; k<matrixLine.size(); k++){
								flComponentCordReal[matrixLine[k]].x += shiftSize;
							}
						}
					}
				}
			}
			
		}else{
			//creating the placement for operators that don not have DSPs and RAMs
			result << tab << "Creating the grid of real coordintes" << endl;
			
			int nrLines = 0;
			vector< vector<string> > subcomponentMatrixLines;
			int maxComponentHeight = 0;
			
			result << tab << "Construct a (almost) 2D version of the virtual grid, organized by lines" << endl;
			
			//construct a (almost) 2D version of the virtual grid, organized by lines
			for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
				if((it->second).y > nrLines)
					nrLines = (it->second).y;
			}
			for(int i=0; i<=nrLines; i++){
				vector<string> tempLevelList;
				
				for(map<string, coordinateType>::iterator it = flComponentCordVirtual.begin(); it != flComponentCordVirtual.end(); it++){
					if((it->second).y == i){
						tempLevelList.push_back(it->first);
					}
				}
				subcomponentMatrixLines.push_back(tempLevelList);
			}
			
			//determine the maximum height of a bounding box of any of 
			//	the sub-components; the goal is to have balanced dimensions 
			//	for the boxes, so make the largest box as close to being 
			//	square -> hence the square root
			for(unsigned int i=0; i<flComponentList.size(); i++){
				Operator* tempOperator;
				int operatorWidth; 
				
				if(subComponents_.find(flComponentList[i]) != subComponents_.end()){
					tempOperator = subComponents_[flComponentList[i]];
				}else{
					tempOperator = flVirtualComponentList[flComponentList[i]];
				}
				
				operatorWidth = ceil(sqrt((tempOperator->estimatedCountLUT)/(double)target_->lutPerSlice));
				if(operatorWidth>maxComponentHeight)
					maxComponentHeight = operatorWidth;
				operatorWidth = ceil(sqrt((tempOperator->estimatedCountFF)/(double)target_->ffPerSlice));
				if(operatorWidth>maxComponentHeight)
					maxComponentHeight = operatorWidth;
			}
			
			//if the modules chain vertically and they run out of space, expand design horizontally
			while(ceil(nrLines*maxComponentHeight*sqrt(1.0/floorplanningRatio)) > target_->topSliceY){
				maxComponentHeight--;
				
				//design cannot be floorplanned with the current constraints
				if(maxComponentHeight == 0){
					cerr << "Error: the design cannot be floorplanned with the current constraints. Please reconsider and re-run." << endl;
					cerr << tab << "number of lines= " << nrLines << " maximum component height=" << maxComponentHeight << " ratio=" << 1.0/floorplanningRatio << endl;
					cerr << tab << "maximum allowed height=" << target_->topSliceY << endl;
					exit(1);
				}
			}
			
			result << tab << "Order the lines of the 2D virtual placement matrix" << endl;
			
			//order the lines of the 2D virtual placement matrix, in the
			//	increasing order of their x coordinates
			for(unsigned int i=0; i<subcomponentMatrixLines.size(); i++){
				vector<string> tempList = subcomponentMatrixLines[i];
				
				for(unsigned int j=0; j<tempList.size()-1; j++)
					for(unsigned int k=j+1; k<tempList.size(); k++){
						if((flComponentCordVirtual[tempList[j]]).x > (flComponentCordVirtual[tempList[k]]).x){
							coordinateType tempCoord = flComponentCordVirtual[tempList[j]];
							flComponentCordVirtual[tempList[j]] = flComponentCordVirtual[tempList[k]];
							flComponentCordVirtual[tempList[k]] = tempCoord;
							
							string tempString = tempList[j];
							tempList[j] = tempList[k];
							tempList[k] = tempString;
						}
					}
					
				subcomponentMatrixLines[i] = tempList;
			}
			
			result << tab << "Generate the real coordinate grid" << endl;
			
			//go through the 2D virtual placement list and generate the
			//	real coordinate list
			int prevCoordX, prevCoordY;
			
			//place each line at a time; lines with lower y coordinates 
			//	are processed first, and inside the lines the elements 
			// are processed in the increasing order of their x coordinates
			for(unsigned int i=0; i<subcomponentMatrixLines.size(); i++){
				vector<string> matrixLine = subcomponentMatrixLines[i];
				
				prevCoordX = 0;
				prevCoordY = ceil(i * maxComponentHeight * sqrt(1.0/floorplanningRatio));
				for(unsigned int j=0; j<matrixLine.size(); j++){
					int lutWidth, ffWidth, componentWidth;
					coordinateType componentLocation, componentDimension;
					Operator* currentComponent;
					
					if(subComponents_.find(matrixLine[j]) != subComponents_.end()){
						currentComponent = subComponents_[matrixLine[j]];
						
						lutWidth = ceil(((double)(currentComponent->estimatedCountLUT)/(target_->lutPerSlice))/maxComponentHeight);		//width in slices
						ffWidth = ceil(((double)(currentComponent->estimatedCountFF)/(target_->ffPerSlice))/maxComponentHeight);		//width in slices
						componentWidth = (lutWidth>ffWidth) ? lutWidth : ffWidth;
						
						componentLocation.x = prevCoordX + 1;
						componentLocation.y = prevCoordY;
						componentDimension.x = ceil(componentWidth * sqrt(1.0/(double)floorplanningRatio));
						componentDimension.y = ceil(maxComponentHeight * sqrt(1.0/(double)floorplanningRatio));
						
						flComponentCordReal[matrixLine[j]] = componentLocation;
						flComponentDimension[matrixLine[j]] = componentDimension;
						
						prevCoordX = prevCoordX + 1 + componentDimension.x;
					}
				}
			}
		}
		
		return result.str();
	}
	
	std::string Operator::createConstraintsFile(){
		ofstream file;
		ostringstream result;
		
		result << "Created output constraints file" << endl;
		
		//create the physical file
		if(target_->getVendor() == "Xilinx")
			file.open("flopoco.ucf");
		else if(target_->getVendor() == "Altera"){
		
		}
		
		result << tab << "Adding constraints" << endl;
		result << tab << "Added constraints to contain the entire operator" << endl;
		file << createPlacementForComponent("root");
		
		//create the placement constraints for each sub-component at a time
		//this time, the boxes for the glue logic aren't gien any bounds, 
		//	but they will have the empty space between the modules, which 
		//	should fit them
		for(unsigned int i=0; i<flComponentList.size(); i++){
			result << tab << "Added constraints for sub-component " << flComponentList[i] << endl;
			file << createPlacementForComponent(flComponentList[i]);
		}
		
		file.close();
		
		if(target_->getVendor() == "Xilinx"){
			cerr << "***Floorplan written to \'flopoco.ucf\' constraints file" << endl;
		}else if(target_->getVendor() == "Altera"){
			
		}
		
		return result.str();
	}
	
	std::string Operator::createPlacementForComponent(std::string moduleName){
		ostringstream constraintString;
		
		constraintString << "";
		
		//create the constraint for the whole operator
		if(moduleName == "root"){
			
			int maxX, maxY, minX, minY;
			map<string, coordinateType>::iterator it = flComponentCordReal.begin();
			
			minX = (it->second).x;
			minY = (it->second).y;
			
			maxX = minX + (flComponentDimension[it->first]).x;
			maxY = minY + (flComponentDimension[it->first]).y;
			
			it++;
			while(it != flComponentCordReal.end()){
				if((it->second).x < minX){
					minX = (it->second).x;
				}
				if((it->second).y < minY){
					minY = (it->second).y;
				}
				if((it->second).x+(flComponentDimension[it->first]).x > maxX){
					maxX = (it->second).x+(flComponentDimension[it->first]).x;
				}
				if((it->second).y+(flComponentDimension[it->first]).y > maxY){
					maxY = (it->second).y+(flComponentDimension[it->first]).y;
				}
				it++;
			}
			
			constraintString << "INST \"*\" AREA_GROUP=\"pblock_root\";" << endl;
			constraintString << "AREA_GROUP \"pblock_root\" RANGE=SLICE_X" 
				<< minX << "Y" << minY << ":SLICE_X" << maxX << "Y" << maxY << ";" << endl;
			constraintString << endl;
				
			return constraintString.str();
		}
		
		if(flVirtualComponentList.find(moduleName) != flVirtualComponentList.end())
			return constraintString.str();
		
		string instanceName = flInstanceNames[moduleName];
		
		if(target_->getVendor() == "Xilinx"){
			//create the constraint
			constraintString << "INST \"" << instanceName << "\" AREA_GROUP=\"pblock_" << instanceName << "\";" << endl;
			//add constraints for function generators and registers
			constraintString << "AREA_GROUP \"pblock_" << instanceName 
				<< "\" RANGE=SLICE_X" << (flComponentCordReal[moduleName]).x << "Y" << (flComponentCordReal[moduleName]).y
				<< ":SLICE_X" << (flComponentCordReal[moduleName]).x + (flComponentDimension[moduleName]).x
				<< "Y" << (flComponentCordReal[moduleName]).y + (flComponentDimension[moduleName]).y << ";" << endl;
			constraintString << "AREA_GROUP \"pblock_" << instanceName << "\" GROUP=OPEN;" << endl;
			constraintString << "AREA_GROUP \"pblock_" << instanceName << "\" PLACE=OPEN;" << endl;
			//add constraints for DSPs
			if((subComponents_[moduleName])->estimatedCountMultiplier != 0){
				vector<int> dspPositions;
				int dspInColumn = (flComponentDimension[moduleName]).y/(target_->dspHeightInLUT);
				
				for(unsigned int i=0; i<target_->multiplierPosition.size(); i++){
					int currentDSPColumn = target_->multiplierPosition[i];
					
					if(currentDSPColumn == (flComponentCordReal[moduleName]).x-1)
						dspPositions.push_back(currentDSPColumn);
					else if((currentDSPColumn>=(flComponentCordReal[moduleName]).x) && (currentDSPColumn<=(flComponentCordReal[moduleName]).x+(flComponentDimension[moduleName]).x))
						dspPositions.push_back(currentDSPColumn);
				}
				
				for(unsigned int i=0; i<dspPositions.size(); i++){
					int currentDSPColumn = dspPositions[i];
					
					constraintString << "AREA_GROUP \"pblock_" << flInstanceNames[moduleName] 
						<< "\" RANGE=DSP48_X" << currentDSPColumn << "Y" << (flComponentCordReal[moduleName]).y/target_->dspHeightInLUT
						<< ":DSP48_X" << currentDSPColumn
						<< "Y" << (flComponentCordReal[moduleName]).y/target_->dspHeightInLUT + dspInColumn << ";" << endl;
				}
			}
			//add constraints for RAMs
			if((subComponents_[moduleName])->estimatedCountMemory != 0){
				vector<int> ramPositions;
				int ramInColumn = (flComponentDimension[moduleName]).y/target_->dspHeightInLUT;
				
				for(unsigned int i=0; i<target_->multiplierPosition.size(); i++){
					int currentRAMColumn = target_->multiplierPosition[i];
					
					if(currentRAMColumn == (flComponentCordReal[moduleName]).x-1)
						ramPositions.push_back(currentRAMColumn);
					else if((currentRAMColumn>=(flComponentCordReal[moduleName]).x) && (currentRAMColumn<=(flComponentCordReal[moduleName]).x+(flComponentDimension[moduleName]).x))
						ramPositions.push_back(currentRAMColumn);
				}
				
				for(unsigned int i=0; i<ramPositions.size(); i++){
					int currentRAMColumn = ramPositions[i];
					
					if(target_->getID() == "Virtex4" || target_->getID() == "Spartan3"){
					constraintString << "AREA_GROUP \"pblock_" << flInstanceNames[moduleName] 
						<< "\" RANGE=RAMB16_X" << currentRAMColumn << "Y" << (flComponentCordReal[moduleName]).y/target_->ramHeightInLUT
						<< ":RAMB16_X" << currentRAMColumn
						<< "Y" << (flComponentCordReal[moduleName]).y/target_->ramHeightInLUT + ramInColumn << ";" << endl;
					}else{
					constraintString << "AREA_GROUP \"pblock_" << flInstanceNames[moduleName] 
						<< "\" RANGE=RAMB36_X" << currentRAMColumn << "Y" << (flComponentCordReal[moduleName]).y/target_->ramHeightInLUT
						<< ":RAMB36_X" << currentRAMColumn
						<< "Y" << (flComponentCordReal[moduleName]).y/target_->ramHeightInLUT + ramInColumn << ";" << endl;	
					}
				}
			}
			//end of constraint
			constraintString << endl;
		}else if(target_->getVendor() == "Altera"){
			//add constraints for function generators
			
			//add constraints for registers
			
			//add constraints for DSPs
			
			//add constraints for RAMs
			
		}
		
		return constraintString.str();
	}
	
	std::string Operator::createFloorplan(){
		ostringstream result;
		
		cerr << "=========================================================================" << endl;
		cerr << "*                          Floorplanning                                *" << endl;
		cerr << "=========================================================================" << endl;
		cerr << "Starting the creation of the floorplan for operator " << uniqueName_ << endl;
		cerr << "***Triggered creation of the virtual arrangement of the sub-components" << endl;
		result << createVirtualGrid();
		cerr << "***Triggered processing of placement constraints" << endl;
		result << processPlacementConstraints();
		cerr << "***Triggered processing of connectivity constraints" << endl;
		result << processConnectivityConstraints();
		cerr << "***Triggered creation of the actual arrangement of the sub-components" << endl;
		result << createPlacementGrid();
		cerr << "***Triggered creation of the constraints file" << endl;
		result << createConstraintsFile();
		cerr << "Finished creating the floorplan for operator " << uniqueName_ << endl;
		cerr << "=========================================================================" << endl;
		
		return result.str();
	}*/
	/////////////////////////////////////////////////////////////////////////////////////////////////
}


