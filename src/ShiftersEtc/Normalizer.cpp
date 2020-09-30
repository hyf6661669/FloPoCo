/*
  A normalizer (leading zero/one counter + shifter + optional sticky bit) for FloPoCo

  Author: Florent de Dinechin

   This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,
  2008-2011.
  All rights reserved.

*/

#include <iostream>
#include <sstream>
#include <vector>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "Normalizer.hpp"

using namespace std;



namespace flopoco{


	Normalizer::Normalizer(OperatorPtr parentOp, Target* target, int wX, int wR, int wCount, bool computeSticky, const int countType) :
		Operator(parentOp, target), wX_(wX), wR_(wR), wCount_(wCount), computeSticky_(computeSticky), countType_(countType) {

		// -------- Parameter set up -----------------
		srcFileName = "Normalizer";
		setCopyrightString("Florent de Dinechin, (2007-2020)");

		REPORT(DETAILED, "wX="<<wX << " wR="<<wR << " wCount="<<wCount << " computeSticky=" << computeSticky  << " countType=" << countType);

		ostringstream name;
		name << "L" << (countType_<0?"ZO":((countType_>0)?"O":"Z")) << "CShifter"
			  << (computeSticky_?"Sticky":"") << "_" << wX_ << "_to_"<<wR_<<"_counting_"<<(1<<wCount_);
		setNameWithFreqAndUID(name.str());

		addInput ("X", wX_);
		if (countType_==-1) addInput ("OZb"); /* if we generate a generic LZOC */
		addOutput("Count", wCount_);
		addOutput("R", wR_);
		if (computeSticky_)   addOutput("Sticky"); /* if we require a sticky bit computation */

		// we consider that wR <= wX. We fix this at the end if not the case
		int wR_true = wR_;
		wR_ = wR > wX_ ? wX_ : wR;


		vhdl << tab << declare(join("level",wCount_), wX_) << " <= I ;"   <<endl;
		if (countType_==-1) vhdl << tab << declare("sozb") << "<= OZb;"<<endl;
		if ((computeSticky_)&&(wR_<wX))   vhdl << tab << declare(join("sticky",wCount_)) << " <= '0' ;"<<endl; //init sticky


		// Now comes the main loop.
		// i is the level index. Level i counts 2^i bits, and shifts by 2^i
		int currLevSize=wX, prevLevSize=0;
		for (int i=wCount_-1; i>=0; i--){
			prevLevSize = currLevSize;

			// level(k) = max ( max (2^k, wR) + 2^k -1) , wX)
			currLevSize = (wR_>intpow2(i)?wR_:intpow2(i));
			currLevSize += (intpow2(i)-1);
			currLevSize = (currLevSize > wX_? wX_: currLevSize);

			// Delay evaluation.
			// As we output the count bits, their computation will not be merged inside the shift
			//REPORT( DEBUG, "currSize="<<currLevSize);

			double countBitDelay = getTarget()->fanoutDelay(currLevSize);
			if (countType>=0)
				countBitDelay += getTarget()->eqConstComparatorDelay( intpow2(i) )  ;
			else
				countBitDelay += getTarget()->eqComparatorDelay( intpow2(i) ) ;
			
			vhdl << tab << declare(countBitDelay, join("count",i))
					 << "<= '1' when " <<join("level",i+1)<<range(prevLevSize-1,prevLevSize - intpow2(i))<<" = "
					 <<"("<<prevLevSize-1<<" downto "<<prevLevSize - intpow2(i)<<"=>"<< (countType_==-1? "sozb": countType_==0?"'0'":"'1'")<<") else '0';"<<endl;

			// The shift will take at most one LUT delay per level. We don't take into account that shift level can be merged: TODO ? It seems non-trivial.
			double shiftDelay = getTarget()->logicDelay(3);
			vhdl << tab << declare(shiftDelay,join("level",i),currLevSize)
					 << "<= " << join("level",i+1)<<"("<<prevLevSize-1<<" downto "<< prevLevSize-currLevSize << ")"
					 << " when " << join("count",i) << "='0' else ";
			int l,r;
			l = prevLevSize - intpow2(i) - 1;
			r = (currLevSize < prevLevSize - intpow2(i) ? (prevLevSize - intpow2(i)) - currLevSize : 0 );
			if (l>=r)
				vhdl << join("level",i+1) << "("<<prevLevSize - intpow2(i) - 1 <<" downto "<< (currLevSize < prevLevSize - intpow2(i) ? (prevLevSize - intpow2(i)) - currLevSize : 0 ) <<")";
			
			if (prevLevSize - intpow2(i) < currLevSize )
				vhdl << (l>=r?" & ":"") << rangeAssign(currLevSize -(prevLevSize - intpow2(i))-1,0,"'0'");
			vhdl << ";"<<endl;

			if ((computeSticky_)&&(wR_<wX)) {

				// Delay computation. Here we try to compute as much of the sticky in each level.
				double levelStickyDelay;
				// n is the size on which we compute the sticky bit
				int n = max( prevLevSize-currLevSize,
										 (currLevSize < prevLevSize - intpow2(i) ? (prevLevSize - int(intpow2(i)) ) - currLevSize : 0 ))  ;
				if ( countType_ == -1 )
					levelStickyDelay= getTarget()->eqComparatorDelay(n);
				else
					levelStickyDelay= getTarget()->eqConstComparatorDelay(n);
				

				vhdl << tab << declare(join("sticky_high_",i)) << "<= '0'";
				if (prevLevSize-currLevSize > 0)
					vhdl << "when " << join("level",i+1)<<"("<<prevLevSize-currLevSize -1 <<" downto "<< 0 <<") = CONV_STD_LOGIC_VECTOR(0,"<< prevLevSize-currLevSize <<") else '1'";
				vhdl << ";"<<endl;

   			vhdl << tab << declare(join("sticky_low_",i)) << "<= '0'";
				if ((currLevSize < prevLevSize - intpow2(i) ? (prevLevSize - intpow2(i)) - currLevSize : 0 ) > 0)
					vhdl << "when " <<join("level",i+1)<<"("<<(currLevSize < prevLevSize - intpow2(i) ? (prevLevSize - intpow2(i)) - currLevSize : 0 ) -1
						  <<" downto "<< 0 <<") = CONV_STD_LOGIC_VECTOR(0,"<< (currLevSize < prevLevSize - intpow2(i) ? (prevLevSize - intpow2(i)) - currLevSize : 0 ) <<") else '1'";
				vhdl << ";"<<endl;

				vhdl << tab << declare(levelStickyDelay, join("sticky",i))
						 << "<= " << join("sticky",i+1) << " or " << join("sticky_high_",i)
						 << " when " << join("count",i) << "='0' else " << join("sticky",i+1) << " or " << join("sticky_low_",i)<<";"<<endl;
			}

			vhdl <<endl;
		}

		//assign back the value to wR_
		wR_ =  wR_true;
		vhdl << tab << "O <= "<< join("level",0)
			  << (wR_<=wX?"":join("&",rangeAssign(wR_-wX-1,0,"'0'")))<<";"<<endl;


		vhdl << tab << declare("sCount",wCount_) <<(wCount_==1?"(0)":"")<<" <= ";
		for (int i=wCount_-1; i>=0; i--){
			vhdl <<join("count",i);
			vhdl << (i>0?" & ":join(";","\n"));
		}

		if((1<<wCount_)-1 > wX_) {
			vhdl << tab << "Count <= sCount;"<<endl;

			/*			vhdl << tab << "Count <= CONV_STD_LOGIC_VECTOR("<<wX_<<","<<wCount_<<") when sCount=CONV_STD_LOGIC_VECTOR("<<intpow2(wCount_)-1<<","<<wCount_<<")"<<endl
				  << tab << tab << "else sCount;"<<endl;*/
		}
		else {
			vhdl << tab << "Count <= sCount;"<<endl;
		}


		if (computeSticky_){
			if (wR_>=wX)
				vhdl << tab << "Sticky <= '0';"<<endl;
			else
				vhdl << tab << "Sticky <= sticky0;"<<endl;
		}
		REPORT( DEBUG, "Leaving Normalizer");

	}

	Normalizer::~Normalizer() {
	}


	int Normalizer::getCountWidth() const{
		return wCount_;
	}



	void Normalizer::emulate(TestCase* tc)
	{
		mpz_class inputValue  = tc->getInputValue("I");

		mpz_class sozb = 42; //dummy value
		if (countType_ == -1)
			sozb = tc->getInputValue("OZb");

		int sticky=0;
		int count =0;
		mpz_class shiftOutputValue = inputValue;


		mpz_class bit = (countType_ == -1) ? sozb : (countType_ == 0 ? 0 : 1); /* what are we counting in the specific case */

		int j=wX_-1;
		while ((count < (1<<wCount_)-1) &&  (j>=0)  && mpz_tstbit(inputValue.get_mpz_t(), j) == bit)   {
			count ++;
			j--;
			shiftOutputValue = shiftOutputValue <<1;
			}

		// Now reformat the output value to its size, and compute the sticky of the remaining bits.
		// The max size of shiftOutputValue is ((1<<wCount)-1) + wX
		mpz_class outputMask = (mpz_class(1) << wR_) -1;
		int numBitsForSticky = wX_ - wR_;
		if(numBitsForSticky >= 0) {// should be the typical use case where we need to compute a sticky bit
			mpz_class stickyMask = (mpz_class(1) << numBitsForSticky) -1;
			mpz_class bitsForSticky = shiftOutputValue & stickyMask;
			sticky = (bitsForSticky==0? 0 :1 );
			shiftOutputValue = (shiftOutputValue >> numBitsForSticky) & outputMask;
		}
		else {
			shiftOutputValue = (shiftOutputValue << numBitsForSticky) & outputMask;
			sticky=0;
		}
			
		
		tc->addExpectedOutput("O", shiftOutputValue);
		tc->addExpectedOutput("Count", count);

		if (computeSticky_)
			tc->addExpectedOutput("Sticky",sticky);
	}


	

	OperatorPtr Normalizer::parseArguments(OperatorPtr parentOp, Target *target, std::vector<std::string> &args) {
		int wX, wR, wCount, countType;
		bool computeSticky;
		UserInterface::parseStrictlyPositiveInt(args, "wX", &wX);
		UserInterface::parseStrictlyPositiveInt(args, "wR", &wR);
		UserInterface::parseStrictlyPositiveInt(args, "wCount", &wCount);
		UserInterface::parseBoolean(args, "computeSticky", &computeSticky);
		UserInterface::parseInt(args, "countType", &countType);
		return new Normalizer(parentOp, target, wX, wR, wCount, computeSticky, countType);
	}


	
	void Normalizer::registerFactory(){
		UserInterface::add("Normalizer", // name
											 "A combined leading zero/one counter and left shifter, useful for floating-point normalization.",
											 "ShiftersLZOCs",  // category
											 "", // see also
											 "wX(int): input size in bits;\
                        wR(int): output size in bits;\
                        wCount(int): size in bits of the count output;\
                        computeSticky(bool)=false: if true and wR<wX, a sticky bit is computed out of the discarded bits;\
                        countType(int)=-1:  0 to count zeroes, 1 to count ones, -1 to have a dynamic OZb input that tells what to count", // This string will be parsed
											 "", // no particular extra doc needed
											 Normalizer::parseArguments
											 ) ;
		
	}


}
