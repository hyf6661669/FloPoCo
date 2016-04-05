// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "IntAdderTree.hpp"
#include "BitHeap/BitHeap.hpp"
#include "IntAddSubCmp/TargetOptCompressor.hpp"
//#include "PrimitiveComponents/Xilinx/Xilinx_TernaryAdd_2State.hpp"

using namespace std;
namespace flopoco {




	IntAdderTree::IntAdderTree(Target* target, int wIn, int noOfInputs, string method) : Operator(target), wIn_(wIn), noOfInputs_(noOfInputs) {
		/* constructor of the IntAdderTree
		   Target is the targeted FPGA : Stratix, Virtex ... (see Target.hpp for more informations)
		   param0 and param1 are some parameters declared by this Operator developpers, 
		   any number can be declared, you will have to modify 
		   -> this function,  
		   -> the prototype of this function (available in IntAdderTree.hpp)
		   -> the lines in main.cpp where the command line arguments are parsed in order to generate this IntAdderTree
		*/
		/* In this constructor we are going to generate an operator that takes as input three bit vectors X,Y,Z of lenght param0, treats them as unsigned integers, sums them and then output the last param1 bit of the sum adding the first bit of the sum (most significant) in front of this output, all the vhdl code needed by this operator has to be generated in this function */

		// definition of the source file name, used for info and error reporting using REPORT 
		srcFileName="IntAdderTree";

//		addHeaderVHDL("use work.xxxx.all");
//		addHeaderVHDL("use work.abs.all");

		// definition of the name of the operator
		ostringstream name;
		name << "IntAdderTree_" << wIn_ << "_" << noOfInputs_;
		setName(name.str());
		// Copyright 
		setCopyrightString("ACME and Co 2010");

		setSequential();
		useNumericStd();

		/* SET UP THE IO SIGNALS
		   Each IO signal is declared by addInput(name,n) or addOutput(name,n) 
		   where name is a string that stands for the name of the variable and 
		   n is an integer (int)   that stands for the length of the corresponding 
		   input/output */

		ostringstream bitHeapName;
		bitHeapName << "_" << target->getID();
		if(target->compressionType() == 0)
			bitHeapName << "_heuristic";
		else
			bitHeapName << "_optimal";

		for(int i=0; i < noOfInputs; i++)
		{
			addInput(join("X",i+1),wIn);
		}

		nextCycle(); //sync inputs (for timing analysis)

		REPORT(DEBUG, "Using compression method: " << method);
		if(method.compare("bitheap") == 0)
		{
			int wOut = ceil(log2(noOfInputs)) + wIn + 2; //2 guard bits are used as extra bits may occur due to inefficient mappings
			BitHeap* bitheap = new BitHeap(this, wOut, true, bitHeapName.str());

			for(int i=0; i < noOfInputs; i++)
			{
				vhdl << tab << declare(join("I",i+1),wIn) << " <= " << join("X",i+1) << ";" << endl;
                bitheap->addUnsignedBitVector(0,join("I",i+1),wIn);
			}

			bitheap->generateCompressorVHDL();

			addOutput("Y" , bitheap->getMaxWeight()+1);

			vhdl << tab << "Y <= " << bitheap->getSumName() << ";" << endl;
		}
		else if(method.compare("add2") == 0)
		{
			for(int i=0; i < noOfInputs; i++)
			{
				vhdl << tab << declare(join("X_0_",i+1),wIn) << " <= " << join("X",i+1) << ";" << endl;
			}
			int inputsAtStage = noOfInputs;

			int s=0;
			while(inputsAtStage > 1)
			{
				int j=1;
				for(int i=0; i < inputsAtStage; i+=2)
				{
					if(i < inputsAtStage-1)
					{
						vhdl 	<< tab << declare(join("X_",s+1,"_",j++),wIn+s+1) 
							<< " <= std_logic_vector(unsigned('0' & " << join("X_",s,"_",i+1) 
							<< ") + unsigned('0' & " << join("X_",s,"_",i+2) << "));" << endl;
					}
					else
					{
						vhdl << tab << declare(join("X_",s+1,"_",j++),wIn+s+1) << " <= '0' & " << join("X_",s,"_",i+1) << ";" << endl;
					}
				}
				inputsAtStage = ceil(inputsAtStage/2.0);
				s++;
				nextCycle();
			}
			addOutput("Y" , wIn+s);

			vhdl << tab << "Y <= X_" << s << "_1;" << endl;
		}
/*
		else if(method.compare("add3") == 0)
		{
			for(int i=0; i < noOfInputs; i++)
			{
				vhdl << tab << declare(join("X_0_",i+1),wIn) << " <= " << join("X",i+1) << ";" << endl;
			}
			int inputsAtStage = noOfInputs;

			int s=0;
			int wordSizeInStage = wIn;
			while(inputsAtStage > 1)
			{
				int j=1;
				for(int i=0; i < inputsAtStage; i+=3)
				{
					if(i < inputsAtStage-2)
					{
						Xilinx_TernaryAdd_2State* add3 = new Xilinx_TernaryAdd_2State(target,wIn+2*(s+1),0,-1);
						addSubComponent(add3);

						inPortMap( add3,"x_i",join("X_",s,"_",i+1) );
						inPortMap( add3,"y_i",join("X_",s,"_",i+2) );
						inPortMap( add3,"z_i",join("X_",s,"_",i+3) );
						outPortMap(add3,"sum_o", join("X_",s+1,"_",j++) );
						vhdl << instance(add3,join("add3_",s,"_",i)) << std::endl;
					}
					else if(i < inputsAtStage-1)
					{
						vhdl << tab << declare(join("X_",s+1,"_",j++),wIn+2*(s+1)) 
						<< " <= std_logic_vector(unsigned(\"00\" & " << join("X_",s,"_",i+1) 
						<< ") + unsigned(\"00\" & " << join("X_",s,"_",i+2) << "));" << endl;
					}
					else
					{
						vhdl << tab << declare(join("X_",s+1,"_",j++),wIn+2*(s+1)) << " <= \"00\" & " << join("X_",s,"_",i+1) << ";" << endl;
					}
				}
				inputsAtStage = ceil(inputsAtStage/3.0);
				s++;
				wordSizeInStage += 2;
				nextCycle();
			}

			addOutput("Y" , wIn+2*s); //currently, this is a worst-case estimate, so some of the MSBs are optimized during synthesis as they are zero

			vhdl << tab << "Y <= X_" << s << "_1;" << endl;
		}
*/
		else
		{
			THROWERROR("method " << method << " unknown!");
		}
		//nextCycle(); //sync output (for timing analysis)

	};

	OperatorPtr IntAdderTree::parseArguments(Target *target, vector<string> &args) {
		int wIn;
		UserInterface::parseInt(args, "wIn", &wIn);
		int noOfInputs;
		UserInterface::parseInt(args, "noOfInputs", &noOfInputs);
		return new IntAdderTree(target, wIn, noOfInputs, "bitheap");
	}


	void IntAdderTree::registerFactory(){
		UserInterface::add("IntAdderTree", // name
											 "An adder tree, built from 2-input adders, 3-input adders or bitheap compression",
											 "BasicInteger", // categories
											 "",
											 "wIn(int): word size of inputs;\
											  noOfInputs(int): no of inputs;",
											 "",
											 IntAdderTree::parseArguments
											 ) ;
	}	
	
	void IntAdderTree::emulate(TestCase * tc) {
		/* This function will be used when the TestBench command is used in the command line
		   we have to provide a complete and correct emulation of the operator, in order to compare correct output generated by this function with the test input generated by the vhdl code */
		/* first we are going to format the entries */

		mpz_class sy=0;
		for(int i=0; i < noOfInputs_; i++)
		{
			sy += tc->getInputValue(join("X",i+1));
		}

		tc->addExpectedOutput("Y",sy);
	}


	void IntAdderTree::buildStandardTestCases(TestCaseList * tcl) {
		// please fill me with regression tests or corner case tests!
	}
	
	
	
}//namespace
