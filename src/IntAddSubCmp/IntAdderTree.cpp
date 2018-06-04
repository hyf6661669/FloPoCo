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
#include "PrimitiveComponents/Xilinx/Xilinx_TernaryAdd_2State.hpp"

using namespace std;
namespace flopoco {




    IntAdderTree::IntAdderTree(Target* target, int wIn, int noOfInputs, string method, bool redundantOutput, bool useVariableColumnCompressors) : Operator(target), wIn_(wIn), noOfInputs_(noOfInputs) {
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
        name << "IntAdderTree_" << method << "_" << wIn_ << "_" << noOfInputs_;
		setName(name.str());
		// Copyright
        setCopyrightString("Marco Kleinlein, Martin Kumm");

		setSequential();
		useNumericStd();

		/* SET UP THE IO SIGNALS
		   Each IO signal is declared by addInput(name,n) or addOutput(name,n)
		   where name is a string that stands for the name of the variable and
		   n is an integer (int)   that stands for the length of the corresponding
		   input/output */

		ostringstream bitHeapName;
		bitHeapName << "_" << target->getID();

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

            bitheap->useVariableColumnCompressors = useVariableColumnCompressors;
            bitheap->generateCompressorVHDL(redundantOutput);

            nextCycle(); //sync output (for timing analysis)
            if(redundantOutput)
            {
                int outputwidth=getSignalByName(join(bitheap->getSumName(),0))->width();
                addOutput("YS" , outputwidth);
                outputwidth=getSignalByName(join(bitheap->getSumName(),1))->width();
                addOutput("YC" , outputwidth);

                vhdl << tab << "YS <= " << bitheap->getSumName() << "0;" << endl;
                vhdl << tab << "YC <= " << bitheap->getSumName() << "1";
                if(getSignalByName(join(bitheap->getSumName(),1))->width() > outputwidth)
                    vhdl << " & (others => ''0'')";
                vhdl << ";" << endl;
            }
            else
            {
                if(bitheap->getMaxHeight() > 1)
                {
                    vhdl << tab << declare("Y_raw", bitheap->getMaxWeight() + 1) << " <= " << bitheap->getSumName() << ";" << endl;
                }
                else
                {
                    vhdl << tab << declare("Y_raw", bitheap->getMaxWeight()) << " <= " << bitheap->getSumName() << ";" << endl;
                }
            }

		}
		else if(method.compare("add2") == 0)
		{
			for(int i=0; i < noOfInputs; i++)
			{
				vhdl << tab << declare(join("X_0_",i+1),wIn) << " <= " << join("X",i+1) << ";" << endl;
			}
			int inputsAtStage = noOfInputs;

			int s=0;
            while(inputsAtStage > 1+redundantOutput) //if redundantOutput=true, check is for >2
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
            vhdl << tab << declare("Y_raw", wIn+s) << " <= X_" << s << "_1;" << endl;
		}
        else if(method.compare("add3") == 0 )
        {
            if( !UserInterface::useTargetSpecificOptimization ){
                THROWERROR( "This component can't be used without target specific optimization." )
            }

            for(int i=0; i < noOfInputs; i++){
                vhdl << tab << declare(join("X_0_",i+1),wIn ) << " <= " << join("X",i+1) << ";" << endl;
            }
            int inputsAtStage = noOfInputs;

            int s=0;
            int wordSizeInStage = wIn;

            while(inputsAtStage > 1+redundantOutput) //if redundantOutput=true, check is for >2
            {
                int j=1;
                for(int i=0; i < inputsAtStage; i+=3)
                {
                    if(i < inputsAtStage-2)
                    {
                        Xilinx_TernaryAdd_2State* add3 = new Xilinx_TernaryAdd_2State(target,wIn+2*(s+1),0,-1);
                        addSubComponent(add3);

                        vhdl << tab << declare(join("X_",s,"_",i+1,"t"),wIn+2*(s+1)) << " <= \"00\" & " << join("X_",s,"_",i+1) << ";" << endl;
                        vhdl << tab << declare(join("X_",s,"_",i+2,"t"),wIn+2*(s+1)) << " <= \"00\" & " << join("X_",s,"_",i+2) << ";" << endl;
                        vhdl << tab << declare(join("X_",s,"_",i+3,"t"),wIn+2*(s+1)) << " <= \"00\" & " << join("X_",s,"_",i+3) << ";" << endl;

                        inPortMap( add3,"x_i",join("X_",s,"_",i+1,"t") );
                        inPortMap( add3,"y_i",join("X_",s,"_",i+2,"t") );
                        inPortMap( add3,"z_i",join("X_",s,"_",i+3,"t") );
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

            vhdl << declare("Y_raw", wIn + 2*s) << " <= X_" << s << "_1;" << endl; //currently, this is a worst-case estimate, so some of the MSBs are optimized during synthesis as they are zero
        }
		else
		{
			THROWERROR("method " << method << " unknown!");
		}

        if(method.compare("add3") == 0 || method.compare("add2") == 0 || method.compare("bitheap") == 0)
        {
            //compute wordsize of output
            mpz_t singleVectorMaxValue;
            mpz_t totalMaxValue;
            mpz_inits(singleVectorMaxValue, totalMaxValue, NULL);

            mpz_pow_ui(singleVectorMaxValue, mpz_class(2).get_mpz_t(), (unsigned) wIn);
            mpz_sub_ui(singleVectorMaxValue, singleVectorMaxValue, 1);  //max value of single input = (2^wIn) - 1
            mpz_mul_ui(totalMaxValue, singleVectorMaxValue, (unsigned) noOfInputs); // totalmaxValue = noOfInputs * ((2^wIn) - 1)
            unsigned int outputLength = mpz_sizeinbase(totalMaxValue, 2);
            //cout << "outputLength is " << outputLength << " and totalMaxValue is " << totalMaxValue << endl;
            mpz_clears(singleVectorMaxValue, totalMaxValue, NULL);

            //get the totalMaxValue lowest bits
            addOutput("Y", outputLength);
            vhdl << tab << "Y <= Y_raw(" << outputLength - 1 << " downto 0);" << endl;
        }

		//nextCycle(); //sync output (for timing analysis)

	};

	OperatorPtr IntAdderTree::parseArguments(Target *target, vector<string> &args) {
		int wIn;
		UserInterface::parseInt(args, "wIn", &wIn);
		int noOfInputs;
		UserInterface::parseInt(args, "noOfInputs", &noOfInputs);
		string method;
		UserInterface::parseString(args, "method", &method);

        bool redundantOutput;
        UserInterface::parseBoolean(args, "redundant", &redundantOutput);

        bool useVariableColumnCompressors;
        UserInterface::parseBoolean(args, "useVCC", &useVariableColumnCompressors);

        return new IntAdderTree(target, wIn, noOfInputs, method, redundantOutput, useVariableColumnCompressors);
	}


	void IntAdderTree::registerFactory(){
		UserInterface::add("IntAdderTree", // name
											 "An adder tree, built from 2-input adders, 3-input adders or bitheap compression",
											 "BasicInteger", // categories
											 "",
											 "wIn(int): word size of inputs;\
											  noOfInputs(int): no of inputs;\
                                              method(string): <bitheap,add2,add3>;\
                                              useVCC(bool)=false: if true, variable column compressors will be included in the bitheap optimization;\
                                              redundant(bool)=false: if true, the output will be in redundant (carry-save) format and the final adder will be omitted",
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
        TestCase *tc;

        tc = new TestCase(this);
        for(int i=0; i < noOfInputs_; i++)
        {
            tc->addInput(join("X",i+1), mpz_class(0));
        }
        emulate(tc);
        tc->addComment("Addition of zeros");
        tcl->add(tc);

        tc = new TestCase(this);
        for(int i=0; i < noOfInputs_; i++)
        {
            tc->addInput(join("X",i+1), mpz_class(1));
        }
        emulate(tc);
        tc->addComment("Addition of ones");
        tcl->add(tc);

        tc = new TestCase(this);
        for(int i=0; i < noOfInputs_; i++)
        {
            tc->addInput(join("X",i+1), ((mpz_class(1) << wIn_) - 1));
        }
        emulate(tc);
        tc->addComment("Addition of max positive value (unsigned)");   //???
        tcl->add(tc);

        tc = new TestCase(this);
        for(int i=0; i < noOfInputs_; i++)
        {
            tc->addInput(join("X",i+1), (mpz_class(1) << (wIn_ - 1)));
        }
        emulate(tc);
        tc->addComment("Addition of min negative value (signed)");   //???
        tcl->add(tc);
    }



}//namespace
