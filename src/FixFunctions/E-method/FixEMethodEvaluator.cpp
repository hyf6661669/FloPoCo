/*
 * FixEMethodEvaluator.cpp
 *
 *  Created on: 7 Dec 2017
 *      Author: mistoan
 */

#include "FixEMethodEvaluator.hpp"

namespace flopoco {

	FixEMethodEvaluator::FixEMethodEvaluator(Target* target, size_t _radix, size_t _maxDigit, int _msbInOut, int _lsbInOut,
		vector<string> _coeffsP, vector<string> _coeffsQ, map<string, double> inputDelays)
	: Operator(target), radix(_radix), maxDigit(_maxDigit),
	  	  n(_coeffsP.size()), m(coeffsQ.size()),
	  	  msbInOut(_msbInOut), lsbInOut(_lsbInOut),
		  maxDegree(n>m ? n : m)
	{
		ostringstream name;

		srcFileName = "FixEMethodEvaluator";
		name << "FixEMethodEvaluator_n_" << n << "_m_" << m
				<< "_msbInOut_" << vhdlize(msbInOut) << "_lsbInOut_" << vhdlize(lsbInOut);
		setName(name.str()+"_uid"+vhdlize(getNewUId()));

		setCopyrightString("Matei Istoan, 2017");

		//safety checks and warnings
		if(n != m)
			REPORT(INFO, "WARNING: degree of numerator and of denominator are different! "
					<< "This will lead to a less efficient implementation.");
		if(maxDigit < (radix-1))
			REPORT(INFO, "WARNING: used digit set is not maximal!");
		if((radix != 2) && (radix != 4) && (radix != 8))
			THROWERROR("FixEMethodEvaluator: radixes higher than 8 currently not supported!");
		if(maxDigit > radix-1)
			THROWERROR("maximum digit larger than the maximum digit in the redundant digit set!");
		if(intlog2(abs((int)maxDigit)) >= (1<<msbInOut))
			THROWERROR("maximum digit not representable on the given input format!");

		//create a copy of the coefficients of P
		copyVector(_coeffsP, &coeffsP, &mpCoeffsP, maxDegree);
		//create a copy of the coefficients Q
		copyVector(_coeffsQ, &coeffsQ, &mpCoeffsQ, maxDegree);

		//compute the number of iterations needed
		nbIter = msbInOut - lsbInOut + 1;
		//the number of iterations is reduced when using a higher radix
		if(radix > 2)
			nbIter = ceil(1.0*nbIter/log2(radix));
		//add an additional number of iterations to compensate for the errors
		g = intlog2(nbIter);
		nbIter += g;

		//set the format of the internal signals
		//	W^Hat
		GenericSimpleSelectionFunction::getWHatFormat(radix, maxDigit, &msbWHat, &lsbWHat);
		dWHat  = new Signal("dWHat", Signal::wire, true, msbWHat, lsbWHat);
		//	W
		msbW = msbWHat;
		lsbW = lsbInOut - g;
		dW   = new Signal("dW", Signal::wire, true, msbW, lsbW);
		//	D
		msbD = intlog2(radix);
		lsbD = 0;
		dD   = new Signal("dD", Signal::wire, true, msbD, lsbD);
		//	X
		msbX = msbInOut;
		lsbX = lsbInOut;
		dX   = new Signal("dX", Signal::wire, true, msbX, lsbX);
		// DiMultX
		msbDiMX = msbX + (int)ceil(log2(maxDigit));
		lsbDiMX = lsbX;
		dDiMX   = new Signal("dDiMX", Signal::wire, true, msbDiMX, lsbDiMX);

		//add the inputs
		addFixInput("X", true, msbInOut, lsbInOut);
		//add the outputs
		addFixOutput("Y", true, msbInOut, lsbInOut);

		//a helper signal
		vhdl << tab << declare("X_std_lv", msbDiMX-lsbDiMX+1) << " <= std_logic_vector(X);" << endl;

		//create the DiMX signals
		addComment(" ---- create the DiMX signals ----", tab);
		//	the multipliers by constants
		IntConstMult *dimxMult[maxDigit+1];

		//multiply by the positive constants
		addComment(" ---- multiply by the positive constants ----", tab);
		for(size_t i=0; i<=maxDigit; i++)
		{
			//create the multiplication between X and the constant given by i
			dimxMult[i] = new IntConstMult(
											target,					//target
											msbInOut-lsbInOut+1, 	//size of X
											mpz_class(i)			//the constant
											);
			addSubComponent(dimxMult[i]);
			inPortMap  (dimxMult[i], "X", "X_std_lv");
			outPortMap (dimxMult[i], "R", join("X_Mult_", i, "_std_lv"));
			vhdl << tab << instance(dimxMult[i], join("ConstMult_", i));

			vhdl << tab << declareFixPoint(join("X_Mult_", i, "_int"), true, msbDiMX, lsbDiMX)
					<< " <= signed(X_Mult_" << i << "_std_lv);" << endl;
			resizeFixPoint(join("X_Mult_", i), join("X_Mult_", i, "_int"), msbDiMX, lsbDiMX, 1);
		}
		//multiply by the negative constants
		addComment(" ---- multiply by the negative constants ----", tab);
		for(size_t i=1; i<=maxDigit; i++)
		{
			vhdl << tab << declareFixPoint(join("X_Mult_", vhdlize(-((int)i)), "_int"), true, msbDiMX, lsbDiMX)
					<< " <= X_Mult_" << i << " xor (others => '1');" << endl;
			vhdl << tab << declareFixPoint(join("X_Mult_", vhdlize(-((int)i))), true, msbDiMX, lsbDiMX)
					<< " <= X_Mult_" << vhdlize(-((int)i)) << "_int + (" << zg(msbDiMX-lsbDiMX) << "&\"1\");" << endl;
		}

		//iteration 0
		//	initialize the elements of the residual vector
		addComment(" ---- iteration 0 ----", tab);
		for(size_t i=0; i<maxDegree; i++)
		{
			vhdl << tab << declareFixPoint(join("W_0_", i), true, msbW, lsbW) << " <= "
					<< signedFixPointNumber(*mpCoeffsP[i], msbW, lsbW, 0) << ";" << endl;
			vhdl << tab << declareFixPoint(join("D_0_", i), true, msbD, lsbD) << " <= "
					<< zg(msbD-lsbD+1, 0) << ";" << endl;
		}

		//create the computation units
		GenericComputationUnit *cu0, *cuI[maxDegree-2], *cuN;

		//compute unit 0
		cu0 = new GenericComputationUnit(
										target,			//target
										radix, 			//radix
										0,				//index
										-1,				//special case
										maxDigit, 		//maximum digit
										dW,				//signal W
										dX,				//signal X
										dD, 			//signal Di
										coeffsQ[0]		//constant q_i
										);
		for(size_t i=1; i<=(maxDegree-2); i++)
		{
			//compute unit i
			cuI[i-1] = new GenericComputationUnit(
												target,			//target
												radix, 			//radix
												i,				//index
												0,				//special case
												maxDigit, 		//maximum digit
												dW,				//signal W
												dX,				//signal X
												dD, 			//signal Di
												coeffsQ[i]		//constant q_i
												);
		}
		//compute unit n
		cu0 = new GenericComputationUnit(
										target,			//target
										radix, 			//radix
										maxDegree-1,	//index
										+1,				//special case
										maxDigit, 		//maximum digit
										dW,				//signal W
										dX,				//signal X
										dD, 			//signal Di
										coeffsQ[nbIter]	//constant q_i
										);

		//create the selection units
		GenericSimpleSelectionFunction *sel;

		sel = new GenericSimpleSelectionFunction(
												target,			//target
												radix,	 		//radix
												maxDigit, 		//maximum digit
												dW 				//signal W
												);

		//iterations 1 to nbIter
		for(size_t iter=1; iter<=nbIter; iter++)
		{
			addComment(join(" ---- iteration ", iter, " ----"), tab);

			//create computation unit index 0
			//	a special case

			//create computation units index 1 to maxDegree-1
			for(size_t i=0; i<maxDegree; i++)
			{

			}

			//create computation unit index maxDegree
			//	a special case

		}
	}


	FixEMethodEvaluator::~FixEMethodEvaluator()
	{
		for(size_t i=0; i<n; i++)
		{
			mpfr_clear(*mpCoeffsP[i]);
		}
		for(size_t i=0; i<m; i++)
		{
			mpfr_clear(*mpCoeffsQ[i]);
		}
	}


	void FixEMethodEvaluator::copyVector(vector<string> originalVector, vector<string> *newVectorS,
			vector<mpfr_t*> *newVectorMP, size_t maxIndex)
	{
		size_t iterLimit = originalVector.size();

		//safety checks
		if(originalVector.size() > maxIndex)
		{
			THROWERROR("copyVector: maxIndex smaller than the size of the original vector");
		}

		for(size_t i=0; i<iterLimit; i++)
		{
			//create a copy as string
			(*newVectorS).push_back(originalVector[i]);

			//create a copy as MPFR
			mpfr_t tmpMpfr;

			mpfr_init2(tmpMpfr, LARGEPREC);
			//	parse the constant using Sollya
			sollya_obj_t node;
			node = sollya_lib_parse_string(originalVector[i].c_str());
			/* If  parse error throw an exception */
			if (sollya_lib_obj_is_error(node))
			{
				THROWERROR("emulate: Unable to parse string "<< originalVector[i] << " as a numeric constant");
			}
			sollya_lib_get_constant(tmpMpfr, node);
			free(node);

			(*newVectorMP).push_back(&tmpMpfr);
		}
		//fill with zeros, if necessary
		for(size_t i=iterLimit; i<maxIndex; i++)
		{
			//create a copy as string
			(*newVectorS).push_back("0");

			//create a copy as MPFR
			mpfr_t tmpMpfr;

			mpfr_init2(tmpMpfr, LARGEPREC);
			mpfr_set_zero(tmpMpfr, 0);

			(*newVectorMP).push_back(&tmpMpfr);
		}
	}


	void FixEMethodEvaluator::emulate(TestCase * tc)
	{

	}

	OperatorPtr FixEMethodEvaluator::parseArguments(Target *target, std::vector<std::string> &args) {
		return nullptr;
	}

	void FixEMethodEvaluator::registerFactory(){
		UserInterface::add("FixEMethodEvaluator", // name
				"Selection function for the E-method.", //description
				"FunctionApproximation", // category
				"",
				"radix(int): the radix of the digit set being used;\
							maxDigit(int): the maximum digit in the redundant digit set;\
							msbIn(int): MSB of the input;\
							lsbIn(int): LSB of the input"
				"",
				"",
				FixEMethodEvaluator::parseArguments,
				FixEMethodEvaluator::unitTest
		) ;

	}

	TestList FixEMethodEvaluator::unitTest(int index)
	{
		// the static list of mandatory tests
		TestList testStateList;
		vector<pair<string,string>> paramList;



		return testStateList;
	}

} /* namespace flopoco */
