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
	  	  n(_coeffsP.size()), m(_coeffsQ.size()),
	  	  msbInOut(_msbInOut), lsbInOut(_lsbInOut),
		  coeffsP(_coeffsP), coeffsQ(_coeffsQ),
		  maxDegree(n>m ? n : m)
	{
		ostringstream name;

		srcFileName = "FixEMethodEvaluator";
		name << "FixEMethodEvaluator_n_" << n << "_m_" << m
				<< "_msbInOut_" << vhdlize(msbInOut) << "_lsbInOut_" << vhdlize(lsbInOut);
		setName(name.str()+"_uid"+vhdlize(getNewUId()));

		useNumericStd_Signed();

		setCopyrightString("Matei Istoan, 2017");

		//safety checks and warnings
		if(n != m)
			REPORT(INFO, "WARNING: degree of numerator and of denominator are different! "
					<< "This will lead to a less efficient implementation.");
		if((radix != 2) && (radix != 4) && (radix != 8))
			THROWERROR("FixEMethodEvaluator: radixes higher than 8 currently not supported!");
		if(maxDigit < (radix-1))
			REPORT(INFO, "WARNING: used digit set is not maximal!");
		if(maxDigit > radix-1)
			THROWERROR("maximum digit larger than the maximum digit in the redundant digit set!");
		if((int)maxDigit > (1<<msbInOut))
			THROWERROR("maximum digit not representable on the given input format!");

		//create a copy of the coefficients of P and Q
		copyVectors();

		REPORT(DEBUG, "coefficients of P, as mp numbers:");
		for(size_t i=0; i<n; i++)
		{
			double tmpD = mpfr_get_d(mpCoeffsP[i], GMP_RNDN);
			REPORT(DEBUG, "" << tmpD);
		}
		REPORT(DEBUG, "coefficients of Q, as mp numbers:");
		for(size_t i=0; i<m; i++)
		{
			double tmpD = mpfr_get_d(mpCoeffsQ[i], GMP_RNDN);
			REPORT(DEBUG, "" << tmpD);
		}

		//safety checks and warnings
		if(mpfr_cmp_si(mpCoeffsQ[0], 1) != 0)
			THROWERROR("element 0 of coeffsQ must be 1!");

		//compute the number of iterations needed
		nbIter = msbInOut - lsbInOut + 1;
		//the number of iterations is reduced when using a higher radix
		if(radix > 2)
			nbIter = ceil(1.0*nbIter/log2(radix));
		//add an additional number of iterations to compensate for the errors
		g = intlog2(nbIter);
		nbIter += g;

		//set the format of the internal signals
		REPORT(DEBUG, "set the format of the internal signals");
		//	W^Hat
		GenericSimpleSelectionFunction::getWHatFormat(radix, maxDigit, &msbWHat, &lsbWHat);
		dWHat  = new Signal("dWHat", Signal::wire, true, msbWHat, lsbWHat);
		//	W
		msbW = msbWHat;
		lsbW = lsbInOut - g;
		dW   = new Signal("dW", Signal::wire, true, msbW, lsbW);
		//	D
		msbD = ceil(log2(radix));
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
		addFixOutput("Y", true, msbInOut, lsbInOut, 2);

		//a helper signal
		vhdl << tab << declare("X_std_lv", msbDiMX-lsbDiMX+1) << " <= std_logic_vector(X);" << endl;

		//create the DiMX signals
		REPORT(DEBUG, "create the DiMX signals");
		addComment(" ---- create the DiMX signals ----", tab);
		//	the multipliers by constants
		IntConstMult *dimxMult[maxDigit+1];

		//multiply by the positive constants
		REPORT(DEBUG, "multiply by the positive constants");
		addComment(" ---- multiply by the positive constants ----", tab);
		for(size_t i=0; i<=maxDigit; i++)
		{
			//create the multiplication between X and the constant given by i
			REPORT(DEBUG, "create the multiplication between X and the constant " << i);
			if(i == 0)
			{
				vhdl << tab << declareFixPoint(join("X_Mult_", i, "_std_lv"), true, msbDiMX, lsbDiMX)
						<< " <= " << zg(msbDiMX-lsbDiMX+1) << ";" << endl;
			}
			else
			{
				dimxMult[i] = new IntConstMult(
												target,					//target
												msbInOut-lsbInOut+1, 	//size of X
												mpz_class(i)			//the constant
												);
				addSubComponent(dimxMult[i]);
				inPortMap  (dimxMult[i], "X", "X_std_lv");
				outPortMap (dimxMult[i], "R", join("X_Mult_", i, "_std_lv"));
				vhdl << tab << instance(dimxMult[i], join("ConstMult_", i));
			}

			vhdl << tab << declareFixPoint(join("X_Mult_", i, "_int"), true, msbDiMX, lsbDiMX)
					<< " <= signed(X_Mult_" << i << "_std_lv);" << endl;
			resizeFixPoint(join("X_Mult_", i), join("X_Mult_", i, "_int"), msbDiMX, lsbDiMX, 1);
		}
		//multiply by the negative constants
		REPORT(DEBUG, "multiply by the negative constants");
		addComment(" ---- multiply by the negative constants ----", tab);
		for(int i=1; i<=(int)maxDigit; i++)
		{
			REPORT(DEBUG, "create the multiplication between X and the constant " << -i);
			vhdl << tab << declareFixPoint(join("X_Mult_", vhdlize(-i)), true, msbDiMX, lsbDiMX)
					<< " <= not(X_Mult_" << vhdlize(-i) << ") + (" << zg(msbDiMX-lsbDiMX) << "&\"1\");" << endl;
		}

		//iteration 0
		//	initialize the elements of the residual vector
		addComment(" ---- iteration 0 ----", tab);
		REPORT(DEBUG, "iteration 0");
		for(size_t i=0; i<maxDegree; i++)
		{
			vhdl << tab << declareFixPoint(join("W_0_", i), true, msbW, lsbW) << " <= "
					<< signedFixPointNumber(mpCoeffsP[i], msbW, lsbW, 0) << ";" << endl;
			vhdl << tab << declareFixPoint(join("D_0_", i), true, msbD, lsbD) << " <= "
					<< zg(msbD-lsbD+1, 0) << ";" << endl;
		}

		//create the computation units
		REPORT(DEBUG, "create the computation units");
		GenericComputationUnit *cu0, *cuI[maxDegree-2], *cuN;

		//compute unit 0
		REPORT(DEBUG, "create the computation unit 0");
		cu0 = new GenericComputationUnit(
										target,			//target
										radix, 			//radix
										maxDigit,		//maximum digit
										0,				//index
										-1,				//special case
										dW,				//signal W
										dX,				//signal X
										dD, 			//signal Di
										coeffsQ[0]		//constant q_i
										);
		addSubComponent(cu0);
		for(size_t i=1; i<=(maxDegree-2); i++)
		{
			//compute unit i
			REPORT(DEBUG, "create the computation unit " << i);
			cuI[i-1] = new GenericComputationUnit(
												target,			//target
												radix, 			//radix
												maxDigit, 		//maximum digit
												i,				//index
												0,				//special case
												dW,				//signal W
												dX,				//signal X
												dD, 			//signal Di
												coeffsQ[i]		//constant q_i
												);
			addSubComponent(cuI[i-1]);
		}
		//compute unit n
		REPORT(DEBUG, "create the computation unit n");
		cuN = new GenericComputationUnit(
										target,					//target
										radix, 					//radix
										maxDigit, 				//maximum digit
										maxDegree-1,			//index
										+1,						//special case
										dW,						//signal W
										dX,						//signal X
										dD, 					//signal Di
										coeffsQ[maxDegree-1]	//constant q_i
										);
		addSubComponent(cuN);

		//create the selection units
		REPORT(DEBUG, "create the selection unit");
		GenericSimpleSelectionFunction *sel;

		sel = new GenericSimpleSelectionFunction(
												target,			//target
												radix,	 		//radix
												maxDigit, 		//maximum digit
												dW 				//signal W
												);
		addSubComponent(sel);

		//iterations 1 to nbIter
		REPORT(DEBUG, "iterations 1 to nbIter");
		for(size_t iter=1; iter<=nbIter; iter++)
		{
			REPORT(DEBUG, "iteration " << iter);
			addComment(join(" ---- iteration ", iter, " ----"), tab);

			//create computation unit index 0
			REPORT(DEBUG, "create computation unit index 0");
			//	a special case
			//inputs
			inPortMap(cu0, "Wi",   join("W_", iter-1, "_0"));
			inPortMap(cu0, "D0",   join("D_", iter-1, "_0"));
			inPortMap(cu0, "Di",   join("D_", iter-1, "_0"));
			inPortMap(cu0, "Dip1", join("D_", iter-1, "_1"));
			inPortMap(cu0, "X",    "X");
			for(int i=(-(int)maxDigit); i<=(int)maxDigit; i++)
			{
				inPortMap(cu0, join("X_Mult_", vhdlize(i)), join("X_Mult_", vhdlize(i)));
			}
			//outputs
			outPortMap(cu0, "Wi_next", join("W_", iter, "_0"));
			//the instance
			vhdl << tab << instance(cu0, join("CU_", iter, "_0"));

			//create computation units index 1 to maxDegree-1
			REPORT(DEBUG, "create computation units index 1 to maxDegree-1");
			for(size_t i=1; i<=(maxDegree-2); i++)
			{
				REPORT(DEBUG, "create computation unit index " << i);
				//inputs
				inPortMap(cuI[i-1], "Wi",   join("W_", iter-1, "_", i));
				inPortMap(cuI[i-1], "D0",   join("D_", iter-1, "_0"));
				inPortMap(cuI[i-1], "Di",   join("D_", iter-1, "_", i));
				inPortMap(cuI[i-1], "Dip1", join("D_", iter-1, "_", i+1));
				inPortMap(cuI[i-1], "X",    "X");
				for(int j=(-(int)maxDigit); j<=(int)maxDigit; j++)
				{
					inPortMap(cuI[i-1], join("X_Mult_", vhdlize(j)), join("X_Mult_", vhdlize(j)));
				}
				//outputs
				outPortMap(cuI[i-1], "Wi_next", join("W_", iter, "_", i));
				//the instance
				vhdl << tab << instance(cuI[i-1], join("CU_", iter, "_", i));
			}

			//create computation unit index maxDegree
			REPORT(DEBUG, "create computation unit index maxDegree");
			//	a special case
			//inputs
			inPortMap(cuN, "Wi",   join("W_", iter-1, "_", maxDegree-1));
			inPortMap(cuN, "D0",   join("D_", iter-1, "_0"));
			inPortMap(cuN, "Di",   join("D_", iter-1, "_", maxDegree-1));
			inPortMap(cuN, "X",    "X");
			//outputs
			outPortMap(cuN, "Wi_next", join("W_", iter, "_", maxDegree-1));
			//the instance
			vhdl << tab << instance(cuN, join("CU_", iter-1, "_", maxDegree-1));

			//create the selection units index 0 to maxDegree-1
			REPORT(DEBUG, "create the selection units index 0 to maxDegree-1");
			for(size_t i=0; i<maxDegree; i++)
			{
				//inputs
				inPortMap(sel,  "W", join("W_", iter, "_", i));
				//outputs
				outPortMap(sel, "D", join("D_", iter, "_", i));
				//the instance
				vhdl << tab << instance(sel, join("SEL_", iter, "_", i));
			}
		}

		//compute the final result
		REPORT(DEBUG, "compute the final result");
		BitHeap *bitheap = new BitHeap(
										this,											// parent operator
										msbW-lsbW+1,									// maximum weight
										false, 											// enable supertiles
										join("Bitheap_"+name.str()+"_", getNewUId())	// bitheap name
										);
		//add the digits of the intermediate computations
		REPORT(DEBUG, "add the digits of the intermediate computations");
		for(int i=(nbIter-1); i>0; i--)
		{
			REPORT(DEBUG, "adding D_" << i << "_" << 0 << " at weight " << (nbIter-1-i)*ceil(log2(radix)));
			bitheap->addSignedBitVector(
										(nbIter-1-i)*ceil(log2(radix)),			//weight
										join("D_", i, "_", 0),					//input signal name
										msbD-lsbD+1								//size
										);
		}
		//add the rounding bit
		REPORT(DEBUG, "add the rounding bit");
		bitheap->addConstantOneBit(g-1);
		//compress the bitheap
		REPORT(DEBUG, "compress the bitheap");
		bitheap->generateCompressorVHDL();

		//retrieve the bits we want from the bit heap
		REPORT(DEBUG, "retrieve the bits we want from the bit heap");
		vhdl << tab << declareFixPoint("sum", true, msbW, lsbW) << " <= signed(" <<
				bitheap->getSumName() << range(msbW-lsbW, 0) << ");" << endl;

		//write the result to the output
		REPORT(DEBUG, "write the result to the output");
		vhdl << tab << "Y <= sum" << range(msbInOut-lsbInOut+g, g) << ";" << endl;

		REPORT(DEBUG, "constructor completed");
	}


	FixEMethodEvaluator::~FixEMethodEvaluator()
	{
		for(size_t i=0; i<n; i++)
		{
			mpfr_clear(mpCoeffsP[i]);
		}
		for(size_t i=0; i<m; i++)
		{
			mpfr_clear(mpCoeffsQ[i]);
		}
	}


	void FixEMethodEvaluator::copyVectors()
	{
		size_t iterLimit = coeffsP.size();

		//copy the coefficients of P
		for(size_t i=0; i<iterLimit; i++)
		{
			//create a copy as MPFR
			mpfr_init2(mpCoeffsP[i], LARGEPREC);
			//	parse the constant using Sollya
			sollya_obj_t node;
			node = sollya_lib_parse_string(coeffsP[i].c_str());
			/* If  parse error throw an exception */
			if (sollya_lib_obj_is_error(node))
			{
				THROWERROR("emulate: Unable to parse string "<< coeffsP[i] << " as a numeric constant");
			}
			sollya_lib_get_constant(mpCoeffsP[i], node);
			free(node);
		}
		//fill with zeros, if necessary
		for(size_t i=iterLimit; i<maxDegree; i++)
		{
			//create a copy as string
			coeffsP.push_back(string("0"));

			//create a copy as MPFR
			mpfr_init2(mpCoeffsP[i], LARGEPREC);
			mpfr_set_zero(mpCoeffsP[i], 0);
		}

		iterLimit = coeffsQ.size();
		//copy the coefficients of Q
		for(size_t i=0; i<iterLimit; i++)
		{
			//create a copy as MPFR
			mpfr_init2(mpCoeffsQ[i], LARGEPREC);
			//	parse the constant using Sollya
			sollya_obj_t node;
			node = sollya_lib_parse_string(coeffsQ[i].c_str());
			/* If  parse error throw an exception */
			if (sollya_lib_obj_is_error(node))
			{
				THROWERROR("emulate: Unable to parse string "<< coeffsQ[i] << " as a numeric constant");
			}
			sollya_lib_get_constant(mpCoeffsQ[i], node);
			free(node);
		}
		//fill with zeros, if necessary
		for(size_t i=iterLimit; i<maxDegree; i++)
		{
			//create a copy as string
			coeffsQ.push_back(string("0"));

			//create a copy as MPFR
			mpfr_init2(mpCoeffsQ[i], LARGEPREC);
			mpfr_set_zero(mpCoeffsQ[i], 0);
		}
	}


	void FixEMethodEvaluator::emulate(TestCase * tc)
	{
		//get the inputs from the TestCase
		mpz_class svX   = tc->getInputValue("X");

		//manage signed digits
		mpz_class big1X      = (mpz_class(1) << (msbInOut-lsbInOut+1));
		mpz_class big1Xp     = (mpz_class(1) << (msbInOut-lsbInOut));

		//handle the signed inputs
		if(svX >= big1Xp)
			svX -= big1X;

		// compute the multiple-precision output
		mpz_class svYd, svYu;
		mpfr_t mpX, mpP, mpQ, mpTmp, mpY;

		//initialize the variables
		mpfr_inits2(LARGEPREC, mpX, mpP, mpQ, mpTmp, mpY, (mpfr_ptr)nullptr);

		//initialize P and Q
		mpfr_set_zero(mpP, 0);
		mpfr_set_zero(mpQ, 0);

		//initialize X
		mpfr_set_z(mpX, svX.get_mpz_t(), GMP_RNDN);
		//	scale X appropriately, by the amount given by lsbInOut
		mpfr_mul_2si(mpTmp, mpTmp, lsbW, GMP_RNDN);

		//compute P
		for(size_t i=0; i<n; i++)
		{
			//compute X^i
			mpfr_pow_si(mpTmp, mpX, i, GMP_RNDN);
			//multiply by coeffsP[i]
			mpfr_mul(mpTmp, mpTmp, mpCoeffsP[i], GMP_RNDN);

			//add the new term to the sum
			mpfr_add(mpP, mpP, mpTmp, GMP_RNDN);
 		}

		//compute Q
		for(size_t i=0; i<m; i++)
		{
			//compute X^i
			mpfr_pow_si(mpTmp, mpX, i, GMP_RNDN);
			//multiply by coeffsQ[i]
			mpfr_mul(mpTmp, mpTmp, mpCoeffsQ[i], GMP_RNDN);

			//add the new term to the sum
			mpfr_add(mpQ, mpQ, mpTmp, GMP_RNDN);
		}

		//compute Y = P/Q
		mpfr_div(mpY, mpP, mpQ, GMP_RNDN);

		//scale the result back to an integer
		mpfr_mul_2si(mpY, mpY, -lsbInOut, GMP_RNDN);

		//round the result
		mpfr_get_z(svYd.get_mpz_t(), mpY, GMP_RNDD);
		mpfr_get_z(svYu.get_mpz_t(), mpY, GMP_RNDU);

		//handle the signed outputs
		if(svYd < 0)
			svYd += big1X;
		if(svYu < 0)
			svYu += big1X;

		//only use the required bits
		svYd &= (big1X-1);
		svYu &= (big1X-1);

		//add this expected output to the TestCase
		tc->addExpectedOutput("Y", svYd);
		tc->addExpectedOutput("Y", svYu);

		//cleanup
		mpfr_clears(mpX, mpP, mpQ, mpTmp, mpY, (mpfr_ptr)nullptr);
	}

	OperatorPtr FixEMethodEvaluator::parseArguments(Target *target, std::vector<std::string> &args) {
		int radix;
		int maxDigit;
		int msbIn;
		int lsbIn;
		vector<string> coeffsP;
		vector<string> coeffsQ;
		string in, in2;

		UserInterface::parseStrictlyPositiveInt(args, "radix", &radix);
		UserInterface::parseStrictlyPositiveInt(args, "maxDigit", &maxDigit);
		UserInterface::parseInt(args, "msbIn", &msbIn);
		UserInterface::parseInt(args, "lsbIn", &lsbIn);
		UserInterface::parseString(args, "coeffsP", &in);
		UserInterface::parseString(args, "coeffsQ", &in2);

		stringstream ss(in);
		string substr;
		while(std::getline(ss, substr, ':'))
		{
			coeffsP.insert(coeffsP.begin(), std::string(substr));
			//coeffsP.push_back(std::string(substr));
		}

		stringstream ss2(in2);
		while(std::getline(ss2, substr, ':'))
		{
			coeffsQ.insert(coeffsQ.begin(), std::string(substr));
			//coeffsQ.push_back(std::string(substr));
		}

		OperatorPtr result = new FixEMethodEvaluator(target, radix, maxDigit, msbIn, lsbIn, coeffsP, coeffsQ);

		return result;
	}

	void FixEMethodEvaluator::registerFactory(){
		UserInterface::add("FixEMethodEvaluator", // name
				"A hardware implementation of the E-method for the evaluation of polynomials and rational polynomials.", //description
				"FunctionApproximation", // category
				"",
				"radix(int): the radix of the digit set being used;\
				 maxDigit(int): the maximum digit in the redundant digit set;\
				 msbIn(int): MSB of the input;\
				 lsbIn(int): LSB of the input;\
				 coeffsP(string): colon-separated list of real coefficients of polynomial P, using Sollya syntax. Example: coeff=\"1.234567890123:sin(3*pi/8)\";\
				 coeffsQ(string): colon-separated list of real coefficients of polynomial Q, using Sollya syntax. Example: coeff=\"1.234567890123:sin(3*pi/8)\""
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
