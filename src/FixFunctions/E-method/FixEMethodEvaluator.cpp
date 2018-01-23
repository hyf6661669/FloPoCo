/*
 * FixEMethodEvaluator.cpp
 *
 *  Created on: 7 Dec 2017
 *      Author: mistoan
 */

#include "FixEMethodEvaluator.hpp"

namespace flopoco {

	FixEMethodEvaluator::FixEMethodEvaluator(Target* target, size_t _radix, size_t _maxDigit, int _msbInOut, int _lsbInOut,
		vector<string> _coeffsP, vector<string> _coeffsQ,
		double _delta, bool _scaleInput, double _inputScaleFactor,
		map<string, double> inputDelays)
	: Operator(target), radix(_radix), maxDigit(_maxDigit),
	  	  n(_coeffsP.size()), m(_coeffsQ.size()),
	  	  msbInOut(_msbInOut), lsbInOut(_lsbInOut),
		  coeffsP(_coeffsP), coeffsQ(_coeffsQ),
		  delta(_delta), scaleInput(_scaleInput), inputScaleFactor(_inputScaleFactor),
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
			THROWERROR("radixes higher than 8 currently not supported!");
		if(maxDigit < (radix-1))
			REPORT(INFO, "WARNING: used digit set is not maximal!");
		if(maxDigit > radix-1)
			THROWERROR("maximum digit larger than the maximum digit in the redundant digit set!");
		if((int)maxDigit > (1<<msbInOut))
			THROWERROR("maximum digit not representable on the given input format!");
		if((delta<0 || delta>1))
			THROWERROR("delta must be in the interval [0, 1)!");

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

		//compute the parameters of the algorithm
		REPORT(DEBUG, "compute the parameters of the algorithm");
		setAlgorithmParameters();

		//check P's coefficients
		checkPCoeffs();
		//check Q's coefficients
		checkQCoeffs();
		//check the ranges on the input
		checkX();

		//compute the number of iterations needed
		nbIter = msbInOut - lsbInOut + 1;
		//the number of iterations is reduced when using a higher radix
		if(radix > 2)
			nbIter = ceil(1.0*nbIter/log2(radix));
		//add an additional number of iterations to compensate for the errors
		g = intlog2(nbIter);
		//nbIter += g;

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

		//scale the input X by the factor delta, if necessary
		if(scaleInput)
		{
			int xScaleSize;

			//a helper signal
			vhdl << tab << declare("X_std_lv", msbX-lsbX+1) << " <= std_logic_vector(X);" << endl;

			//scale the input
			//	multiply by inputScaleFactor
			//		set by default to 1/2*alpha, if not otherwise specified by the user
			FixRealKCM *scaleMult = new FixRealKCM(
												 	 target,							//target
													 true,								//signed input
													 msbX,								//msbIn
													 lsbX,								//lsbIn
													 lsbX-g,							//lsbOut
													 std::to_string(inputScaleFactor),	//the constant
													 1.0								//target ulp error
												 	 );
			addSubComponent(scaleMult);
			inPortMap (scaleMult, "X", "X_std_lv");
			outPortMap(scaleMult, "R", "X_scaled_int");
			vhdl << tab << instance(scaleMult, "ScaleConstMult");

			//extend the result of the multiplication to the original signal, if necessary
			//	inputScaleFactor<1, so the scaled signal's msb <= original signal's msb
			xScaleSize = getSignalByName("X_scaled_int")->width();
			if(xScaleSize < (msbX-lsbX+1+g))
			{
				//extension required
				ostringstream xScaleExtension;

				for(int i=xScaleSize; i<(msbX-lsbX+1+g); i++)
					xScaleExtension << "X_scaled_int(" << xScaleSize-1 << ") & ";
				vhdl << tab << declare("X_scaled", msbX-lsbX+1+g) << " <= "
						<< xScaleExtension.str() << "X_scaled_int;" << endl;
			}
			else
			{
				vhdl << tab << declare("X_scaled", msbX-lsbX+g+1) << " <= X_scaled_int;" << endl;
			}

			//update the size of the signals based on X
			//	X
			msbX = msbInOut;
			lsbX = lsbInOut-g;
			dX   = new Signal("dX", Signal::wire, true, msbX, lsbX);
			// DiMultX
			msbDiMX = msbX + (int)ceil(log2(maxDigit));
			lsbDiMX = lsbX;
			dDiMX   = new Signal("dDiMX", Signal::wire, true, msbDiMX, lsbDiMX);
		}
		else
		{
			//no scaling required
			//	just copy the input to an intermediate signal
			vhdl << tab << declare("X_scaled", msbX-lsbX+1) << " <= std_logic_vector(X);" << endl;
		}

		//a helper signal
		vhdl << tab << declare("X_scaled_std_lv", msbX-lsbX+1) << " <= X_scaled;" << endl;
		//a helper signal
		vhdl << tab << declareFixPoint("X_scaled_signed", true, msbX, lsbX) << " <= signed(X_scaled);" << endl;

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
												msbX-lsbX+1,		 	//size of X
												mpz_class(i)			//the constant
												);
				addSubComponent(dimxMult[i]);
				inPortMap  (dimxMult[i], "X", "X_scaled_std_lv");
				outPortMap (dimxMult[i], "R", join("X_Mult_", i, "_std_lv"));
				vhdl << tab << instance(dimxMult[i], join("ConstMult_", i));
			}

			vhdl << tab << declareFixPoint(join("X_Mult_", i, "_int"), true,
					getSignalByName(join("X_Mult_", i, "_std_lv"))->width()+lsbX-1, lsbX)
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
					<< " <= -X_Mult_" << i << ";" << endl;
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

		//iteration 1
		//	the elements of the residual vector are the ones at the previous iteration
		//	shifted by log2(radix) positions to the left
		if(nbIter >= 1)
		{
			addComment(" ---- iteration 1 ----", tab);
			REPORT(DEBUG, "iteration 1");
			for(size_t i=0; i<maxDegree; i++)
			{
				//TODO: this can be more precise than the original algorithm
				//		instead of shifting the previous signal, just generate
				//		the signal from the constant multiplied by two

				//create the residual vector
				/*
				vhdl << tab << declareFixPoint(join("W_1_", i), true, msbW, lsbW) << " <= "
						<< "W_0_" << i << range(msbW-lsbW-ceil(log2(radix)), 0)
						<< " & " << zg(ceil(log2(radix))) << ";" << endl;
				*/
				mpfr_t mpTmp;

				mpfr_init2(mpTmp, LARGEPREC);
				mpfr_mul_2ui(mpTmp, mpCoeffsP[i], 1, GMP_RNDN);

				vhdl << tab << declareFixPoint(join("W_1_", i), true, msbW, lsbW) << " <= "
						<< signedFixPointNumber(mpTmp, msbW, lsbW, 0) << ";" << endl;

				//create the selection unit
				/*
				//inputs
				inPortMap(sel,  "W", join("W_1_", i));
				//outputs
				outPortMap(sel, "D", join("D_1_", i));
				//the instance
				vhdl << tab << instance(sel, join("SEL_1_", i));
				*/
				vhdl << tab << declareFixPoint(join("D_1_", i), true, msbD, lsbD) << " <= "
						<< signedFixPointNumber(mpTmp, msbD, lsbD, 0) << ";" << endl;

				mpfr_clear(mpTmp);
			}
		}

		//iterations 1 to nbIter
		REPORT(DEBUG, "iterations 1 to nbIter");
		for(size_t iter=2; iter<=nbIter; iter++)
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
			inPortMap(cu0, "X",    "X_scaled_signed");
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
				inPortMap(cuI[i-1], "X",    "X_scaled_signed");
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
			inPortMap(cuN, "X",    "X_scaled_signed");
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
			/*
			REPORT(DEBUG, "adding D_" << i << "_" << 0 << " at weight " << (nbIter-1-i)*ceil(log2(radix)));
			bitheap->addSignedBitVector(
										(nbIter-1-i)*ceil(log2(radix)),			//weight
										join("D_", i, "_", 0),					//input signal name
										msbD-lsbD+1								//size
										);
			*/
			REPORT(DEBUG, "adding D_" << i << "_" << 0 << " at weight " << (nbIter+g-1-i)*ceil(log2(radix)));
			bitheap->addSignedBitVector(
										(nbIter+g-1-i)*ceil(log2(radix)),		//weight
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
		REPORT(DEBUG, "retrieve the bits from the bit heap");
		vhdl << tab << declareFixPoint("sum", true, msbW, lsbW) << " <= signed(" <<
				bitheap->getSumName() << range(msbW-lsbW, 0) << ");" << endl;

		//write the result to the output
		REPORT(DEBUG, "write the result to the output, only the bits that we want");
		vhdl << tab << "Y <= sum" << range(msbInOut-lsbInOut+g, g) << ";" << endl;

		REPORT(DEBUG, "constructor completed");
	}


	FixEMethodEvaluator::~FixEMethodEvaluator()
	{
		for(size_t i=0; i<maxDegree; i++)
		{
			mpfr_clear(mpCoeffsP[i]);
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


	void FixEMethodEvaluator::setAlgorithmParameters()
	{
		//set xi
		xi    = 0.5  * (1+delta);
		//set alpha
		alpha = 0.25 * (1-delta);
		//set the scale factor
		//	compute the value of the scale factor in two cases:
		//	1) when scaling is required, and it hasn't been set by the user
		//	2) when scaling isn't required, for the analysis of the parameters
		if(((scaleInput == true) && (inputScaleFactor == -1)) || (scaleInput == false))
		{
			//if the scale factor needs to be automatically be set,
			//	then set if to the maximum admissible limit
			//	we need to first determine the maximum value of the coefficients of Q,
			//	and then subtract it from alpha
			mpfr_t mpTmp, mpAlpha;

			mpfr_inits2(LARGEPREC, mpTmp, mpAlpha, (mpfr_ptr)nullptr);

			mpfr_set_zero(mpTmp, 0);
			for(size_t i=1; i<m; i++)
			{
				if(mpfr_cmp(mpCoeffsQ[i], mpTmp) > 0)
					mpfr_set(mpTmp, mpCoeffsQ[i], GMP_RNDN);
			}

			mpfr_set_d(mpAlpha, alpha, GMP_RNDN);

			mpfr_sub(mpTmp, mpAlpha, mpTmp, GMP_RNDN);
			inputScaleFactor = mpfr_get_d(mpTmp, GMP_RNDN);

			mpfr_clears(mpTmp, mpAlpha, (mpfr_ptr)nullptr);
		}
	}


	void FixEMethodEvaluator::checkPCoeffs()
	{
		mpfr_t mpXi;

		mpfr_init2(mpXi, LARGEPREC);
		mpfr_set_d(mpXi, xi, GMP_RNDN);
		//check that the coefficients of P are smaller than xi
		for(size_t i=0; i<maxDegree; i++)
		{
			if(mpfr_cmp(mpCoeffsP[i], mpXi) > 0)
				THROWERROR("checkPCoeff: coefficient coeffsP[" << i << "]=" << coeffsP[i]
					<< " does not satisfy the constraints");
		}

		mpfr_clear(mpXi);
	}


	void FixEMethodEvaluator::checkQCoeffs()
	{
		mpfr_t mpAlpha, mpLimit;

		mpfr_inits2(LARGEPREC, mpAlpha, mpLimit, (mpfr_ptr)nullptr);

		//the value of alpha
		mpfr_set_d(mpAlpha, alpha, GMP_RNDN);

		//the maximum value allowed for the coefficients
		mpfr_set_d(mpLimit, inputScaleFactor, GMP_RNDN);
		mpfr_sub(mpLimit, mpAlpha, mpLimit, GMP_RNDN);

		//check that coeffsQ[0]=1
		if(mpfr_cmp_ui(mpCoeffsQ[0], 1) != 0)
			THROWERROR("checkQCoeff: coefficient coeffsQ[0]=" << coeffsQ[0]
				<< " should be 1");
		//check that the coefficients of Q are smaller than alpha/2
		for(size_t i=1; i<maxDegree; i++)
		{
			if(mpfr_cmp(mpCoeffsQ[i], mpLimit) > 0)
				THROWERROR("checkQCoeff: coefficient coeffsQ[" << i << "]=" << coeffsQ[i]
					<< " does not satisfy the constraints");
		}

		mpfr_clears(mpAlpha, mpLimit, (mpfr_ptr)nullptr);
	}


	void FixEMethodEvaluator::checkX()
	{
		mpfr_t mpAlpha, mpLimit, mpX, mpTmp;

		mpfr_inits2(LARGEPREC, mpAlpha, mpLimit, mpX, mpTmp, (mpfr_ptr)nullptr);

		//the value of alpha
		mpfr_set_d(mpAlpha, alpha, GMP_RNDN);

		//compute the maximum value of Q's coefficients
		mpfr_set_zero(mpLimit, 0);
		for(size_t i=1; i<m; i++)
		{
			if(mpfr_cmp(mpCoeffsQ[i], mpLimit) > 0)
				mpfr_set(mpLimit, mpCoeffsQ[i], GMP_RNDN);
		}
		mpfr_sub(mpLimit, mpAlpha, mpLimit, GMP_RNDN);

		//check that the largest value that X can take is smaller than alpha/2
		//	largest value X can take
		mpfr_set_ui(mpX, 1, GMP_RNDN);
		mpfr_mul_2si(mpX, mpX, msbInOut, GMP_RNDN);
		//	need to subtract 1 ulp
		mpfr_set_ui(mpTmp, 1, GMP_RNDN);
		mpfr_mul_2si(mpTmp, mpTmp, lsbInOut, GMP_RNDN);
		//	get the actual largest value X can take
		mpfr_sub(mpX, mpX, mpTmp, GMP_RNDN);
		//scale the input
		mpfr_mul_d(mpX, mpX, inputScaleFactor, GMP_RNDN);

		//now perform the test
		if(mpfr_cmp(mpX, mpLimit) > 0)
			THROWERROR("checkX: input format for X, with msb=" << msbInOut
					<< " and lsb=" << lsbInOut << " does not satisfy the constraints");

		mpfr_clears(mpAlpha, mpLimit, mpX, mpTmp, (mpfr_ptr)nullptr);
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
		double dP, dQ, dY, dTmp;

		//initialize the variables
		mpfr_inits2(LARGEPREC, mpX, mpP, mpQ, mpTmp, mpY, (mpfr_ptr)nullptr);

		//initialize P and Q
		mpfr_set_zero(mpP, 0);
		mpfr_set_zero(mpQ, 0);

		//initialize X
		mpfr_set_z(mpX, svX.get_mpz_t(), GMP_RNDN);
		dTmp = mpfr_get_d(mpX, GMP_RNDN);
		//	scale X appropriately, by the amount given by lsbInOut
		mpfr_mul_2si(mpX, mpX, lsbInOut, GMP_RNDN);
		dTmp = mpfr_get_d(mpX, GMP_RNDN);

		//if required, scale the input
		if(scaleInput == true)
			mpfr_mul_d(mpX, mpX, inputScaleFactor, GMP_RNDN);

		//compute P
		for(int i=0; i<(int)n; i++)
		{
			//compute X^i
			mpfr_pow_si(mpTmp, mpX, i, GMP_RNDN);
			dTmp = mpfr_get_d(mpTmp, GMP_RNDN);
			//multiply by coeffsP[i]
			mpfr_mul(mpTmp, mpTmp, mpCoeffsP[i], GMP_RNDN);
			dTmp = mpfr_get_d(mpTmp, GMP_RNDN);

			//add the new term to the sum
			dP = mpfr_get_d(mpP, GMP_RNDN);
			mpfr_add(mpP, mpP, mpTmp, GMP_RNDN);
			dP = mpfr_get_d(mpP, GMP_RNDN);
 		}

		//compute Q
		for(int i=0; i<(int)m; i++)
		{
			//compute X^i
			mpfr_pow_si(mpTmp, mpX, i, GMP_RNDN);
			dTmp = mpfr_get_d(mpTmp, GMP_RNDN);
			//multiply by coeffsQ[i]
			mpfr_mul(mpTmp, mpTmp, mpCoeffsQ[i], GMP_RNDN);
			dTmp = mpfr_get_d(mpTmp, GMP_RNDN);

			//add the new term to the sum
			dQ = mpfr_get_d(mpQ, GMP_RNDN);
			mpfr_add(mpQ, mpQ, mpTmp, GMP_RNDN);
			dQ = mpfr_get_d(mpQ, GMP_RNDN);
		}

		//compute Y = P/Q
		mpfr_div(mpY, mpP, mpQ, GMP_RNDN);
		dY = mpfr_get_d(mpY, GMP_RNDN);

		//scale the result back to an integer
		mpfr_mul_2si(mpY, mpY, -lsbInOut+msbInOut, GMP_RNDN);

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
		double delta;
		bool scaleInput;
		double inputScaleFactor;
		string in, in2;

		UserInterface::parseStrictlyPositiveInt(args, "radix", &radix);
		UserInterface::parseStrictlyPositiveInt(args, "maxDigit", &maxDigit);
		UserInterface::parseInt(args, "msbIn", &msbIn);
		UserInterface::parseInt(args, "lsbIn", &lsbIn);
		UserInterface::parseString(args, "coeffsP", &in);
		UserInterface::parseString(args, "coeffsQ", &in2);
		UserInterface::parseFloat(args, "delta", &delta);
		UserInterface::parseBoolean(args, "scaleInput", &scaleInput);
		UserInterface::parseFloat(args, "inputScaleFactor", &inputScaleFactor);

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

		OperatorPtr result = new FixEMethodEvaluator(target, radix, maxDigit, msbIn, lsbIn,
				coeffsP, coeffsQ, delta, scaleInput, inputScaleFactor);

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
				 coeffsQ(string): colon-separated list of real coefficients of polynomial Q, using Sollya syntax. Example: coeff=\"1.234567890123:sin(3*pi/8)\";\
				 delta(real)=0.5: the value for the delta parameter in the E-method algorithm;\
				 scaleInput(bool)=false: flag showing if the input is to be scaled by the factor delta;\
				 inputScaleFactor(real)=-1: the factor by which the input is scaled"
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
