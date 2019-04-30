/*

  A class that manages fixed-point polynomial approximation of multiple functions for FloPoCo

	At the moment, only on [0,1], no possibility to have it on [-1,1].
	Rationale: no point really as soon as we are not in the degenerate case alpha=0.

  Authors: Matei Istoan

  This file is part of the FloPoCo project

  Initial software.
  Copyright © INSA-Lyon, ENS-Lyon, INRIA, CNRS, UCBL, IMDEA Software

  All rights reserved.

*/


/*
	 The function is assumed to have inputs in [0,1]

	 Stylistic remark: use index i for the subintervals, and j for the degree

*/

#include "MultiFunctApprox.hpp"

namespace flopoco {

	MultiFunctApprox::MultiFunctApprox(vector<FixFunction*> functs_, double targetAccuracy_, int degree_) :
			degree(degree_), targetAccuracy(targetAccuracy_)
	{
		srcFileName="MultiFunctApprox";
		uniqueName_=join(srcFileName, "_", Operator::getNewUId());

		functs.insert(functs.begin(), functs_.begin(), functs_.end());
		needToFreeF = false;

		build();
	}


	MultiFunctApprox::MultiFunctApprox(vector<string> sollyaStringList_, double targetAccuracy_, int degree_) :
			degree(degree_), targetAccuracy(targetAccuracy_)
	{
		srcFileName="MultiFunctApprox";
		uniqueName_=join(srcFileName, "_", Operator::getNewUId());

		for(size_t i=0; i<sollyaStringList_.size(); i++)
		{
			FixFunction *f;
			f = new FixFunction(sollyaStringList_[i], false /* on [0,1]*/);
			functs.push_back(f);
		}
		needToFreeF = true;

		build();
	}


	MultiFunctApprox::~MultiFunctApprox()
	{
		if(needToFreeF)
		{
			for(size_t i=0; i<functs.size(); i++)
				free(functs[i]);
		}

		for(size_t i=0; i<functApprox.size(); i++)
			free(functApprox[i]);
	}


	void MultiFunctApprox::build()
	{
		// initialize the internal variables
		alpha = INT_MIN;
		LSB = INT_MAX;
		for(int i=0; i<degree; i++)
			MSB.push_back(INT_MIN);
		approxErrorBound = -1;

		// create the initial approximations
		//	at the same time, the maximum alpha, the minimum LSB and the maximum MSBs for each degree
		//	no need to determine the maximum degree, as we're trying to force the same one, not getting
		//	an approximation of that degree is an exception, handled later on
		REPORT(INFO,"Generating the initial approximations");
		for(size_t i=0; i<functs.size(); i++)
		{
			REPORT(INFO,"Current function: " << functs[i]->getDescription());

			PiecewisePolyApprox *approx = new PiecewisePolyApprox(functs[i], targetAccuracy, degree);
			functApprox.push_back(approx);

			//bookkeeping
			if(approx->alpha > alpha)
				alpha = approx->alpha;
			if(approx->LSB < LSB)
				LSB = approx->LSB;
			if(approx->MSB.size() > MSB.size())
			{
				for(size_t i=MSB.size(); i<approx->MSB.size(); i++)
					MSB.push_back(INT_MIN);
			}
			for(size_t i=0; i<approx->MSB.size(); i++)
				if(approx->MSB[i] > MSB[i])
					MSB[i] = approx->MSB[i];
			if(approx->approxErrorBound > approxErrorBound)
				approxErrorBound = approx->approxErrorBound;
		}
		REPORT(INFO,"Initial function generation complete");

		// check if all approximations have the same LSB
		//	in case they do not, recompute the approximations, and force the minimum global LSB on all of them
		REPORT(INFO,"Generating the updated approximations");
		for(size_t i=0; i<functs.size(); i++)
		{
			PiecewisePolyApprox *approx = functApprox[i];

			if(approx->LSB > LSB)
				approx->reBuild(LSB);
		}
		REPORT(INFO,"Updated functions generation complete");

		// now do some reporting
		createApproximationsReport();
	}


	void MultiFunctApprox::createApproximationsReport()
	{
		vector<int> totalMemorySizes;
		int totalOutputSize;

		REPORT(INFO,"Parameters and details of the multi-function approximations: ");
		REPORT(INFO,"\t ***********");

		for(size_t i=0; i<functApprox.size(); i++)
		{
			REPORT(INFO,"\t Function: " << functApprox[i]->getFunctionDescription());
			REPORT(INFO,"\t Parameters of the approximation polynomials: ");
			REPORT(INFO,"\t   Degree=" << functApprox[i]->degree	<< "\t alpha=" << functApprox[i]->alpha
					<< "\t maxApproxErrorBound=" << functApprox[i]->approxErrorBound << "\t common coeff LSB="  << functApprox[i]->LSB);

			totalOutputSize=0;
			for(int j=0; j<=degree; j++)
			{
				int size = functApprox[i]->MSB[j]-functApprox[i]->LSB + (functApprox[i]->coeffSigns[j]==0? 1 : 0);
				totalOutputSize += size ;
				REPORT(INFO,"\t\t  MSB["<<j<<"] = " << std::setw(3) << functApprox[i]->MSB[j] << "\t size=" << std::setw(3) << size
						<< (functApprox[i]->coeffSigns[j]==0? "\t variable sign " : "\t constant sign ")
						<< std::setw(2) << functApprox[i]->coeffSigns[j]);
			}
			totalMemorySizes.push_back(totalOutputSize * (1<<functApprox[i]->alpha));

			REPORT(INFO, "\t   Total size of the table for function " << functApprox[i]->getFunctionDescription()
					<< " is " << (1<<functApprox[i]->alpha) << " x " << totalOutputSize << " bits");
			REPORT(INFO,"\t ***********");
		}

		totalOutputSize = 0;
		for(size_t i=0; i<totalMemorySizes.size(); i++)
			totalOutputSize += totalMemorySizes[i];
		REPORT(INFO, "Total size of the tables is " << totalOutputSize << " bits");
	}


	OperatorPtr MultiFunctApprox::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args)
	{
		vector<string> functsList;
		double ta;
		int d;

		UserInterface::parseColonSeparatedStringList(args, "functsList", &functsList);
		UserInterface::parseFloat(args, "targetAcc", &ta);
		UserInterface::parseInt(args, "d", &d);

		MultiFunctApprox *mfa = new MultiFunctApprox(functsList, ta, d);
		cout << "Accuracy is " << mfa->approxErrorBound << " ("<< log2(mfa->approxErrorBound) << " bits)";

		return NULL;
	}


	void MultiFunctApprox::registerFactory()
	{
		UserInterface::add("MultiFunctApprox", // name
				"Helper/Debug feature, does not generate VHDL. A set of uniformly segmented piecewise polynomial approximations of the functions given by functsList,\
 accurate to targetAcc on [0,1)",
				"FunctionApproximation",
				"",
				"functsList(string): set of functions to be evaluated, expressed between double-quotes and separated by colons, for instance \"exp(x):sqrt(x):1/x\";\
				targetAcc(real): the target approximation error of the polynomial WRT the function;\
				d(int): the degree to use",
				"",
				MultiFunctApprox::parseArguments
		);
	}

} /* namespace flopoco */
