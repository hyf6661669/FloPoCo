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
	}


	void MultiFunctApprox::build()
	{

	}


	OperatorPtr MultiFunctApprox::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args)
	{
		vector<string> functsList;
		double ta;
		int d;

		UserInterface::parseColonSeparatedStringList(args, "functs", &functsList);
		UserInterface::parseFloat(args, "targetAcc", &ta);
		UserInterface::parseInt(args, "d", &d);

		MultiFunctApprox *ppa = new MultiFunctApprox(functsList, ta, d);
		cout << "Accuracy is " << ppa->approxErrorBound << " ("<< log2(ppa->approxErrorBound) << " bits)";

		return NULL;
	}


	void MultiFunctApprox::registerFactory()
	{
		UserInterface::add("MultiFunctApprox", // name
				"Helper/Debug feature, does not generate VHDL. A set of uniformly segmented piecewise polynomial approximations of the functions given by functsList,\
				accurate to targetAcc on [0,1)",
				"FunctionApproximation",
				"",
				"\
				functsList(string): set of functions to be evaluated, expressed between double-quotes and separated by colons, for instance \"exp(x):sqrt(x):1/x\";\
				targetAcc(real): the target approximation error of the polynomial WRT the function;\
				d(int): the degree to use",
				"",
				PiecewisePolyApprox::parseArguments
		);
	}

} /* namespace flopoco */
