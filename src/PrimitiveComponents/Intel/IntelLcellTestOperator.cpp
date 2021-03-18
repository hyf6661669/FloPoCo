//
// Created by Martin Kumm on 05.03.21.
//
#include "IntelLcellTestOperator.hpp"
#include "Intel_LCELL.hpp"

using namespace std;

namespace flopoco
{
	IntelLcellTestOperator::IntelLcellTestOperator(Operator *parentOp, Target *target, int param1) : Primitive(parentOp, target)
	{
		srcFileName="IntelLcellTestOperator";

		// definition of the name of the operator
		ostringstream name;
		name << "IntelLcellTestOperator_" << param1 << "_" << param1;
		setName(name.str()); // See also setNameWithFrequencyAndUID()

		addInput("A");
		addInput("B");
		addOutput("Y");

//		vhdl << "Y <= A & B;" << endl;

		string lut_mask="x\"abcdabcdacbdabcd\"";
		Intel_LCELL *lcell = new Intel_LCELL(this,target,lut_mask,false,false);
		inPortMap("dataa","A");
		inPortMap("datab","A");
		inPortMap("datac","A");
		inPortMap("datad","A");
		inPortMap("datae","A");
		inPortMap("dataf","A");
		inPortMap("datag","A");
		inPortMap("sharein","A");
		inPortMap("cin","A");
		outPortMap("combout","Y");
		outPortMap("sumout","Y");
		outPortMap("cout","Y");
		outPortMap("shareout","Y");

		vhdl << lcell->primitiveInstance("lcell1") << endl;

	}

	OperatorPtr IntelLcellTestOperator::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args)
	{
		int param1;

		UserInterface::parseInt(args, "param1", &param1);

		return new IntelLcellTestOperator(parentOp,target,param1);
	}

	void IntelLcellTestOperator::registerFactory(){
		UserInterface::add("IntelLcellTestOperator", // name
						   "Test Operator for Intel ALM Primitives",
						   "Primitives", // categories
						   "",
						   "param1(int): some int parameter",
						   "",
						   IntelLcellTestOperator::parseArguments
		) ;
	}
}