#pragma once

#include "../Primitive.hpp"

namespace flopoco
{

	class IntelLcellTestOperator : Primitive
	{
	public:

		IntelLcellTestOperator(Operator *parentOp, Target *target, int param1);

		/** Factory method */
		static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target , vector<string> &args);
		/** Register the factory */
		static void registerFactory();

	};

} //namespace