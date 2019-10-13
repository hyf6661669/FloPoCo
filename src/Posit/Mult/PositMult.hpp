#ifndef POSITMULT_HPP
#define POSITMULT_HPP
#include "../../Operator.hpp"
#include "../../utils.hpp"

namespace flopoco {
	class PositMult : public Operator {
		public:
			PositMult(Target* target, Operator* parentOp, int width, int wES);

			// destructor
			// ~PositMult();

			/** Factory method that parses arguments and calls the constructor */
			static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target , vector<string> &args);

			/** Factory register method */ 
			static void registerFactory();

			/** emulate() function to be shared by various implementations */
			void emulate(TestCase * tc);

		private:
			int width;
			int wES;

	};


}//namespace
#endif