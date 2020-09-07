#ifndef POSITADD2_HPP
#define POSITADD2_HPP
#include "../../Operator.hpp"
#include "../../utils.hpp"

namespace flopoco {
	class PositAdd2 : public Operator {
		public:
			PositAdd2(Target* target, Operator* parentOp, int width, int wES);

			// destructor
			// ~PositAdd2();

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