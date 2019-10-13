#ifndef POSITDECODER_HPP
#define POSITDECODER_HPP
#include "../../Operator.hpp"
#include "../../utils.hpp"

namespace flopoco {
	class PositDecoder : public Operator {
		public:
	        PositDecoder(Target* target, Operator* parentOp, int width, int wES);
	  
	  		// destructor
    		// ~PositDecoder();

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