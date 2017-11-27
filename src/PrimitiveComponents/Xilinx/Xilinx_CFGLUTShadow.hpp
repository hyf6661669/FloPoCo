#ifndef Xilinx_CFGLUT5_H
#define Xilinx_CFGLUT5_H

#include "PrimitiveComponents/Primitive.hpp"
#include "PrimitiveComponents/Xilinx/Xilinx_Primitive.hpp"

namespace flopoco {
	// new operator class declaration
    class Xilinx_CFGLUTShadow : public Operator {

	public:
		// definition of some function for the operator    

		// constructor, defined there with two parameters (default value 0 for each)
        Xilinx_CFGLUTShadow(Target* target, bool simpleInterface=true);

	};


}//namespace

#endif
