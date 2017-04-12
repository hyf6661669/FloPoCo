#include "BaseMultiplierLUT.hpp"

namespace flopoco {


BaseMultiplierLUT::BaseMultiplierLUT(bool isSignedX, bool isSignedY, int wX, int wY){

    srcFileName = "BaseMultiplierLUT";
    uniqueName_ = "BaseMultiplierLUT";

    this->wX = wX;
    this->wY = wY;

    BaseMultiplier::BaseMultiplier(isSignedX,isSignedY);
}

Operator* BaseMultiplierLUT::generateOperator()
{
}

bool BaseMultiplierLUT::shapeValid(int x, int y)
{
}
	
}   //end namespace flopoco
