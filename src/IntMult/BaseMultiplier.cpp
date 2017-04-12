#include "BaseMultiplier.hpp"

namespace flopoco {


BaseMultiplier::BaseMultiplier(bool isSignedX, bool isSignedY)
{
    srcFileName = "BaseMultiplier";
    uniqueName_ = "BaseMultiplier";

    this->isSignedX = isSignedX;
    this->isSignedY = isSignedY;
}


	
}   //end namespace flopoco
