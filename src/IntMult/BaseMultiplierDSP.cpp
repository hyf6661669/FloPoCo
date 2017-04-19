#include "BaseMultiplierDSP.hpp"
#include "../PrimitiveComponents/Xilinx/Xilinx_LUT6.hpp"
#include "../PrimitiveComponents/Xilinx/Xilinx_CARRY4.hpp"
#include "../PrimitiveComponents/Xilinx/Xilinx_LUT_compute.h"

namespace flopoco {


BaseMultiplierDSP::BaseMultiplierDSP(bool isSignedX, bool isSignedY, int wX, int wY, bool flipXY) : BaseMultiplier(isSignedX,isSignedY)
{

    srcFileName = "BaseMultiplierDSP";
    uniqueName_ = "BaseMultiplierDSP" + std::to_string(wX) + "x" + std::to_string(wY);

    this->flipXY = flipXY;

    if(!flipXY)
    {
        this->wX = wX;
        this->wY = wY;
    }
    else
    {
        this->wX = wY;
        this->wY = wX;
    }

}

Operator* BaseMultiplierDSP::generateOperator(Target* target)
{
    return new BaseMultiplierDSPOp(target, isSignedX, isSignedY, wX, wY, flipXY);
}
	

BaseMultiplierDSPOp::BaseMultiplierDSPOp(Target* target, bool isSignedX, bool isSignedY, int wX, int wY, bool flipXY) : Operator(target)
{
    ostringstream name;
    name << "BaseMultiplierDSP";
    setName(name.str());

    string in1,in2;

    if(!flipXY)
    {
        this->wX = wX;
        this->wY = wY;
        in1 = "X";
        in2 = "Y";
    }
    else
    {
        this->wX = wY;
        this->wY = wX;
        in1 = "Y";
        in2 = "X";
    }

    addInput(in1, wX, true);
    addInput(in2, wY, true);

    if((isSignedX == true) || (isSignedY == true)) throw string("unsigned inputs currently not supported by BaseMultiplierDSPOp, sorry");

    wR = wX + wY; //ToDo: addjust this for signed
    addOutput("R", wR);

    vhdl << tab << "R <= std_logic_vector(unsigned(X) * unsigned(Y));" << endl;
}

}   //end namespace flopoco

