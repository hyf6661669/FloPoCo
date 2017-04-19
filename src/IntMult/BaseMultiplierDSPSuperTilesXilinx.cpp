#include "BaseMultiplierDSPSuperTilesXilinx.hpp"
#include "../PrimitiveComponents/Xilinx/Xilinx_LUT6.hpp"
#include "../PrimitiveComponents/Xilinx/Xilinx_CARRY4.hpp"
#include "../PrimitiveComponents/Xilinx/Xilinx_LUT_compute.h"

namespace flopoco {


BaseMultiplierDSPSuperTilesXilinx::BaseMultiplierDSPSuperTilesXilinx(bool isSignedX, bool isSignedY, TILE_SHAPE shape) : BaseMultiplier(isSignedX,isSignedY)
{

    srcFileName = "BaseMultiplierDSPSuperTilesXilinx";
    uniqueName_ = "BaseMultiplierDSPSuperTilesXilinx";

    this->flipXY = flipXY;
    this->shape = shape;

	switch(shape)
	{
		case SHAPE_A:
			wX = 31;
			wY = 34;
			break;
		case SHAPE_B:
			wX = 41;
			wY = 41;
			break;
		//...
		default:
			throw string("Error in " + srcFileName + ": shape " + shape + " unknown");		
	}
	
    if(flipXY)
    {
    	int tmp = wX;
        wX = wY;
        wY = tmp;
    }

}

Operator* BaseMultiplierDSPSuperTilesXilinx::generateOperator(Target* target)
{
    return new BaseMultiplierDSPSuperTilesXilinxOp(target, isSignedX, isSignedY, wX, wY, flipXY);
}
	

BaseMultiplierDSPSuperTilesXilinxOp::BaseMultiplierDSPSuperTilesXilinxOp(Target* target, bool isSignedX, bool isSignedY, int wX, int wY, bool flipXY) : Operator(target)
{
    ostringstream name;
    name << "BaseMultiplierDSPSuperTilesXilinx";
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

    if((isSignedX == true) || (isSignedY == true)) throw string("unsigned inputs currently not supported by BaseMultiplierDSPSuperTilesXilinxOp, sorry");

	declare("D1", 41); //temporary output of the first DSP
	declare("D2", 41); //temporary output of the second DSP

	switch(shape)
	{
		case SHAPE_A:
		    wR = 42;		    
			//total operation is: (X(30 downto 7) * Y(16 downto 0) << 17) + (X(23 downto 0) * Y(33 downto 17) << 17)
			//realized as:        (X(30 downto 7) * Y(16 downto 0) + (X(23 downto 0) * Y(33 downto 17)) << 17
			vhdl << tab << "D1 <= std_logic_vector(unsigned(X(23 downto 0)) * unsigned(Y(33 downto 17)));" << endl;	    
			vhdl << tab << "D2 <= std_logic_vector(unsigned(X(30 downto 7)) * unsigned(Y(16 downto 0)));" << endl;
			vhdl << tab << "R <= std_logic_vector(unsigned('0' & D1) + unsigned('0' & D2)) & \"00000000000000000\";" << endl;
			break;
		case SHAPE_B:
			//total operation is: (X(23 downto 0) * Y(40 downto 23) << 17) + (X(23 downto 0) * Y(33 downto 17) << 17)
		    wR = 42;
			break;
		case SHAPE_C:
		    wR = 42;
			break;
		case SHAPE_D:
		    wR = 42;
			break;
		case SHAPE_E:
		    wR = 58;
			break;
		//...
		default:
			throw string("Error in " + srcFileName + ": shape " + shape + " unknown");		
	}
    addOutput("R", wR);
	
    vhdl << tab << "R <= std_logic_vector(unsigned(X) * unsigned(Y));" << endl;
}

}   //end namespace flopoco

