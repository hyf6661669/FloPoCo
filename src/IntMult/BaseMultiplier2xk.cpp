#include "BaseMultiplier2xk.hpp"
#include "../PrimitiveComponents/Xilinx/Xilinx_LUT6.hpp"
#include "../PrimitiveComponents/Xilinx/Xilinx_CARRY4.hpp"
#include "../PrimitiveComponents/Xilinx/Xilinx_LUT_compute.h"

namespace flopoco {


BaseMultiplier2xk::BaseMultiplier2xk(bool isSignedX, bool isSignedY, int k) : BaseMultiplier(isSignedX,isSignedY)
{

    srcFileName = "BaseMultiplier2xk";
    uniqueName_ = "BaseMultiplier2xk";

    this->wX = 2;
    this->wY = k;

}

Operator* BaseMultiplier2xk::generateOperator(Target* target)
{
    return new BaseMultiplier2xkOp(target, isSignedX, isSignedY, wY);
}
	

BaseMultiplier2xkOp::BaseMultiplier2xkOp(Target* target, bool isSignedX, bool isSignedY, int width) : Operator(target)
{
    ostringstream name;
    name << "BaseMultiplier2x" << width;
    setName(name.str());

    addInput("X", 2, true);
    addInput("Y", width, true);
    addOutput("R", width+1, 1, true);

    int needed_cc = ( width / 4 ) + ( width % 4 > 0 ? 1 : 0 ); //no. of required carry chains

    declare( "cc_s", needed_cc * 4 );
    declare( "cc_di", needed_cc * 4 );

    for(int i=0; i < width; i++)
    {
        //LUT content of the LUTs:
        lut_op lutop_o6 = (lut_in(0) & lut_in(1)) ^ (lut_in(2) & lut_in(3)); //xor of two partial products
        lut_op lutop_o5 = lut_in(0) & lut_in(1); //and of first partial product
        lut_init lutop( lutop_o5, lutop_o6 );

        Xilinx_LUT6_2 *cur_lut = new Xilinx_LUT6_2( target );
        cur_lut->setGeneric( "init", lutop.get_hex() );

        inPortMap(cur_lut,"i0","X" + of(2*i));
        inPortMap(cur_lut,"i1","Y" + of(2*i));
        inPortMap(cur_lut,"i2","X" + of(2*i+1));
        inPortMap(cur_lut,"i3","Y" + of(2*i+1));
        inPortMapCst(cur_lut, "i4","'0'");
        inPortMapCst(cur_lut, "i5","'1'");

        outPortMap(cur_lut,"o5","cc_di" + of(i),false);
        outPortMap(cur_lut,"o6","cc_s" + of(i),false);

        vhdl << cur_lut->primitiveInstance( join("lut",i), this ) << endl;
    }
}

}   //end namespace flopoco

