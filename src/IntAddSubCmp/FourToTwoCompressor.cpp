#include <iostream>
#include <sstream>
#include <string>
#include "gmp.h"
#include "mpfr.h"
#include <vector>
#include <gmpxx.h>
#include <stdio.h>
#include <stdlib.h>
#include "FourToTwoCompressor.hpp"
#include "../PrimitiveComponents/Xilinx/Xilinx_LUT6.hpp"
#include "../PrimitiveComponents/Xilinx/Xilinx_LUT_compute.h"

using namespace std;
namespace flopoco{

    FourToTwoCompressor::FourToTwoCompressor(Target * target,int width) : VariableColumnCompressor(target)
	{
        setWidth(width);

        ostringstream name;
        name << "Compressor_4_to_2_width_" << width;
        setName(name.str());

        stringstream xs;
        int w=0;
        for(unsigned i=0;i<height.size();i++)
        {
            w += height[i];
            addInput(join("X",i), height[i]);

            if(i!=0)
            {
                xs<<"& X"<<height.size()-i-1<<" ";
            }
            else
            {
                xs<<"X"<<height.size()-1<<" ";
            }
        }

        xs<<";\n";
        vhdl << tab << declare("X", w) << " <=" << xs.str();

        addOutput("R", wOut);


        vhdl << "--FourToTwoCompressor with width " << width << endl;

        lut_op lutop_o6 = ~( lut_in( 0 )^lut_in( 2 ) ) & ~( lut_in( 1 )^lut_in( 3 ) );
        lut_op lutop_o5 = ( ~( lut_in( 3 ) )&lut_in( 1 ) )
                       | ( ~( lut_in( 3 )^lut_in( 1 ) ) & ( ~( lut_in( 2 ) )&lut_in( 0 ) ) );
        lut_init lutop( lutop_o5, lutop_o6 );

        Xilinx_LUT6_2 *cur_lut = new Xilinx_LUT6_2( target );
        cur_lut->setGeneric( "init", lutop.get_hex() );

        inPortMap(cur_lut,"i1",declare("x1"));
//        outPortMap(cur_lut,"o",declare("x2"));

        addToGlobalOpList(cur_lut);
    }
	
    FourToTwoCompressor::~FourToTwoCompressor()
	{
	}
	
    void FourToTwoCompressor::setWidth(int width)
    {
        this->width = width;

        //adjust size of basic compressor to the size of a variable column compressor of this specific size:
        height.resize(width);
        for(int i=0; i < width-1; i++)
        {
            height[i] = 4;
        }
        height[width-1] = 2;

        wOut=width+1;

        outputs.resize(wOut);
        outputs[0] = 1;
        for(int i=1; i < wOut-2; i++)
        {
            outputs[i] = 2;
        }
        outputs[wOut-2] = 1;
        outputs[wOut-1] = 1;
    }
}
	


	
