#include <iostream>
#include <sstream>
#include <string>
#include "gmp.h"
#include "mpfr.h"
#include <vector>
#include <gmpxx.h>
#include <stdio.h>
#include <stdlib.h>
#include "CompressorTypeB.hpp"
#include "../PrimitiveComponents/Xilinx/Xilinx_LUT6.hpp"
#include "../PrimitiveComponents/Xilinx/Xilinx_CARRY4.hpp"
#include "../PrimitiveComponents/Xilinx/Xilinx_LUT_compute.h"

using namespace std;
namespace flopoco{

    CompressorTypeB::CompressorTypeB(Target * target, int width_, int case3_) : Operator(target)
	{
        case3 = case3_;
        setWidth(width_);

        ostringstream name;
        name << "CompressorTypeB" << "_width_" << width;
        setName(name.str());

        addInput("A1",width);
        addInput("A2",width);
        addInput("A3",width);
        addInput("S",2);
        addInput("B0",width);

        addOutput("Y", width+1);


        int needed_cc = ( width / 4 ) + ( width % 4 > 0 ? 1 : 0 ); //no. of required carry chains

        cout << "no of required carry-chains for width=" << width << " is " << needed_cc << endl;

        declare( "cc_s", needed_cc * 4 );
        declare( "cc_di", needed_cc * 4 );
        declare( "cc_co", needed_cc * 4 );
        declare( "cc_o", needed_cc * 4 );



        //init unused carry-chain inputs to zero:
        if(needed_cc*4 > width)
        {
            vhdl << tab << "cc_s(" << needed_cc*4-1 << " downto " << width << ") <= (others => '0');" << endl;
            vhdl << tab << "cc_di(" << needed_cc*4-1 << " downto " << width << ") <= (others => '0');" << endl;
        }    
        vhdl << tab << "cc_di(" << width-1 << " downto 0) <= B0";
        vhdl << endl;

        for(int i=0; i < width; i++)
        {
            declare( join("X",i), 6 );
            vhdl << tab;
            vhdl << join("X",i) << " << ";
            vhdl <<                "A1" << of(i) << " & ";
            vhdl <<                "A2" << of(i) << " & ";
            vhdl <<                "A3" << of(i) << " & ";
            vhdl <<                "S"  << of(0) << " & ";
            vhdl <<                "S"  << of(1) << " & ";
            vhdl <<                "B0" << of(i)  << ";" << std::endl;
        }
        vhdl << endl;

        //create LUTs, except the last LUT:
        //LUT content of the LUTs exept the last LUT:

        //lut_in(0) A1
        //lut_in(1) A2
        //lut_in(2) A3
        //lut_in(3) S0
        //lut_in(4) S1
        //lut_in(5) B0

        lut_op lutop_o6;
        switch(case3)
        {
            //---------------------------------------------------------------------------------------------------------\ /
            case 3: lutop_o6 = ((lut_in(0) & ~lut_in(4) & ~lut_in(3)) | (lut_in(1) & ~lut_in(4) & lut_in(3)) | (~lut_in(0) & lut_in(4) & lut_in(3))) ^ lut_in(5); break; // ((A1 if S=00) or (A2 if S=01) or (A3 if S=10) or (-A1 if S=11)) + B1
            case 4: lutop_o6 = ((lut_in(0) & ~lut_in(4) & ~lut_in(3)) | (lut_in(1) & ~lut_in(4) & lut_in(3)) | (~lut_in(1) & lut_in(4) & lut_in(3))) ^ lut_in(5); break; // ((A1 if S=00) or (A2 if S=01) or (A3 if S=10) or (-A2 if S=11)) + B1
            case 5: lutop_o6 = ((lut_in(0) & ~lut_in(4) & ~lut_in(3)) | (lut_in(1) & ~lut_in(4) & lut_in(3)) | (~lut_in(2) & lut_in(4) & lut_in(3))) ^ lut_in(5); break; // ((A1 if S=00) or (A2 if S=01) or (A3 if S=10) or (-A3 if S=11)) + B1
            case 6: lutop_o6 = ((lut_in(0) & ~lut_in(4) & ~lut_in(3)) | (lut_in(1) & ~lut_in(4) & lut_in(3))                                       ) ^ lut_in(5); break; // ((A1 if S=00) or (A2 if S=01) or (A3 if S=10) or (  0 if S=11)) + B1
        default: std::cout << "ERROR: no valid case selected!" << std::endl << "valid options are:" << std::endl <<  "3 => -A1" << std::endl << "4 => -A2" <<std::endl << "5 => -A3" << endl << "6 =>   0" << std::endl; exit(-1);
        }
        lut_init lutop(lutop_o6);
        for(int i=0; i < width; i++)
        {

            Xilinx_LUT6 *cur_lut = new Xilinx_LUT6( target );
            cur_lut->setGeneric( "init", lutop.get_hex() );

            inPortMap(cur_lut,"i0",join("X",i) + of(0));
            inPortMap(cur_lut,"i1",join("X",i) + of(1));
            inPortMap(cur_lut,"i2",join("X",i) + of(2));
            inPortMap(cur_lut,"i3",join("X",i) + of(3));
            inPortMap(cur_lut,"i4",join("X",i) + of(4));
            inPortMap(cur_lut,"i5",join("X",i) + of(5));

            outPortMap(cur_lut,"o","cc_s" + of(i),false);

            vhdl << cur_lut->primitiveInstance( join("lut",i), this ) << endl;
            //        addToGlobalOpList(cur_lut);
        }

        for( int i = 0; i < needed_cc; i++ )
        {
            Xilinx_CARRY4 *cur_cc = new Xilinx_CARRY4( target );

            inPortMapCst( cur_cc, "cyinit", "'0'" );
            if( i == 0 )
            {
                switch(case3)
                {
                    case 6: vhdl << tab << declare("subtract") << " <= '0'; -- there is no subtraction case. S='11' is addition with 0." << std::endl << std::endl; break;
                    default: vhdl << tab << declare("subtract") << " <= S(0) and S(1); -- in case of subtraction a +1 is needed" << std::endl << std::endl;
                }
                inPortMapCst( cur_cc, "ci", "subtract" ); //carry-in can not be used as AX input is blocked!!
            }
            else
            {
                inPortMap( cur_cc, "ci", "cc_co" + of( i * 4 - 1 ) );
            }
            inPortMap( cur_cc, "di", "cc_di" + range( i * 4 + 3, i * 4 ) );
            inPortMap( cur_cc, "s", "cc_s" + range( i * 4 + 3, i * 4 ) );
            outPortMap( cur_cc, "co", "cc_co" + range( i * 4 + 3, i * 4 ), false);
            outPortMap( cur_cc, "o", "cc_o" + range( i * 4 + 3, i * 4 ), false);

            stringstream cc_name;
            cc_name << "cc_" << i;
            vhdl << cur_cc->primitiveInstance( cc_name.str(), this );
        }

        vhdl << endl;

        vhdl << tab << "R0 <= cc_co(" << width-1 << ") & cc_o(" << width-1 << " downto 0);" << endl;

    }
	
    CompressorTypeB::~CompressorTypeB()
	{
	}
	
    void CompressorTypeB::setWidth(int width_)
    {
        width = width_;

        //adjust size of basic compressor to the size of a variable column compressor of this specific size:
        height.resize(width);

        //no of bits at MSB position
        if(useLastColumn)
            height[0] = 2;
        else
            height[0] = 0;

        for(int i=1; i < width-1; i++)
        {
            height[i] = 4;
        }
        height[width-1] = 4;

        wOut=width+1;

        outputs.resize(wOut);

        if(useLastColumn)
            outputs[0] = 1; //there is one output at MSB
        else
            outputs[0] = 0; //there is no output at MSB

        for(int i=1; i < wOut-1; i++)
        {
            outputs[i] = 2;
        }
        outputs[wOut-1] = 1; //there is one output at LSB
    }
}
	


	
