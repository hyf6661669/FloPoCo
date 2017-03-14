#include <iostream>
#include <sstream>
#include <string>
#include "gmp.h"
#include "mpfr.h"
#include <vector>
#include <gmpxx.h>
#include <stdio.h>
#include <stdlib.h>
#include "VariableColumnCompressor.hpp"


using namespace std;
namespace flopoco{

    VariableColumnCompressor::VariableColumnCompressor(Target * target) : BasicCompressor(target)
	{
//        vhdl << "--VariableColumnCompressor" << endl;
	}
	
	VariableColumnCompressor::~VariableColumnCompressor()
	{
	}
	
    std::ostream& operator<<(std::ostream& o, const VariableColumnCompressor& vcc ) // output
    {
        o << "(";
        for(unsigned j=0; j < vcc.height.size()-1; j++)
        {
            o << vcc.height[j] << ",";
        }
        o << vcc.height[vcc.height.size()-1] << ";";

        for(unsigned j=0; j < vcc.getOutputSize()-1; j++)
        {
            o << vcc.outputs[j] << ",";
        }
        o << vcc.outputs[vcc.outputs.size()-1] << ")";

        return o;
    }


}

	
