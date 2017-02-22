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
        vhdl << "--VariableColumnCompressor" << endl;
	}
	
	VariableColumnCompressor::~VariableColumnCompressor()
	{
	}
	
}
	


	
