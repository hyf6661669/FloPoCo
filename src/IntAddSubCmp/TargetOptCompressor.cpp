#include <iostream>
#include <sstream>
#include <string>
#include "gmp.h"
#include "mpfr.h"
#include <vector>
#include <gmpxx.h>
#include <stdio.h>
#include <stdlib.h>
#include "TargetOptCompressor.hpp"


using namespace std;
namespace flopoco{

	TargetOptCompressor::TargetOptCompressor(Target * target, vector<int> h)
	:BasicCompressor(target)
	{
		int w=0;
		int param=0;

		setSequential(); //a pipelined compressors is sequential

		while(h[h.size()-1]==0)
		{
			h.erase(h.end()-1); 		//what happens here?
		}

		for(int i=h.size()-1; i>=0;i--)		//what happens here?
			height.push_back(h[i]);

		ostringstream name;
		name << "Compressor_";

		ostringstream compname;
		compname << "gpc";

		for(unsigned i=0; i<height.size();i++)
		{
			w=w+height[i];
			param=param+intpow2(height.size()-i-1)*height[i];
			if(height[height.size()-i-1] > 0)
			{
				stringstream inputName;
				inputName << "X" << i;
				addInput(inputName.str(),height[height.size()-i-1]);
			}
			name << height[i];
			compname << '_' << height[i];
		}
		wOut=intlog2(param);
		outputs.resize(wOut);
		for(unsigned i=0; i < outputs.size(); i++) outputs[i]=1;

		name << "_" << wOut;
		setName(name.str());

		setCopyrightString("Martin Kumm, 2014");

		stringstream outputName;
		outputName << "R";
		addOutput(outputName.str(),wOut);

		compname << '_' << wOut;

		//this component currently simply wrapps to the GPC library
		vhdl << endl;
        vhdl << tab << declare("ce_i_const") << " <= \'1\'; --workaround for constants. Vivado does not accept ce_i => \'1\' in portmap" << endl;
        //vhdl << tab << "ce_i_const <= '1'" << endl; //workaround for constants. Vivado does not accept ce_i => '1' in portmap
		vhdl << tab << "gpc_inst : entity work." << compname.str() << endl;
		vhdl << tab << tab << "generic map (use_output_ff => false)" << endl;
		vhdl << tab << tab << "port map (clk_i => clk, rst_i => rst, ce_i => ce_i_const, ";

		char var = 'u';
		for(unsigned i=0; i<height.size();i++)
		{
			if(height[height.size()-i-1] > 0)
			{
				vhdl << var << "_i => X" << i << ", ";
			}
			var++;
		}
		vhdl << "s_o => R);" << endl;
		vhdl << endl;
	}
	
	TargetOptCompressor::~TargetOptCompressor()
	{
	}

/*
	unsigned TargetOptCompressor::getOutputSize()
	{
		return outputheight.size();
	}
*/
}
	


	
