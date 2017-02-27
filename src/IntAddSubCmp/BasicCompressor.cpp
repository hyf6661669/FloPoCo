#include <iostream>
#include <sstream>
#include <string>
#include "gmp.h"
#include "mpfr.h"
#include <vector>
#include <gmpxx.h>
#include <stdio.h>
#include <stdlib.h>
#include "BasicCompressor.hpp"


using namespace std;
namespace flopoco{

// personalized parameter
//string BasicCompressor::operatorInfo = "UserDefinedInfo list parameter;

/* Uni KS start */
BasicCompressor::BasicCompressor(Target * target)
:Operator(target)
{
	areaCost = 1.0;
}
/* Uni KS stop */

BasicCompressor::BasicCompressor(Target * target, vector<int> h)
:Operator(target)
{
	ostringstream name;
	stringstream nm;
	/// compressors are always combinational independently by the target technology
	setCombinatorial();

	int w=0;
	int param=0;

	while(h[h.size()-1]==0)
	{
		h.erase(h.end()-1);
	}

	for(int i=h.size()-1; i>=0;i--)
		height.push_back(h[i]);

	name << "Compressor_";

	for(unsigned i=0; i<height.size();i++)
	{
		w += height[i];
		param=param+intpow2(height.size()-i-1)*height[i];
		name<<height[i];
	}

	wOut=intlog2(param);
	/* Uni KS start */
	outputs.resize(wOut);
	for(int i=0; i < outputs.size(); i++) outputs[i]=1;

	//set area cost of compressor:
	string targetID = target->getID();

	if((targetID == "Virtex5") || (targetID == "Virtex6") || (targetID == "Virtex7") || (targetID == "Spartan6"))
	{
		//count a LUT5 with shared inputs as half of the cost of a LUT6:
		if(w < target->lutInputs())
		{
			areaCost = ceil((double) wOut*0.5); //less LUT inputs are used and two outputs can be computed per LUT
		}
		else
		{
			areaCost = (double) wOut; //all LUT inputs are used and one output is computed per LUT
		}
	}
	else
	{
		areaCost = (double) wOut; //all LUT inputs are used and one output is computed per LUT
	}	
	/* Uni KS stop */

	name << "_" << wOut;
	setName(name.str());
	setCopyrightString("Bogdan Popa, Illyes Kinga, 2012");



	stringstream xs;


	for(unsigned i=0;i<height.size();i++)
	{
		addInput(join("X",i), h[i]);

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

	addOutput("R", wOut);

	vhdl << tab << declare("X", w) << " <=" << xs.str();

	vhdl << tab << "with X select R <= \n";

	for (mpz_class i = 0; i < (1 << w); i++)
	{

		mpz_class ppcnt=0;//popcnt(i);
		mpz_class ii=i;
		for(unsigned j=0;j<h.size();j++)
		{


			ppcnt+=popcnt(ii-((ii>>h[j])<<h[j]))*intpow2(j);
			ii=ii>>h[j];
		}

		vhdl << tab << tab << "\"" << unsignedBinary(ppcnt,wOut) << "\" when \""
		<< unsignedBinary(i,w) << "\", \n";


		}



		vhdl << tab << tab << "\"" << std::string (wOut, '-') << "\" when others;\n" << endl;

		REPORT(DEBUG, "Generated " << name.str());


	}




	BasicCompressor::~BasicCompressor(){
	}

	unsigned BasicCompressor::getColumnSize(int column)
	{
		if (column>=(signed)height.size())
			return 0;
		else
			return height[height.size()-column-1];
	}

    unsigned BasicCompressor::getOutputSize() const
	{
		return wOut;
	}

	/* Uni KS start */
	unsigned BasicCompressor::getNumberOfColumns()
	{
		return height.size();
	}
	/* Uni KS stop */

	void BasicCompressor::emulate(TestCase * tc)
	{
		mpz_class r=0;

		for(unsigned i=0;i<height.size();i++)
			{
				mpz_class sx = tc->getInputValue(join("X",i));
				mpz_class p= popcnt(sx);
				r += p<<i;
			}



		tc->addExpectedOutput("R", r);
	}

    std::ostream& operator<<(std::ostream& o, const BasicCompressor& bc ) // output
    {
        o << "(";
        for(unsigned j=0; j < bc.height.size()-1; j++)
        {
            o << bc.height[j] << ",";
        }
        o << bc.height[bc.height.size()-1] << ";" << bc.getOutputSize() << ")";
        return o;
    }

//	OperatorPtr BasicCompressor::parseArguments(Target *target, const vector<string> &args) {
//		return new BasicCompressor(target);
//	}
//
//	void BasicCompressor::registerFactory(){
//		UserInterface::add("BasicCompressor", // name
//											 "",
//											 "operator; floating point; floating-point adders", // categories
//											 "",
//											 "",
//											 BasicCompressor::parseArguments
//											 ) ;
//
//	}
}




