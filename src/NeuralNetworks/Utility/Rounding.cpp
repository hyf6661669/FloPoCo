// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

//for reading the .txt file
#include <fstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "Rounding.hpp"

#include "PrimitiveComponents/GenericMux.hpp"
#include "NeuralNetworks/NeuralNetwork.hpp"

using namespace std;
namespace flopoco {




    Rounding::Rounding(Target* target, unsigned int wordSizeFrom, unsigned int fractionFrom, unsigned int wordSizeTo, unsigned int fractionTo, roundingTypeEnum round) : Operator(target) {

        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="Rounding";

        // author
        setCopyrightString("Nicolai Fiege, 2018");

        // definition of the name of the operator
        ostringstream name;
        name << "Rounding_" << wordSizeFrom << "_" << fractionFrom << "_to_" << wordSizeTo << "_" << fractionTo << "_";
        switch(round){
            case  roundingTypeEnum::truncation: name << "trunc";  break;
            case  roundingTypeEnum::saturation: name << "satur"; break;
            default: THROWERROR("Undefined rounding type");
        }
        setName(name.str());

        // std logic arith is THE DEVIL
        useNumericStd();

        //add in/out
        addInput("X",wordSizeFrom);
        addOutput("R",wordSizeTo);
        switch(round){
            case  roundingTypeEnum::truncation: this->buildTruncation(target,wordSizeFrom,fractionFrom,wordSizeTo,fractionTo);  break;
            case  roundingTypeEnum::saturation: this->buildSaturation(target,wordSizeFrom,fractionFrom,wordSizeTo,fractionTo); break;
            default: THROWERROR("Undefined rounding type");
        }
    }

    void Rounding::buildTruncation(Target* target, unsigned int wordSizeFrom, unsigned int fractionFrom, unsigned int wordSizeTo, unsigned int fractionTo)
    {
        int MSBFrom = wordSizeFrom - fractionFrom;
        int MSBTo = wordSizeTo-fractionTo;
        int from = wordSizeFrom-(MSBFrom-MSBTo)-1;
        int downto = fractionFrom-fractionTo;
        if(downto>0)
        {
            this->vhdl << tab << declare("roundingNumber",downto) << " <= X(" << downto-1 << ")";
            if(downto>1)
            {
                this->vhdl << " & (" << downto-2 << " downto 0 => '0')";
            }
            this->vhdl <<";" << endl;
            this->vhdl << tab << declare("R_beforeRound",wordSizeFrom) << " <= std_logic_vector(signed(X)+"
                       << "signed(roundingNumber));" << endl;


            this->vhdl << tab << "R <= R_beforeRound(R_beforeRound'length-1) & R_beforeRound(" << from-1 << " downto " << downto << ");" << endl;
        }
        else if(downto==0)
        {
            this->vhdl << tab << "R <= X(X'length-1) & X(" << from-1 << " downto " << downto << ");" << endl;
        }
        else
        {
            stringstream e;
            e << "Error while rounding, requested fraction is bigger than actual fraction! This should normally not happen";
            THROWERROR(e.str());
        }
    }

    void Rounding::buildSaturation(Target* target, unsigned int wordSizeFrom, unsigned int fractionFrom, unsigned int wordSizeTo, unsigned int fractionTo)
    {
        mpz_class xMax_d = mpz_class(1);
        xMax_d = xMax_d << (wordSizeTo-1);
        mpz_class xMin_d = mpz_class(-1) * xMax_d;

        xMax_d = xMax_d - 1;


        Rounding* truncOp = new Rounding(target,wordSizeFrom,fractionFrom,wordSizeTo,fractionTo,roundingTypeEnum::truncation);
        addSubComponent(truncOp);
        inPortMap(truncOp,"X","X");
        outPortMap(truncOp,"R","X_truncated",true);
        this->vhdl << instance(truncOp,"Truncation_instance");

        this->vhdl << declare("X_min",wordSizeTo) << " <= \"" << NeuralNetwork::convertMpzToString(mpz_class(xMin_d),wordSizeTo) << "\";" << endl;
        this->vhdl << declare("X_max",wordSizeTo) << " <= \"" << NeuralNetwork::convertMpzToString(mpz_class(xMax_d),wordSizeTo) << "\";" << endl;

        this->vhdl << declare("X_smaller_min",1) << " <= \"1\" when signed(X(X'length-1 downto " << (fractionFrom-fractionTo) << ")) < signed(X_min) else \"0\";" << endl;
        this->vhdl << declare("X_bigger_max",1) << " <= \"1\" when signed(X(X'length-1 downto " << (fractionFrom-fractionTo) << ")) > signed(X_max) else \"0\";" << endl;

        GenericMux* minMux = new GenericMux(target,wordSizeTo,2);
        addSubComponent(minMux);
        inPortMap(minMux,minMux->getSelectName(),"X_smaller_min");
        inPortMap(minMux,minMux->getInputName(0),"X_truncated");
        inPortMap(minMux,minMux->getInputName(1),"X_min");
        outPortMap(minMux,minMux->getOutputName(),"X_temp",true);
        this->vhdl << instance(minMux,"Min_Mux");

        GenericMux* maxMux = new GenericMux(target,wordSizeTo,2);
        addSubComponent(maxMux);
        inPortMap(maxMux,maxMux->getSelectName(),"X_bigger_max");
        inPortMap(maxMux,maxMux->getInputName(0),"X_temp");
        inPortMap(maxMux,maxMux->getInputName(1),"X_max");
        outPortMap(maxMux,maxMux->getOutputName(),"R",false);
        this->vhdl << instance(maxMux,"Max_Mux");
    }

}//namespace flopoco
