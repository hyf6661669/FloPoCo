// general c++ library for manipulating streams
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"
#include "Xilinx_Primitive.hpp"

using namespace std;
namespace flopoco {
    Xilinx_Primitive::Xilinx_Primitive(Target* target) : Operator(target) {
        stringstream copyr;
        copyr << ">>Universität Kassel" << endl;
        copyr << ">>Fachgebiet Digitaltechnik" << endl;
        copyr << ">>Marco Kleinlein" << endl;
        setCopyrightString(copyr.str());
        setCombinatorial();
    }

    Xilinx_Primitive::~Xilinx_Primitive() {}

    void Xilinx_Primitive::setGeneric( string name, string value){
        generics_.insert( std::make_pair(name,value) );
    }

    void Xilinx_Primitive::setGeneric(string name, const long value){
        setGeneric(name,std::to_string(value) );
    }

    std::map<string, string> &Xilinx_Primitive::getGenerics(){
        return generics_;
    }

    void Xilinx_Primitive::outputVHDL(ostream &o, string name){
        o << "library UNISIM;";
        o << "use UNISIM.Vcomponents.all;";
    }

    void Xilinx_Primitive::outputVHDLComponent(ostream &o, string name){
    }

    string Xilinx_Primitive::primitiveInstance(string instanceName){
        Xilinx_Primitive *op = this;
        ostringstream o;
        // TODO add checks here? Check that all the signals are covered for instance

        o << tab << instanceName << ": " << op->getName();
        if (op->isSequential())
            o << "  -- pipelineDepth="<< op->getPipelineDepth() << " maxInDelay=" << getMaxInputDelays(op->getInputDelayMap());
        o << endl;

        if( !getGenerics().empty() ){
            o << tab << tab << "generic map ( ";
            std::map<string,string>::iterator it = getGenerics().begin();
            o << it->first << " => " << it->second;
            for( ++it; it!=getGenerics().end();++it  ){
                o << "," << endl << tab << tab << it->first << " => " << it->second;

            }
            o << ")" << endl;
        }

        o << tab << tab << "port map ( ";
        // build vhdl and erase portMap_
        map<string,string>::iterator it;
        if(op->isSequential()) {
            o << "clk  => clk";
            o << "," << endl << tab << tab << "           rst  => rst";
            if (op->isRecirculatory()) {
                o << "," << endl << tab << tab << "           stall_s => stall_s";
            };
            if (op->hasClockEnable()) {
                o << "," << endl << tab << tab << "           ce => ce";
            };
        }


        for (it=op->getPortMap().begin()  ; it != op->getPortMap().end(); it++ ) {
            bool outputSignal = false;
            for ( int k = 0; k < int(op->getIOList()->size()); k++){
                if ((op->getIOList()->at(k)->type() == Signal::out) && ( op->getIOList()->at(k)->getName() == (*it).first )){
                    outputSignal = true;
                }
            }

            bool parsing = vhdl.isParsing();

            if ( outputSignal && parsing){
                vhdl.flush(getCurrentCycle());
                vhdl.disableParsing(true);
            }

            if (it!=op->getPortMap().begin() || op->isSequential())
                o << "," << endl <<  tab << tab << "           ";

            // The following code assumes that the IO is declared as standard_logic_vector
            // If the actual parameter is a signed or unsigned, we want to automatically convert it
            Signal* rhs;
            string rhsString;
            // The following try was intended to distinguish between variable and constant
            // but getSignalByName doesn't catch delayed variables
            try{
                //cout << "its = " << (*it).second << "  " << endl;
                rhs = getDelayedSignalByName((*it).second);
                if (rhs->isFix() && !outputSignal){
                    rhsString = std_logic_vector((*it).second);
                }
                else {
                    rhsString = (*it).second;
                }

            }
            catch(string e) {
                //constant here
                rhsString=(*it).second;
            }

            o << (*it).first << " => " << rhsString;

            if ( outputSignal && parsing ){
                vhdl << o.str();

                vhdl.flush(getCurrentCycle());
                o.str("");
                vhdl.disableParsing(!parsing);
            }
            //op->portMap_.erase(it);
        }

        o << ");" << endl;



/*
        //Floorplanning related-----------------------------------------
        floorplan << manageFloorplan();
        flpHelper->addToFlpComponentList(op->getName());
        flpHelper->addToInstanceNames(op->getName(), instanceName);
        //--------------------------------------------------------------

*/
        // add the operator to the subcomponent list: still explicit in origin/master
        //
        // subComponents_.push_back(op);
        return o.str();
    }
}//namespace
