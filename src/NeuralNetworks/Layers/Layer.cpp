// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of this class
#include "Layer.hpp"


using namespace std;
namespace flopoco {

Layer::Layer(Target *target, LayerArguments* la, NeuralNetwork* parent_) : Operator(target), parent(parent_)
{
    this->myArguments=la;
    this->inputMemoryParallelAccess=false;
    this->outputMemoryParallelAccess=false;
    this->addressWidth = 0;
    this->outputWidth = 0;
    this->outputHeight = 0;
}

Layer::Layer(Target *target, NeuralNetwork* parent_) : Operator(target), parent(parent_)
{
    this->myArguments = new LayerArguments();
    this->inputMemoryParallelAccess=false;
    this->outputMemoryParallelAccess=false;
    this->addressWidth=0;
    this->outputHeight=0;
    this->outputWidth=0;
}

string Layer::getOutputSignalName(int feature)
{
    stringstream e;
    e << "Layer.getOutputSignalName: this should never be called" << endl;
    THROWERROR(e.str());
    return "";
}

string flopoco::Layer::getInputSignalName(int feature)
{
    stringstream e;
    e << "Layer.getInputSignalName: this should never be called" << endl;
    THROWERROR(e.str());
    return "";
}

void Layer::generateVHDLCode(Target* target)
{
    stringstream e;
    e << "Layer.generateVHDLCode: this should never be called" << endl;
    THROWERROR(e.str());
}

int flopoco::Layer::getWidthByPortName(string name)
{
    vector<flopoco::Signal*>* ioL = this->getIOList();
    for(auto it : (*ioL))
    {
        if(it->getName()==name)
        {
            return it->width();
        }
    }
    stringstream e;
    e << "Port with name '" << name << "' does not exist in this Layer!";
    THROWERROR(e.str());
    return -1;
}



}//namespace flopoco
