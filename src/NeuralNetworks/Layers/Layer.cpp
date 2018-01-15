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

Layer::Layer(Target *target, LayerArguments* la) : Operator(target)
{
//    this->wordSize=-1;
//    this->fraction=-1;
//    this->weightWordSize=-1;
//    this->weightFraction=-1;

//    this->horizontalSize=-1;
//    this->verticalSize=-1;
//    this->numberOfInputFeatures=-1;
//    this->numberOfOutputFeatures=-1;

//    this->inputFeaturesParallel=false;;
//    this->outputFeaturesParallel=false;

//    this->activationFunction="No Activation Function";

//    this->layerType="No Type";

    this->myArguments=la;
}

Layer::Layer(Target *target) : Operator(target)
{
    this->myArguments = new LayerArguments();
}

string Layer::getOutputSignalName(int feature)
{
    cout << "Layer.getOutputSignalName: this should never be called" << endl;
    return "";
}

string flopoco::Layer::getInputSignalName(int feature)
{
    cout << "Layer.getInputSignalName: this should never be called" << endl;
    return "";
}

void Layer::generateVHDLCode(Target* target)
{
    cout << "Layer.generateVHDLCode: this should never be called" << endl;
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
