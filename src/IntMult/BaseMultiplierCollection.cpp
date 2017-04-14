#include "BaseMultiplierCollection.hpp"
#include "BaseMultiplierLUT.hpp"

namespace flopoco {


BaseMultiplierCollection::BaseMultiplierCollection(Target* target){

    srcFileName = "BaseMultiplierCollection";
    uniqueName_ = "BaseMultiplierCollection";
	
    this->target = target;

    //first simple test, shape 0 is a 3x3 mult.:
    BaseMultiplier* bm = new BaseMultiplierLUT(false,false,3,3); //3x3 LUT-based multiplier
    baseMultipliers.push_back(bm);
}

BaseMultiplier* BaseMultiplierCollection::getBaseMultiplier(int shape)
{
    if(shape < baseMultipliers.size())
        return baseMultipliers[shape];
    else
        return nullptr;
}

BaseMultiplierCollection::~BaseMultiplierCollection()
{
    for(BaseMultiplier* bm : baseMultipliers)
    {
        delete bm;
    }
}

	
}   //end namespace flopoco
