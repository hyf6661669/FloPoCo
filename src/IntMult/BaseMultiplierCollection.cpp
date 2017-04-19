#include "BaseMultiplierCollection.hpp"
#include "BaseMultiplierLUT.hpp"
#include "BaseMultiplier2xk.hpp"

namespace flopoco {


BaseMultiplierCollection::BaseMultiplierCollection(Target* target){

    srcFileName = "BaseMultiplierCollection";
    uniqueName_ = "BaseMultiplierCollection";
	
    this->target = target;

    //first simple test, shape 0 is a 3x3 mult.:
    baseMultipliers.push_back(new BaseMultiplierLUT(false,false,3,3)); //3x3 LUT-based multiplier

    baseMultipliers.push_back(new BaseMultiplier2xk(false, false, 10)); //2x10 LUT/carry-chain-based multiplier
}

BaseMultiplier* BaseMultiplierCollection::getBaseMultiplier(int shape)
{
    if(shape < ((int) baseMultipliers.size()))
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
