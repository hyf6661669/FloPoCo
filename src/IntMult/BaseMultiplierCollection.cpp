#include "BaseMultiplierCollection.hpp"
#include "BaseMultiplierLUT.hpp"
#include "BaseMultiplier2xk.hpp"

namespace flopoco {


BaseMultiplierCollection::BaseMultiplierCollection(Target* target){

    srcFileName = "BaseMultiplierCollection";
    uniqueName_ = "BaseMultiplierCollection";
	
    this->target = target;

    //create logic-based multipliers:
    baseMultipliers.push_back(new BaseMultiplierLUT(false,false,1,1)); //1x1 LUT-based multiplier (an AND gate)
    baseMultipliers.push_back(new BaseMultiplierLUT(false,false,1,2)); //1x2 LUT-based multiplier (two AND gates)
    baseMultipliers.push_back(new BaseMultiplierLUT(false,false,2,1)); //2x1 LUT-based multiplier (two AND gates)
    baseMultipliers.push_back(new BaseMultiplierLUT(false,false,2,3)); //2x3 LUT-based multiplier (three LUT6)
    baseMultipliers.push_back(new BaseMultiplierLUT(false,false,3,2)); //3x2 LUT-based multiplier (three LUT6)
    baseMultipliers.push_back(new BaseMultiplierLUT(false,false,3,3)); //3x3 LUT-based multiplier (six LUT6)

    for(int k=2; k < 32; k++) //ToDo: adjust limits
    {
        baseMultipliers.push_back(new BaseMultiplier2xk(false, false, k, false)); //2xk LUT/carry-chain-based multiplier
        baseMultipliers.push_back(new BaseMultiplier2xk(false, false, k, true));  //kx2 LUT/carry-chain-based multiplier
    }

    int i=0;
    cout << "The following multiplier shapes were generated:" << endl;

    for(BaseMultiplier* bm : baseMultipliers)
    {
        cout << "shape " << i++ << ": " << bm->getXWordSize() << "x" << bm->getYWordSize() << " of type " << bm->getName() << endl;
    }
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
