#include "BaseMultiplierCollection.hpp"


namespace flopoco {


BaseMultiplierCollection::BaseMultiplierCollection(Target* target){

    srcFileName = "BaseMultiplierCollection";
    uniqueName_ = "BaseMultiplierCollection";
	
	target_ = target;

    BaseMultiplier* bm;// = new
    baseMultipliers.push_back(bm);
}

BaseMultiplierCollection::~BaseMultiplierCollection()
{
    for(BaseMultiplier* bm : baseMultipliers)
    {
        delete bm;
    }
}

	
}   //end namespace flopoco
