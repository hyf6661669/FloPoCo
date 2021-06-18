//
// Created by Annika Oeste on 17.06.21.
//

#include "GenericAsic.hpp"

namespace flopoco{
    GenericAsic::GenericAsic(): Target()	{
        id_             		= "GenericAsic";
        vendor_         		= "Asic";
    }

    GenericAsic::~GenericAsic() {};

    // TODO: return correct values instead of dummy values in methods
    double GenericAsic::logicDelay(int inputs){ return 0.0;}
    double GenericAsic::adderDelay(int size, bool addRoutingDelay){ return 0.0;}
    double GenericAsic::adder3Delay(int size){ return 0.0;}

    double GenericAsic::eqComparatorDelay(int size){ return 0.0;}
    double GenericAsic::eqConstComparatorDelay(int size){ return 0.0;}


    double GenericAsic::DSPMultiplierDelay(){ return 0.0;}
    double GenericAsic::DSPAdderDelay(){ return 0.0;}
    double GenericAsic::DSPCascadingWireDelay(){ return 0.0;}
    double GenericAsic::DSPToLogicWireDelay(){ return 0.0;}
    double GenericAsic::LogicToDSPWireDelay(){ return 0.0;}
    void   GenericAsic::delayForDSP(MultiplierBlock* multBlock, double currentCp, int& cycleDelay, double& cpDelay) {

    }

    double GenericAsic::RAMDelay() { return 0.0;}
    double GenericAsic::LogicToRAMWireDelay() { return 0.0;}

    double GenericAsic::carryPropagateDelay(){ return 0.0;}
    double GenericAsic::addRoutingDelay(double d) { return 0.0;}
    double GenericAsic::fanoutDelay(int fanout){ return 0.0;}
    double GenericAsic::lutDelay(){ return 0.0;}
    double GenericAsic::ffDelay(){ return 0.0;}

    bool   GenericAsic::suggestSubaddSize(int &x, int wIn) { return false;}
    bool   GenericAsic::suggestSubadd3Size(int &x, int wIn){ return false;}
    bool   GenericAsic::suggestSlackSubaddSize(int &x, int wIn, double slack){ return false;}
    bool   GenericAsic::suggestSlackSubcomparatorSize(int &x, int wIn, double slack, bool constant){ return false;}

    long   GenericAsic::sizeOfMemoryBlock() { return 0;}
    double GenericAsic::lutConsumption(int lutInputSize) { return 0.0;}
}