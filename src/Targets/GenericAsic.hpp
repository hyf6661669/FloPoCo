//
// Created by Annika Oeste on 17.06.21.
//

#include "../Target.hpp"

namespace flopoco{
    class GenericAsic: public Target {
    public:
        /** The default constructor. */
        GenericAsic();
        /** The destructor */
        ~GenericAsic();

        /** overloading the virtual functions of Target
		 * @see the target class for more details
		 */

        double logicDelay(int inputs);

        double adderDelay(int size, bool addRoutingDelay=true);

        double adder3Delay(int size);
        double eqComparatorDelay(int size);
        double eqConstComparatorDelay(int size);


        double DSPMultiplierDelay();
        double DSPAdderDelay();
        double DSPCascadingWireDelay();
        double DSPToLogicWireDelay();
        double LogicToDSPWireDelay();
        void   delayForDSP(MultiplierBlock* multBlock, double currentCp, int& cycleDelay, double& cpDelay);

        double RAMDelay();
        double LogicToRAMWireDelay();

        double carryPropagateDelay();
        double addRoutingDelay(double d);
        double fanoutDelay(int fanout = 1);
        double lutDelay();
        double ffDelay();

        bool   suggestSubmultSize(int &x, int &y, int wInX, int wInY);
        bool   suggestSubaddSize(int &x, int wIn);
        bool   suggestSubadd3Size(int &x, int wIn);
        bool   suggestSlackSubaddSize(int &x, int wIn, double slack);
        bool   suggestSlackSubadd3Size(int &x, int wIn, double slack);
        bool   suggestSlackSubcomparatorSize(int &x, int wIn, double slack, bool constant);

        long   sizeOfMemoryBlock();
        DSP*   createDSP();
        int    getEquivalenceSliceDSP();
        double lutConsumption(int lutInputSize);
    };
}

