#ifndef Xilinx_TernaryAdd_2State_H
#define Xilinx_TernaryAdd_2State_H

#include "Operator.hpp"
#include "utils.hpp"

#include "Xilinx_LUT_compute.h"


namespace flopoco {
    class Xilinx_TernaryAdd_2State : public Operator {
      public:
        short mapping[3];
        short state_type;
        bool negating;

        Xilinx_TernaryAdd_2State( Target *target, int wIn, short state, short state2 = -1 );
        ~Xilinx_TernaryAdd_2State() {}

        void insertCarryInit( short state, short state2 );
        void computeState( short state, short state2 );
        string computeLUT( short state, short state2 );
	};
}//namespace

#endif
