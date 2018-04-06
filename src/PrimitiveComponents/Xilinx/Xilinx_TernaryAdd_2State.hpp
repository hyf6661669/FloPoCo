#ifndef Xilinx_TernaryAdd_2State_H
#define Xilinx_TernaryAdd_2State_H

#include "Operator.hpp"
#include "utils.hpp"

#include "Xilinx_LUT_compute.h"


namespace flopoco {
    class Xilinx_TernaryAdd_2State : public Operator {
        int wIn_;
        short state_,state2_;
      public:
        short mapping[3];
        short state_type;

        Xilinx_TernaryAdd_2State( Target *target,const int &wIn,const short &state,const short &state2 = -1 );
        ~Xilinx_TernaryAdd_2State() {}

        void insertCarryInit( );
        void computeState( );
        string computeLUT( );

        static OperatorPtr parseArguments( Target *target, vector<string> &args ){
            if( target->getVendor() != "Xilinx" )
                throw std::runtime_error( "Can't build xilinx primitive on non xilinx target" );

            int wIn;
            int mode,mode2;
            UserInterface::parseInt(args,"wIn",&wIn );
            UserInterface::parseInt(args,"mode",&mode);
            UserInterface::parseInt(args,"mode2",&mode2);
            return new Xilinx_TernaryAdd_2State(target,wIn,mode,mode2);
        }

        static void registerFactory(){
            UserInterface::add( "xilinx_ternary_addsub", // name
                                "A ternary adder/subtractor build of xilinx primitives. For each configuration (selectable by optional input sel_i) the signs of up to two inputs may be different.", // description, string
                                "Primitives", // category, from the list defined in UserInterface.cpp
                                "",
                                "wIn(int): The wordsize of the adder; \
                                mode(int)=0: Input negation bitmask for 1st configuration, x_i is negated when (mode & 1)==1, y_i is negated when (mode & 2)==1, z_i is negated when (mode & 4)==1; \
                                mode2(int)=-1: Input negation bitmask for 2nd configuration, x_i is negated when (mode2 & 1)==1, y_i is negated when (mode2 & 2)==1, z_i is negated when (mode2 & 4)==1",
                                "",
                                Xilinx_TernaryAdd_2State::parseArguments
                              );
        }

        // Operator interface
    public:
        virtual void emulate(TestCase *tc)
        {
            mpz_class x = tc->getInputValue("x_i");
            mpz_class y = tc->getInputValue("y_i");
            mpz_class z = tc->getInputValue("z_i");
            mpz_class s = 0;
            mpz_class sel = 0;
            if( state_ != state2_  ){
                sel = tc->getInputValue("sel_i");
            }
            short signs = 0;
            if(sel==0){
                signs = state_;
            }else{
                signs = state2_;
            }

            if(0x1&signs)
                s -= x;
            else
                s += x;

            if(0x2&signs)
                s -= y;
            else
                s += y;

            if(0x4&signs)
                s -= z;
            else
                s += z;

            mpz_class mask = ((1<<wIn_)-1);
            s = s & mask;
            tc->addExpectedOutput("sum_o", s);
        }
	};
}//namespace

#endif
