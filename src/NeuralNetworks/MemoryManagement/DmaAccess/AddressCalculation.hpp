//
// Created by nfiege on 10/12/18.
//

#ifndef FLOPOCO_ADDRESSCALCULATION_HPP
#define FLOPOCO_ADDRESSCALCULATION_HPP

#include "Operator.hpp"

/* This file contains a lot of useful functions to manipulate vhdl */
#include "utils.hpp"
#include <string>

namespace flopoco {
    class AddressCalculation : public Operator{

    public:
        /*!
         * @brief create circuit to calculate memory address: addr = startAddress + (innerCounter * #bytesPerAccess) + (outerCounter * #bytesPerFeature)
         * @param target
         * @param startAddress base address for the layer
         * @param innerCounterWidth
         * @param outerCounterWidth
         * @param bytesPerAccess the number of bytes for the given input feature:
         * a negative number indicates that this is NOT A CONSTANT and instead specifies the word size of the input
         * this is the case for dense layers, where the number of bytes per access is not constant
         * @param bytesPerFeature the number of bytes per output feature
         */
        AddressCalculation(Target *target, unsigned int startAddress, unsigned int innerCounterWidth,
                           unsigned int outerCounterWidth, int bytesPerAccess, int bytesPerFeature);

        static std::string uintTo32BitString(unsigned int i);
    };
}


#endif //FLOPOCO_ADDRESSCALCULATION_HPP
