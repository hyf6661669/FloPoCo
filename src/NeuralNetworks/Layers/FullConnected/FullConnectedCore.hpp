//
// Created by nfiege on 22/06/18.
//

#ifndef FLOPOCO_FULLCONNECTEDCORE_HPP
#define FLOPOCO_FULLCONNECTEDCORE_HPP

#include "Operator.hpp"

#include "utils.hpp"
#include <string>
#include "NeuralNetworks/Multiplicators/NeuralNetworkMultiplication.hpp"

namespace flopoco {

    class FullConnectedCore : public Operator {

    public:
        FullConnectedCore(Target* target, multiplicationType multType, int wordSize, int fraction, int weightWordSize, int weightFraction, int secondInputWidth, char roundingMode=0x00);

        ~FullConnectedCore() {}

        int getOutputWordSize() const { return this->outputWordSize; }

    private:
        int outputWordSize;
    };
} // namespace flopoco

#endif //FLOPOCO_FULLCONNECTEDCORE_HPP
