//
// Created by nfiege on 28/06/18.
//

#ifndef FLOPOCO_NEURALNETWORKMULTIPLICATION_HPP
#define FLOPOCO_NEURALNETWORKMULTIPLICATION_HPP

#include "Operator.hpp"

#include "utils.hpp"
#include <string>
namespace flopoco {
    typedef enum {
        simple,
        none
    } multiplicationType;


    class NeuralNetworkMultiplication : public Operator {

    public:
        NeuralNetworkMultiplication(Target *target, int wIn1_, int wIn2_, multiplicationType mType=multiplicationType::simple);

        int getOutputWordSize() const {return this->wOut;}

    private:
        int wIn1;
        int wIn2;
        int wOut;

        void buildSimpleMultiplier(Target* target);
    };
}

#endif //FLOPOCO_NEURALNETWORKMULTIPLICATION_HPP
