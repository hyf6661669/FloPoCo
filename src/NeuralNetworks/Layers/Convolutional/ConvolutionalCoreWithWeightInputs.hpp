//
// Created by nfiege on 18/05/18.
//

#ifndef FLOPOCO_CONVOLUTIONALCOREWITHWEIGHTINPUTS_HPP
#define FLOPOCO_CONVOLUTIONALCOREWITHWEIGHTINPUTS_HPP

#include "Operator.hpp"
#include "utils.hpp"
#include <string>
#include <vector>
#include "NeuralNetworks/Multiplicators/NeuralNetworkMultiplication.hpp"

using namespace std;

namespace flopoco {
    class ConvolutionalCoreWithWeightInputs : public Operator {
    public:
        ConvolutionalCoreWithWeightInputs(Target* target, multiplicationType multType, unsigned int wIn, unsigned int wF, unsigned int wWeightsIn, unsigned int wFWeights, unsigned int size, char roundingType=0x00, string id="0");

        unsigned int getOutputWordSize() const {return this->outputWordSize;}
        unsigned int getOutputFraction() const {return this->outputFraction;}

    private:
        unsigned int outputWordSize;
        unsigned int outputFraction;
        unsigned int wordSize;
        unsigned int fraction;
        unsigned int weightWordSize;
        unsigned int weightFraction;

        void roundOutput(Target* target, unsigned int wordSizeFrom, unsigned int fractionFrom, string roundFromName, char roundingType);
    };

}
#endif //FLOPOCO_CONVOLUTIONALCOREWITHWEIGHTINPUTS_HPP
