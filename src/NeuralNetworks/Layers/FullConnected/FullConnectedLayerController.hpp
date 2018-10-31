//
// Created by nfiege on 25/06/18.
//

#ifndef FLOPOCO_FULLCONNECTEDLAYERCONTROLLER_HPP
#define FLOPOCO_FULLCONNECTEDLAYERCONTROLLER_HPP

#include "Operator.hpp"

#include "utils.hpp"
#include <string>

namespace flopoco {


    class FullConnectedLayerController : public Operator {
    public:
        FullConnectedLayerController(Target* target, unsigned int numberOfNeurons, unsigned int weightsPerNeuron, unsigned int weightRAMDepth, unsigned int fullConnectedCorePipelineDepth, int fullConnectedCoreAdderPipelineDepth=-1);

        ~FullConnectedLayerController() {}

        bool thisHasMultiplexerSwitchOutput() const {return this->hasMultiplexerSwitchOutput;}

    private:
        bool hasMultiplexerSwitchOutput;
    };
}

#endif //FLOPOCO_FULLCONNECTEDLAYERCONTROLLER_HPP
