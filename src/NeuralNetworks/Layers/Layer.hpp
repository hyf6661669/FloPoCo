#ifndef LAYER_H
#define LAYER_H

#include <string>
#include <vector>

#include "Operator.hpp"
#include "NeuralNetworks/LayerArguments.hpp"
#include "NeuralNetworks/NeuralNetwork.hpp"

using namespace std;


namespace flopoco {

    class NeuralNetwork;

    class Layer : public Operator {

    public:

        Layer(Target* target, LayerArguments* la, NeuralNetwork* parent);
        Layer(Target* target, NeuralNetwork* parent);

        LayerArguments* myArguments;

        int getWidthByPortName(string name);
		
        virtual string getOutputSignalName(int feature);
        virtual string getInputSignalName(int feature);
        virtual string getIntermediateResultName(unsigned int featureNumber){return this->getOutputSignalName(featureNumber)+"_intermediate_result";}

        static inline string getValidDataName(unsigned int featureNumber){return "validData_"+to_string(featureNumber);}
        static inline string getGetNewDataName(unsigned int featureNumber){return "getNewData_"+to_string(featureNumber);}

        static inline string getWeightInputName(){return "weights_i";}
        static inline string getWeightsValidInputName(){return "weightsValid_i";}
        static inline string getLastWeightsInputName(){return "lastWeights_i";}
        static inline string getNewReadAccessOutputName(){return "newReadAccess_o";}
        static inline string getNextStartAddressOutputName(){return "nextStartAddress_o";}
        static inline string getNumberOfBytesOutputName(){return "numberOfBytes_o";}

        bool getInputMemoryParallelAccess() const { return this->inputMemoryParallelAccess; }
        bool getOutputMemoryParallelAccess() const { return this->outputMemoryParallelAccess; }
        int getAddressWidth() const { return this->addressWidth; }

    protected:
        virtual void generateVHDLCode(Target* target);

        bool inputMemoryParallelAccess;
        bool outputMemoryParallelAccess;
        int addressWidth;
        int outputWidth;
        int outputHeight;

        NeuralNetwork* parent;

    private:

        
    };


}//namespace

#endif
