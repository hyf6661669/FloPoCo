#ifndef LAYER_H
#define LAYER_H

#include <string>
#include <vector>

#include "Operator.hpp"
#include "NeuralNetworks/LayerArguments.hpp"

using namespace std;


namespace flopoco {

    class Layer : public Operator {

    public:

        Layer(Target* target, LayerArguments* la);
        Layer(Target* target);

//        int wordSize;
//		int fraction;
//		int weightWordSize;
//		int weightFraction;

//        int horizontalSize;
//        int verticalSize;
//        int numberOfInputFeatures;
//        int numberOfOutputFeatures;

//        bool inputFeaturesParallel;
//        bool outputFeaturesParallel;

//        string activationFunction;

//        string layerType;

        LayerArguments* myArguments;
		
        virtual string getOutputSignalName(int feature);
        virtual string getInputSignalName(int feature);

    protected:
        virtual void generateVHDLCode(Target* target);

    private:

        
    };


}//namespace

#endif
