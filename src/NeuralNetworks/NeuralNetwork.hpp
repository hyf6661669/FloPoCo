#ifndef NEURALNETWORK_H
#define NEURALNETWORK_H
/* Each Operator declared within the flopoco framework has 
   to inherit the class Operator and overload some functions listed below*/
#include "Operator.hpp"

/* This file contains a lot of useful functions to manipulate vhdl */
#include "utils.hpp"
#include <string>
#include <vector>

#include "NeuralNetworks/Layers/Layer.hpp"

#include "LayerArguments.hpp"

/*  All flopoco operators and utility functions are declared within
    the flopoco namespace.
    You have to use flopoco:: or using namespace flopoco in order to access these
    functions.
*/

namespace flopoco {

	// new operator class declaration
    class NeuralNetwork : public Operator {

    public:
        NeuralNetwork(Target* target, map <string, LayerArguments*> layers_, vector <pair <string, string>> edges_, string name_ = "NeuralNetwork");
        NeuralNetwork(Target* target, string name_ = "NeuralNetwork");

		// destructor
        ~NeuralNetwork() {}


        void addLayer(string uniqueLayerName, LayerArguments* arguments);
        void addEdge(string from, string to);
        void buildNeuralNetwork(Target* target);

    private:
        string name;
        map <string, LayerArguments*> layerArgs; //map: layer name -> layer arguments
        //map <string, Layer*> layers; //map: layer name -> layer arguments
        vector <pair <string, string>> edges; //pair: layer1 name, layer2 name: connection from layer1 to layer2

        map <string, Layer*> buildAllLayers(Target* target, Operator* globalOp);
        Layer* buildSpecificLayer(Target* target, string nameOfThisLayer);

        void buildAllEdges(Target* target, map<string, Layer*> &layerPtrs);
        void buildSpecificEdge(Target* target, string from, string to, Layer* fromOp, Layer* toOp);
        void checkEdges();
        void checkConnection(Layer* from, Layer* to);
        bool checkIfLayerExists(string layerName);

        Layer* buildConvLayer(Target* target, string convLayer);
        Layer* buildFullConnectedLayer(Target* target, string FCLayer);
        Layer* buildPoolingLayer(Target* target, string poolingLayer);

        Layer* buildInput(Target* target, string nextLayerName, Layer* nextLayer);
        Layer* buildOutput(Target* target, string lastLayerName, Layer* lastLayer);
        Operator* buildGlobalController(Target* target);

        void instEverything(map<string, Layer*> &instMe, Operator* globalOp);
	};


}//namespace

#endif
