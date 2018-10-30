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

#include "NeuralNetworks/Multiplicators/NeuralNetworkMultiplication.hpp"

/*  All flopoco operators and utility functions are declared within
    the flopoco namespace.
    You have to use flopoco:: or using namespace flopoco in order to access these
    functions.
*/

namespace flopoco {

	class Layer;


	// new operator class declaration
    class NeuralNetwork : public Operator {

    public:
        NeuralNetwork(Target* target, map <string, LayerArguments*> layers_, vector <pair <string, string>> edges_, string name_ = "NeuralNetwork", bool serial = false, char roundingType_=0x00, char multiplicationTypeChar_=0x00, bool onlyOutputClassIndex_=false);
        NeuralNetwork(Target* target, string name_ = "NeuralNetwork", char roundingType_=0x00, char multiplicationTypeChar_=0x00, bool onlyOutputClassIndex_=false);

		// destructor
        ~NeuralNetwork() {}


        void addLayer(string uniqueLayerName, LayerArguments* arguments);
        void addEdge(string from, string to);

        void setLayers(map <string, LayerArguments*> layers_);
        void setEdges(vector <pair <string, string>> edges_);

        void buildNeuralNetwork(Target* target, bool serial=false);
        void buildOnlyOneLayer(Target* target, string layerName);

        static unsigned int calculateNewOutputSize(unsigned int oldSize, unsigned int padding1, unsigned int padding2, unsigned int coreSize, unsigned int stride);

        multiplicationType getMultiplicationType() const;

		static string convertMpzToString(mpz_class number, unsigned int width);

    private:
        string name;
        int fullConnectedWeightsPerDMAAccess;
        char roundingType;
        char multiplicationTypeChar;
		bool outputsAreValidAlreadyGenerated;
		bool onlyOutputClassIndex;

        map <string, LayerArguments*> layerArgs; //map: layer name -> layer arguments
        vector <pair <string, string>> edges; //pair: layer1 name, layer2 name: connection from layer1 to layer2

        LayerArguments* getLayerArguments(string layerName);

        map <string, Layer*> buildAllLayers(Target* target, Operator* globalOp);
        Layer* buildSpecificLayer(Target* target, string nameOfThisLayer);

        static vector<Layer*> needDMA(map<string, Layer*> * layers);
        void checkAndHandleDMA(Target* target, map<string, Layer*> * layers);
        vector<Layer*> sortDMALayersByPriority(vector<Layer*> layersToSort);
        void createHeaderForWeightsArray(unsigned int startaddress ,unsigned int numberOfDMALayers);
		vector<unsigned int> getConvolutionalWeightsForHeaderFile(LayerArguments* lA, ofstream &file);
		vector<unsigned int> getFullConnectedWeightsForHeaderFile(LayerArguments* lA, ofstream &file);
		unsigned int getMaskToCutOffBits(unsigned int width);

        void checkForErrors();

        void handleSerialCalculation(bool calcSerial);

        void buildAllEdges(Target* target, map<string, Layer*> *layerPtrs);
        void buildSpecificEdge(Target* target, string from, string to, Layer* fromOp, Layer* toOp);
        void checkEdges();
        void checkConnection(string from, string to);
        bool checkIfLayerExists(string layerName);

        Layer* buildConvLayer(Target* target, string convLayer);
        Layer* buildFullConnectedLayer(Target* target, string FCLayer);
        Layer* buildPoolingLayer(Target* target, string poolingLayer);

        Layer* buildInput(Target* target, string nextLayerName, Layer* nextLayer);
        Layer* buildOutput(Target* target, string lastLayerName, Layer* lastLayer);
        Operator* buildGlobalController(Target* target);

        string findLastLayerName(string actualLayerName) const;
        string findNextLayerName(string actualLayerName) const;

        pair<bool, bool> getParallelMemoryAccesses(string layerName) const;

        void instEverything(map<string, Layer*> *instMe, Operator* globalOp);
	};


}//namespace

#endif
