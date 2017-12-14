// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "NeuralNetwork.hpp"

// include all layers (including memory access)
#include "BlockRam.hpp"
#include "ConvolutionalLayer.hpp"
#include "FullConnectedLayer.hpp"
#include "ReLU.hpp"


using namespace std;
namespace flopoco {

    NeuralNetwork::NeuralNetwork(Target* target, map <string, LayerArguments> layers_, map <int, int> edges_, string name_) :
        Operator(target), layers(layers_), edges(edges_), name(name_) {

        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="NeuralNetwork";

        // author
        setCopyrightString("Nicolai Fiege, 2017");

        // definition of the name of the operator
        setName(name);
		
        // fill orderedLayers
        for(auto layerIt = layers.begin(); layerIt!=layers.end(); ++layerIt)
        {
            if((*layerIt).second.getNumber()<0)
            {
                stringstream e;
                e << "The layer '" << (*layerIt).first << "' has number<0 (" << (*layerIt).second.getNumber() << ")";
                THROWERROR(e);
            }
            orderedLayers[(*layerIt).second.getNumber()]=(*layerIt).first;
        }

        // throw errors and warnings
        if(orderedLayers.size()!=layers.size())
        {
            cout << "NeuralNetwork.constructor: Warning: Something is fishy in the argument for 'layers_'" << endl;
        }

        if(layers[orderedLayers[0]].getLayerType()!="Input")
        {
            stringstream e;
            e << "The first layer is not an Input-Layer";
            THROWERROR(e);
        }

        if(layers[orderedLayers[orderedLayers.size()-1]].getLayerType()!="Output")
        {
            stringstream e;
            e << "The last layer is not an Output-Layer";
            THROWERROR(e);
        }
		

        // create vhdl code for this network
        buildNeuralNetwork();
    }

    void NeuralNetwork::buildNeuralNetwork()
    {
        addInput("a",8);
        addOutput("b",8);

        vhdl << "b <= a;" << endl;

    }

}//namespace flopoco
