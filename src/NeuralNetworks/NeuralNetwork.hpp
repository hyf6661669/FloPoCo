#ifndef NEURALNETWORK_H
#define NEURALNETWORK_H
/* Each Operator declared within the flopoco framework has 
   to inherit the class Operator and overload some functions listed below*/
#include "Operator.hpp"

/* This file contains a lot of useful functions to manipulate vhdl */
#include "utils.hpp"
#include <string>
#include <vector>

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
        NeuralNetwork(Target* target, map <string, LayerArguments> layers_, map <int, int> edges_, string name_ = "NeuralNetwork");

		// destructor
        ~NeuralNetwork() {}
		
		map <string, LayerArguments> layers; //map: layer name -> layer arguments
        map <int, string> orderedLayers; //map: layer number -> layer name
        map <int, int> edges; //map: layer1 number -> layer2 number: connection from layer1 to layer2
		
		string name;

        void buildNeuralNetwork();

	private:
	};


}//namespace

#endif
