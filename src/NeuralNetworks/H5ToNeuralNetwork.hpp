/* Each Operator declared within the flopoco framework has 
   to inherit the class Operator and overload some functions listed below*/
#include "Operator.hpp"

/* This file contains a lot of useful functions to manipulate vhdl */
#include "utils.hpp"
#include <string>

/*  All flopoco operators and utility functions are declared within
    the flopoco namespace.
    You have to use flopoco:: or using namespace flopoco in order to access these
    functions.
*/
#include "NeuralNetwork.hpp"

namespace flopoco {

	// new operator class declaration
	class H5ToNeuralNetwork : public Operator {

    public:
		H5ToNeuralNetwork(Target* target, string pathToH5);

		// destructor
        ~H5ToNeuralNetwork() {}


		// Below all the functions needed to test the operator
		/* the emulate function is used to simulate in software the operator
		   in order to compare this result with those outputed by the vhdl opertator */
        //void emulate(TestCase * tc);

		/* function used to create Standard testCase defined by the developper */
        //void buildStandardTestCases(TestCaseList* tcl);


		/* function used to bias the (uniform by default) random test generator
		   See FPExp.cpp for an example */
		// TestCase* buildRandomTestCase(int i);

		/** Factory method that parses arguments and calls the constructor */
        static OperatorPtr parseArguments(Target *target , vector<string> &args);
		
		/** Factory register method */ 
        static void registerFactory();

        NeuralNetwork* neuralNetwork; //a pointer to the actual neural network that we want to build

	private:
		string pathToH5;
        vector<vector<vector<double>>> getRandomWeights(int a, int b, int c, int wordSize, int fractionalPart);
	};


}//namespace
