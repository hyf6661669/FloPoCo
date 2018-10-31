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
	class KerasToNeuralNetwork : public Operator {

    public:
        KerasToNeuralNetwork(Target* target, string pathToFolder_, bool serial=true, int wordSize_=16, int fraction_=0, int weightWordSize_=-1, int weightFraction_=-1, char roundingType=0x00, bool onlyOutputClassIndex_=false);

		// destructor
        ~KerasToNeuralNetwork() {}


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
		unsigned int layerIDCounter;
        string pathToFolder;
        vector<vector<vector<double>>> getRandomWeights(int a, int b, int c, int wordSize, int fractionalPart);

        NeuralNetwork* buildTestForZedboard(Target* target);
        NeuralNetwork* buildOneLayerForDebugging(Target* target);
        NeuralNetwork* buildFromCSVFiles(Target* target, bool serial=false, char roundingType=0x00);

        void processToplevelFile(string pathToToplevelFile);
        void readCSV(string layerName);

        void readHeaderOfCSV(string pathToFile, string layerName);
        void readWeightsOfCSV(string pathToFile, string layerName);
        void readBiasOfCSV(string pathToFile, string layerName);

        void setUndefinedInputShapes(unsigned int layerPosition);

		vector<string> orderLayerNames(vector<pair<string, string> > edges);

        mpz_class readWeightsForFullConnected(ifstream& file, LayerArguments* lA, string firstLineAfterHeader);
        mpz_class readWeightsForConvolutional(ifstream& file, LayerArguments* lA, string firstLineAfterHeader);

        string networkName;
        vector<pair<string, string> > edges;
        map<string, LayerArguments*> layerArgs;
        vector<string> orderedLayerNames;

        int wordSize;
        int fraction;
        int weightWordSize;
        int weightFraction;

        bool onlyOutputClassIndex;

		vector<double> flatten4DTensor(vector<vector<vector<vector<double>>>> t);
	};


}//namespace
