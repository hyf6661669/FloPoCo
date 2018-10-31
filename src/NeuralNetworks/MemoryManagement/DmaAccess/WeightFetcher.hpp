#ifndef WEIGHTFETCHER_H
#define WEIGHTFETCHER_H
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

namespace flopoco {

    // new operator class declaration
    class WeightFetcher : public Operator {

    public:
        WeightFetcher(Target* target, unsigned int weightsPerAccess_, unsigned int weightWidth, unsigned int startAddress_, unsigned int innerCounterMax, unsigned int outerCounterMax, unsigned int weightsPerOuterCounterStep=0);

        // destructor
        ~WeightFetcher() {}

        // communication with DMA controller
        static string getDataInPortName() {return "Data_In";}
        static string getValidInPortName() {return "Valid_In";}
        static string getLastInPortName() {return "Last_In";}
        static string getNewReadAccessPortName() {return "New_Read_Access_Out";}
        static string getNumberOfBytesPortName() {return "Number_Of_Bytes_Out";}
        static string getStartAddressPortName() {return "Start_Address_Out";}

        // communication with the rest of the layer
        static string getNextWeightsValidPortName() {return "Next_Weights_Valid_Out";}
        static string getInnerCounterPortName() {return "Inner_Counter_In";}
        static string getOuterCounterPortName() {return "Outer_Counter_In";}
        static string getStartNewReadPortName() {return "Start_New_Read_In";}
        static string getWeightsArrivedAtConvCorePortName() {return "Weights_Arrived_In";}
        string getNextWeightsPortName(unsigned int weightNumber) const;

    private:
        unsigned int weightsPerAccess;
        unsigned int startAddress;
        map<unsigned int, unsigned  int> getLUTData(unsigned int outerCounterMax, unsigned int innerCounterMax, unsigned int dataBeatsPerAccess, unsigned int dataBeatsPerLastAccess) const;
    };


}//namespace

#endif