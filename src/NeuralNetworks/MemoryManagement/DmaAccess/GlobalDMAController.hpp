#ifndef GLOBALDMACONTROLLER_H
#define GLOBALDMACONTROLLER_H
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
    class GlobalDMAController : public Operator {

    public:
        GlobalDMAController(Target* target, unsigned int numberOfInputs_);

        // destructor
        ~GlobalDMAController() {}

        // signals to communicate with layers
        string getNewReadAccessPortName(unsigned int portNumber) const;
        string getDataPortName(unsigned int portNumber) const;
        string getNewStartAddressPortName(unsigned int portNumber) const;
        string getLastDataPortName(unsigned int portNumber) const;
        string getNumberOfTotalBytesPortName(unsigned int portNumber) const;

        // signals to read data stream from DMA (AXIS)
        static string getDMAReadyPortName() { return "ReadyToDMA"; }
        static string getDMAValidPortName() { return "ValidFromDMA"; }
        static string getDMADataPortName() { return "DataFromDMA"; }
        static string getDMALastPortName() { return "LastFromDMA"; }

        // signals for programming DMA (AXI4-Lite)
        static string getNewDMAAccessPortName() { return "NewDMAAccess"; }
        static string getDMAStartAddressPortName() { return "DMAStartAddress"; }
        static string getDMANumberOfBytesPortName() { return "DMANoOfBytes"; }
        static string getDMAFinishedProgrammingPortName() { return "DMAFinishedProgramming"; }
        static string getDMAReadyToBeProgrammedPortName() { return "DMAReadyToBeProgrammed"; }

    private:
        unsigned int numberOfInputs;
    };


}//namespace

#endif