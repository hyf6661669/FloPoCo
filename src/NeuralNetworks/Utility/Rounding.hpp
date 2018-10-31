#ifndef ROUNDING_H
#define ROUNDING_H
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
    typedef enum {
        truncation,
        saturation
    } roundingTypeEnum;

    // new operator class declaration
    class Rounding : public Operator {

    public:
        Rounding(Target* target, unsigned int wordSizeFrom, unsigned int fractionFrom, unsigned int wordSizeTo, unsigned int fractionTo, roundingTypeEnum round);

        // destructor
        ~Rounding() {}

    private:
        void buildTruncation(Target* target, unsigned int wordSizeFrom, unsigned int fractionFrom, unsigned int wordSizeTo, unsigned int fractionTo);
        void buildSaturation(Target* target, unsigned int wordSizeFrom, unsigned int fractionFrom, unsigned int wordSizeTo, unsigned int fractionTo);
    };


}//namespace

#endif