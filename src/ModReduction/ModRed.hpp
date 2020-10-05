#ifndef FLOPOCO_MODRED_H
#define FLOPOCO_MODRED_H
#include <vector>
#include <sstream>
#include <gmp.h>
#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"
#include "Table.hpp"
#include "BitHeap/BitHeap.hpp"

namespace flopoco {
    class ModRed : public Operator {

    public:

        /**
         * The ModRed constructor
         * @param[in] target           the target device
         * @param[in] wIn              input width of vector the modulo is calculated for
         * @param[in] mod              constant modulus
         **/
        ModRed(Operator *parentOp, Target *target, int wIn, int modulus);

        /**
         * The emulate function.
         * @param[in] tc               a test-case
         */
        void emulate(TestCase *tc);

        static TestList unitTest(int index);

        void buildStandardTestCases(TestCaseList *tcl);

        /** Factory method that parses arguments and calls the constructor */
        static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args);

        /** Factory register method */
        static void registerFactory();

    private:
        int wIn, mod;

    };
}



#endif //FLOPOCO_MODRED_H
