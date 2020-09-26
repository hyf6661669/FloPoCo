//
// Created by Annika Oeste on 20.06.20.
//

#include "../Operator.hpp"
#include "BitHeap.hpp"
#include "../utils.hpp"

namespace flopoco {
    class ModuloBitHeapOperator: public Operator {
    public:

            int wIn; // length of the input
            int modulo;

    public:

            ModuloBitHeapOperator(OperatorPtr parentOp, Target* target, int wIn = 0, int modulo = 1);

            ~ModuloBitHeapOperator() {};

            void emulate(TestCase * tc);

            void buildStandardTestCases(TestCaseList* tcl);

            static TestList unitTest(int index);

            /** Factory method that parses arguments and calls the constructor */
            static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target , vector<string> &args);

            /** Factory register method */
            static void registerFactory();

    protected:

        vector<vector<vector<long long> > > bits;

        int stage;

        vector<vector<long long> > rangeBH;

        void applyPseudoCompressors();

        int getMaxHeight();

        int getMSBInStage();

        int getModuloMSB();

        int reqBitsForRange(int min, int max);

        int getLeadingZero(int value);
    };
}//namespace





