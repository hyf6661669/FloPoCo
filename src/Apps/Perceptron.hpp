#include "Operator.hpp"

namespace flopoco{

    /** The FPDiv class */
    class Perceptron : public Operator
    {
    private:
        int inputLsb;
        int inputMsb;
        int expLsb;
        int sumLsb;
        int sumMsb;
        int prevLayerWidth;

    public:
        Perceptron(OperatorPtr parentOp, Target* target, int inputLsb, int inputMsb, int expLsb, int sumMsb, int sumLsb,
            int prevLayerWidth);
        ~Perceptron() {};

        void emulate(TestCase * tc);
        void buildStandardTestCases(TestCaseList* tcl);

        static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target , vector<string> &args);
        static void registerFactory();
    };
}
