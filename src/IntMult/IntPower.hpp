#include "Operator.hpp"
#include "utils.hpp"
#include "GenericBinaryPolynomial.hpp"

namespace flopoco {

class IntPower : public GenericBinaryPolynomial {
  public:
    IntPower(Target* target, size_t wIn, size_t n, std::map<string,double> inputDelays = emptyDelayMap);

    ~IntPower() {};

    void emulate(TestCase * tc);
    void buildStandardTestCases(TestCaseList* tcl);
    void buildRandomTestCases(TestCaseList* tcl, int n);
    TestCase* buildRandomTestCases(int i);

    static OperatorPtr parseArguments( Target *target, vector<string> &args );
    static void registerFactory();

protected:
    size_t wIn, n;
};

}

