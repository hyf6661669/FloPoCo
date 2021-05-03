#pragma once

#include "../Primitive.hpp"

namespace flopoco
{

    class IntelRCCM : Operator
    {
    public:

        IntelRCCM(Operator *parentOp, Target *target, int wIn, string type);

        /** Factory method */
        static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target , vector<string> &args);
        /** Register the factory */
        static void registerFactory();

        void emulate (TestCase* tc);
        static TestList unitTest(int);

    private:
        int wIn;
        string type;

        void addLCell(int id, Target *target, string lut_mask, string dataa="'0'", string datab="'0'",
                      string datac="'0'", string datad="'0'", string datae="'0'", string dataf="'0'",
                      string cin="'0'", string sumout="open", string cout="open");

    };

} //namespace