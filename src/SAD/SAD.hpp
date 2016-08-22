#include "Operator.hpp"
#include "utils.hpp"

namespace flopoco {
    class SAD : public Operator {
      public:
        unsigned int wordsize;
        int total_dimension;
        enum SAD_mode {
            MODE_STD,
            MODE_KUMM,
            MODE_STRAIGHT_FORW,
            MODE_PARALLEL
        };

      public:
        SAD(Target *target, const int &wIn, const int &dimension_x, const int &dimension_y = 0, SAD_mode mode = MODE_STD, const string &addertree_type = "add2" );

        ~SAD() {};

        void emulate( TestCase *tc );
        void buildStandardTestCases( TestCaseList *tcl );

        static OperatorPtr parseArguments( Target *target, vector<string> &args );

        static void registerFactory();
	};


}//namespace
