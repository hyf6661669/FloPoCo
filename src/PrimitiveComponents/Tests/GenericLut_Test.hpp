#ifndef GenericLut_Test_H
#define GenericLut_Test_H

#include "Operator.hpp"
#include "utils.hpp"

#include <vector>
#include <string>
#include "PrimitiveComponents/bool_eq.hpp"

namespace flopoco {

	// new operator class declaration
    class GenericLut_Test : public Operator {
      public:
        GenericLut_Test( Target *target );
        ~GenericLut_Test() {};

        void build_map();


        // Operator interface
      public:
        virtual void emulate( TestCase *tc );

        static OperatorPtr parseArguments( Target *target, vector<string> &args );

        static void registerFactory();
    };


}//namespace

#endif
