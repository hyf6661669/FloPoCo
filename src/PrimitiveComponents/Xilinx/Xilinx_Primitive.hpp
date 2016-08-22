#ifndef Xilinx_PRIMITIVE_H
#define Xilinx_PRIMITIVE_H

#include "Operator.hpp"
#include "utils.hpp"
#include <map>

namespace flopoco {

    class UniKs {
      public:
        enum Authors {
            AUTHOR_MKUMM = 0x1,
            AUTHOR_KMOELLER = 0x2,
            AUTHOR_JKAPPAUF = 0x4,
            AUTHOR_MKLEINLEIN = 0x8
        };

        static std::string getAuthorsString( const int &authors );


        // this is bad and i should feel bad
        static void addUnisimLibrary( OperatorPtr op );
    };

    // new operator class declaration
    class Xilinx_Primitive : public Operator {
        std::map<std::string, std::string> generics_;
      public:

        // constructor, defined there with two parameters (default value 0 for each)
        Xilinx_Primitive( Target *target );

        // destructor
        ~Xilinx_Primitive();

        static void checkTargetCompatibility( Target *target );

        void setGeneric( std::string name, string value );
        void setGeneric( string name, const long value );
        std::map<std::string, std::string> &getGenerics();

        std::string primitiveInstance( string instanceName );
        // Operator interface
      public:
        virtual void outputVHDL( ostream &o, string name );
        virtual void outputVHDLComponent( ostream &o, string name );

    };
}//namespace

#endif
