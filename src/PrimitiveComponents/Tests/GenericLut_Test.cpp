// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

// include the header of the Operator
#include "GenericLut_Test.hpp"
#include "PrimitiveComponents/GenericLut.hpp"
#include "PrimitiveComponents/Xilinx/Xilinx_Primitive.hpp"

using namespace std;

flopoco::GenericLut_Test::GenericLut_Test( Target *target ) : Operator( target ) {
    setCopyrightString( UniKs::getAuthorsString( UniKs::AUTHOR_MKLEINLEIN ) );
    setCombinatorial();
    stringstream tname;
    tname << "GenericLut_Test";
    srcFileName = "GenericLut_Test";
    setName( tname.str() );
    build_map();
}


void flopoco::GenericLut_Test::build_map() {
    const uint max_pairs = 20;
    const uint wIn = 6;
    const uint wOut = 2;
    std::map<uint, uint> values;

    for ( uint i = 0; i < max_pairs; ++i ) {
        const uint rand1 = rand() % ( 1 << wIn );
        const uint rand2 = rand() % ( (1 << wOut) -1 ) + 1;
        values[ rand1 ] = rand2;
    }

    addInput(  "i", wIn );
    addOutput( "o", wOut );
    GenericLut *lut = new GenericLut( this->getTarget(), "testlut", values );
    addSubComponent( lut );

    for( uint i = 0; i < wIn; ++i ) {
        inPortMap( lut, join( "i", i ), "i" + of( i ) );
    }

    for( uint i = 0; i < wOut; ++i ) {
        outPortMap( lut, join( "o", i ), "o" + of( i ),false );
    }

    vhdl << instance( lut, "some_lut" );
}

void flopoco::GenericLut_Test::emulate( flopoco::TestCase *tc ) {
}

flopoco::OperatorPtr flopoco::GenericLut_Test::parseArguments( flopoco::Target *target, vector<string> &args ) {
    return new GenericLut_Test( target );
}

void flopoco::GenericLut_Test::registerFactory() {
    UserInterface::add( "genericlut_test", // name
                        "A test operator for generic lut.", // description, string
                        "Misc", // category, from the list defined in UserInterface.cpp
                        "",
                        "",
                        "Nope.",
                        GenericLut_Test::parseArguments
                      );
}
