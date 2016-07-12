// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

// include the header of the Operator
#include "GenericLut.hpp"
#include "Xilinx/Xilinx_Primitive.hpp"

using namespace std;

namespace flopoco {
    GenericLut::GenericLut( Target *target, const std::string &name, const std::map<unsigned int, unsigned int> &pairs, const unsigned int &wIn, const unsigned int &wOut ) : Operator( target ) {
        setCopyrightString( UniKs::getAuthorsString( UniKs::AUTHOR_MKLEINLEIN ) );
        Xilinx_Primitive::checkTargetCompatibility( target );
        setCombinatorial();
        stringstream tname;
        tname << "GenericLut_" << name;
        srcFileName = "GenericLut";
        setName( tname.str() );
        wOut_ = wOut;
        wIn_ = wIn;

        if( wOut_ == 0 ) {
            std::map<unsigned int, unsigned int>::const_iterator el = std::max_element( pairs.begin(), pairs.end(),
            []( const std::pair<unsigned int, unsigned int> &p1, const std::pair<unsigned int, unsigned int> &p2 ) {
                return p1.second < p2.second;
            } );
            unsigned int x = el->second;
            unsigned int c = 0;

            while( x > 0 ) {
                c++;
                x >>= 1;
            }

            wOut_ = c;
        }

        if( wIn_ == 0 ) {
            std::map<unsigned int, unsigned int>::const_iterator el = std::max_element( pairs.begin(), pairs.end(),
            []( const std::pair<unsigned int, unsigned int> &p1, const std::pair<unsigned int, unsigned int> &p2 ) {
                return p1.first < p2.first;
            } );
            unsigned int x = el->first;
            unsigned int c = 0;

            while( x > 0 ) {
                c++;
                x >>= 1;
            }

            wIn_ = c;
        }

        for( unsigned int i = 0; i < wIn_; ++i ) {
            addInput( join( "i", i ), 1 );
        }

        for( unsigned int i = 0; i < wOut_; ++i ) {
            addOutput( join( "o", i ), 1 );
        }

        for( unsigned int out = 0; out < wOut_; ++out ) {
            const unsigned int mask = ( 1 << out );
            bool_eq eq;

            for( std::map<unsigned int, unsigned int>::const_iterator it = pairs.begin();
                 it != pairs.end(); ++it ) {
                if( it->second & mask ) {
                    bool_eq part;

                    for( unsigned int in = 0; in < wIn_; ++in ) {
                        if( it->first & ( 1 << in ) ) {
                            part &= bool_eq::in( in );
                        } else {
                            part &= ~bool_eq::in( in );
                        }
                    }

                    eq |= part;
                }
            }

            equations_.push_back( eq );
        }
    }

    void GenericLut::build() {
        if( UserInterface::useTargetSpecificOptimization ) {
            //build_lut
            throw std::runtime_error( "Target specific optimizations for GenericLut not implemented yet." );
        } else {
            build_select();
        }
    }

    void GenericLut::build_select() {
        vhdl << tab << declare( "t_in", wIn_ ) << " <= ";

        for( uint i = 0; i < wIn_; ++i ) {
            if ( i > 0 ) {
                vhdl << " & ";
            }

            vhdl << join( "i", i );
        }

        vhdl << ";" << std::endl;
        declare( "t_out", wOut_ );
        vhdl << tab << "with t_in select t_out <= " << std::endl;
        const uint max_val = ( 1 << wIn_ );

        for( uint i = 0; i < max_val; ++i ) {
            std::vector<bool> input_vec( wIn_ );
            std::vector<bool> output_vec( wOut_ );

            for( uint j = 0; j < wIn_; ++j ) { // eingangswert
                input_vec[j] = i & ( 1 << j );
            }

            for( uint j = 0; j < wOut_; ++j ) { // ausgangswert
                output_vec[j] = equations_[j].eval( input_vec );
            }

            // schreiben
            vhdl << tab << tab << "\"";

            for( std::vector<bool>::reverse_iterator it = output_vec.rbegin(); it != output_vec.rend(); ++it ) {
                vhdl << ( *it ? "1" : "0" );
            }

            vhdl << "\" when \"";

            for( std::vector<bool>::reverse_iterator it = input_vec.rbegin(); it != input_vec.rend(); ++it ) {
                vhdl << ( *it ? "1" : "0" );
            }

            vhdl << "\"," << std::endl;
        }

        vhdl << tab << tab << "\"";

        for( int i = 0; i < wOut_; ++i ) {
            vhdl << "0";
        }

        vhdl << "\" when others;" << std::endl << std::endl;

        for( unsigned int i = 0; i < wIn_; ++i ) {
            vhdl << tab << join( "o", i ) << " <= " << "t_out" << of( i ) << std::endl;
        }
    }

    GenericLut::GenericLut( Target *target, const std::string &name, const std::vector<bool_eq> &equations ): Operator( target ), equations_( equations ) {
        setCopyrightString( UniKs::getAuthorsString( UniKs::AUTHOR_MKLEINLEIN ) );
        Xilinx_Primitive::checkTargetCompatibility( target );
        setCombinatorial();
        stringstream tname;
        tname << "GenericLut_" << name;
        srcFileName = "GenericLut";
        setName( tname.str() );
        wIn_ = 0;

        for( unsigned int i = 0; i < equations_.size(); ++i ) {
            const std::vector<unsigned long> inputs = equations_.at( i ).getInputs();
            std::vector<unsigned long>::const_iterator max = std::max_element( inputs.begin(), inputs.end() );

            if( *max > wIn_ ) {
                wIn_ = *max;
            }
        }

        wOut_ = equations_.size();
        build();
    }

}//namespace
