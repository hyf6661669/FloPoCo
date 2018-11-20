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
            addInput( join( "i", i ), 1, false );
        }

        for( unsigned int i = 0; i < wOut_; ++i ) {
            addOutput( join( "o", i ), 1, 1, false  );
        }

        ////////Build Lut////////
        if( UserInterface::useTargetSpecificOptimization ) {
            //build_lut
            throw std::runtime_error( "Target specific optimizations for GenericLut not implemented yet." );
        } else {

            declare( "t_in", wIn_ );

            for( uint i = 0; i < wIn_; ++i ) {

                vhdl << tab << "t_in(" << i << ") <= " << join( "i", i ) << ";" << endl;
            }

            declare( "t_out", wOut_ );
            vhdl << tab << "with t_in select t_out <= " << std::endl;

            for(auto mapIt=pairs.begin(); mapIt!=pairs.end(); ++mapIt)
            {
                unsigned int inputNumber=mapIt->first;
                unsigned int outputNumber=mapIt->second;

                vhdl << tab << tab << "\"";

                //build binary output out of output number
                for(unsigned int i=0; i<this->wOut_; i++)
                {
                    vhdl << ((outputNumber>>(this->wOut_-1-i))&1);
                }

                vhdl << "\" when \"";
                //build binary input out of input number
                for(unsigned int i=0; i<this->wIn_; i++)
                {
                    vhdl << ((inputNumber>>(this->wIn_-1-i))&1);
                }
                vhdl << "\"," << endl;

            }
            //set default case
            vhdl << tab << tab << "\"";
            for(unsigned int i=0; i<wOut_; i++)
            {
                vhdl << "0";
            }
            vhdl << "\" when others;" << endl << endl;

            //handle critical path
            this->manageCriticalPath(target->localWireDelay() + target->lutDelay());

            //build output signals
            for( unsigned int i = 0; i < wOut_; ++i ) {
                vhdl << tab << join( "o", i ) << " <= " << "t_out" << of( i ) << ";" << std::endl;
            }
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
        declare( "t_in", wIn_ );

        for( uint i = 0; i < wIn_; ++i ) {

            vhdl << "t_in(" << i << ") <= " << join( "i", i ) << ";" << endl;
        }

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

            bool some_set = false;

            for( uint i = 0; i < wOut_; ++i ) {
                if( output_vec[i] ) {
                    some_set = true;
                    break;
                }
            }

            // write
            if( some_set ) {
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
        }

        vhdl << tab << tab << "\"";

        for( uint i = 0; i < wOut_; ++i ) {
            vhdl << "0";
        }

        vhdl << "\" when others;" << std::endl << std::endl;

        for( unsigned int i = 0; i < wOut_; ++i ) {
            vhdl << tab << join( "o", i ) << " <= " << "t_out" << of( i ) << ";" << std::endl;
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
