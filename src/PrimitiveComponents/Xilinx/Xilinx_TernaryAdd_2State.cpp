// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "Xilinx_TernaryAdd_2State.hpp"
#include "Xilinx_TernaryAdd_2State_slice.hpp"
#include "PrimitiveComponents/Xilinx/Xilinx_Primitive.hpp"
#include "Xilinx_LUT_compute.h"

using namespace std;
namespace flopoco {

    Xilinx_TernaryAdd_2State::Xilinx_TernaryAdd_2State( Target *target,const int &wIn,const short &state,const short &state2 )
        : Operator( target ), wIn_(wIn), state_(state),state2_(state2) {
        setCopyrightString( UniKs::getAuthorsString( UniKs::AUTHOR_MKLEINLEIN ) );
        if( state2_ == -1 ) {
            state2_ = state_;
        }

        srcFileName = "Xilinx_TernaryAdd_2State";
        stringstream namestr;
        namestr << "Xilinx_TernaryAdd_2State_ws" << wIn << "_s" << ( state_ & 0x7 );

        if( state2_ != state_ ) {
            namestr << "_s" << ( state2_ & 0x7 );
        }

        setName( namestr.str() );
        string lut_content = computeLUT( );
        setCombinatorial();
        //addGeneric("input_word_size","integer","10");
        //addGeneric("use_carry_in","boolean","false");
        addInput( "x_i", wIn );
        addInput( "y_i", wIn );
        addInput( "z_i", wIn );

        if( state_ == state2_ ) {
            vhdl << declare( "sel_i" ) << " <= '0';" << std::endl;
        } else {
            addInput( "sel_i" );
        }

        addOutput( "sum_o", wIn );
        const uint num_slices = ( ( wIn - 1 ) / 4 ) + 1;
        declare( "x", wIn );
        declare( "y", wIn );
        declare( "z", wIn );
        declare( "bbus", wIn+1 );
        declare( "carry_cct" );
        declare( "carry", num_slices );

        if( wIn <= 0 ) {
            throw std::runtime_error( "An adder with wordsize 0 is not possible." );
        }

        insertCarryInit( );
        vhdl << tab << "x <= x_i;" << std::endl;
        vhdl << tab << "y <= y_i;" << std::endl;
        vhdl << tab << "z <= z_i;" << std::endl;

        if( num_slices == 1 ) {
            Xilinx_TernaryAdd_2State_slice *single_slice = new Xilinx_TernaryAdd_2State_slice( target, wIn, true, lut_content );
            addSubComponent( single_slice );
            inPortMap( single_slice, "x_in", "x" );
            inPortMap( single_slice, "y_in", "y" );
            inPortMap( single_slice, "z_in", "z" );
            inPortMap( single_slice, "sel_in", "sel_i" );
            inPortMap( single_slice, "bbus_in", "bbus" + range( wIn, 0 ) );
            inPortMap( single_slice, "carry_in", "carry_cct" );
            outPortMap( single_slice, "bbus_out", "bbus" + range( wIn, 0 ), false );
            outPortMap( single_slice, "carry_out", "carry" + range(0, 0 ), false );
            outPortMap( single_slice, "sum_out", "sum_o" + range( wIn, 0 ), false );
            vhdl << instance( single_slice, "slice_i" );
        }

        if ( num_slices > 1 ) {
            for( uint i = 0; i < num_slices; ++i ) {
                if( i == 0 ) {  // FIRST SLICE
                    Xilinx_TernaryAdd_2State_slice *first_slice = new Xilinx_TernaryAdd_2State_slice( target, 4, true, lut_content );
                    addSubComponent( first_slice );
                    inPortMap( first_slice, "x_in", "x" + range( 3, 0 ) );
                    inPortMap( first_slice, "y_in", "y" + range( 3, 0 ) );
                    inPortMap( first_slice, "z_in", "z" + range( 3, 0 ) );
                    inPortMap( first_slice, "sel_in", "sel_i" );
                    inPortMap( first_slice, "bbus_in", "bbus" + range( 3, 0 ) );
                    inPortMap( first_slice, "carry_in", "carry_cct" );
                    outPortMap( first_slice, "bbus_out", "bbus"  + range( 4, 1 ), false );
                    outPortMap( first_slice, "carry_out", "carry" + of( 0 ), false );
                    outPortMap( first_slice, "sum_out", "sum_o" + range( 3, 0 ), false );
                    vhdl << instance( first_slice, join( "slice_", i ) ) << endl;
                } else if( i == (num_slices - 1) ) { // LAST SLICE
                    Xilinx_TernaryAdd_2State_slice *last_slice = new Xilinx_TernaryAdd_2State_slice( target, wIn - ( 4 * i ), false, lut_content );
                    addSubComponent( last_slice );
                    inPortMap( last_slice, "x_in", "x" + range( wIn - 1, 4 * i ) );
                    inPortMap( last_slice, "y_in", "y" + range( wIn - 1, 4 * i ) );
                    inPortMap( last_slice, "z_in", "z" + range( wIn - 1, 4 * i ) );
                    inPortMap( last_slice, "sel_in", "sel_i" );
                    inPortMap( last_slice, "bbus_in", "bbus" + range( wIn - 1, 4 * i ) );
                    inPortMap( last_slice, "carry_in", "carry" + of( i - 1 ) );
                    outPortMap( last_slice, "bbus_out", "bbus" + range( wIn, 4 * i + 1 ), false );
                    outPortMap( last_slice, "carry_out", "carry" + of( i ), false );
                    outPortMap( last_slice, "sum_out", "sum_o" + range( wIn - 1, 4 * i ), false );
                    vhdl << instance( last_slice, join( "slice_", i ) ) << endl;
                } else {
                    Xilinx_TernaryAdd_2State_slice *full_slice = new Xilinx_TernaryAdd_2State_slice( target, 4, false, lut_content );
                    addSubComponent( full_slice );
                    inPortMap( full_slice, "x_in", "x" + range( ( 4 * i ) + 3, 4 * i ) );
                    inPortMap( full_slice, "y_in", "y" + range( ( 4 * i ) + 3, 4 * i ) );
                    inPortMap( full_slice, "z_in", "z" + range( ( 4 * i ) + 3, 4 * i ) );
                    inPortMap( full_slice, "sel_in", "sel_i" );
                    inPortMap( full_slice, "bbus_in", "bbus" + range( ( 4 * i ) + 3, 4 * i ) );
                    inPortMap( full_slice, "carry_in", "carry" + of( i + 1 ) );
                    outPortMap( full_slice, "bbus_out", "bbus" + range( ( 4 * i ) + 4, 4 * i + 1 ), false );
                    outPortMap( full_slice, "carry_out", "carry" + of( i ), false );
                    outPortMap( full_slice, "sum_out", "sum_o" + range( ( 4 * i ) + 3, 4 * i ), false );
                    vhdl << instance( full_slice, join( "slice_", i ) ) << endl;
                }
            }
        }
    };

    string Xilinx_TernaryAdd_2State::computeLUT() {
        lut_op add_o5_co;
        lut_op add_o6_so;

        for( int i = 0; i < 32; ++i ) {
            int s = 0;

            for( int j = 0; j < 3; j++ ) {
                if( i & ( 0x8 ) ) {
                    if( ( ( state2_ & ( 1 << ( 2 - j ) ) ) && !( i & ( 1 << j ) ) ) || ( !( state2_ & ( 1 << ( 2 - j ) ) ) && ( i & ( 1 << j ) ) ) ) {
                        s++;
                    }
                } else {
                    if( ( ( state_ & ( 1 << ( 2 - j ) ) ) && !( i & ( 1 << j ) ) ) || ( !( state_ & ( 1 << ( 2 - j ) ) ) && ( i & ( 1 << j ) ) ) ) {
                        s++;
                    }
                }
            }

            if( s & 0x2 ) {
                lut_op co_part;

                for( int j = 0; j < 5; j++ ) {
                    if ( i & ( 1 << j ) ) {
                        co_part = co_part & lut_in( j );
                    } else {
                        co_part = co_part & ( ~lut_in( j ) );
                    }
                }

                add_o5_co = add_o5_co | co_part;
            }

            if( i & 0x10 ) {
                s++;
            }

            if( s & 0x1 ) {
                lut_op so_part;

                for( int j = 0; j < 5; j++ ) {
                    if ( i & ( 1 << j ) ) {
                        so_part = so_part & lut_in( j );
                    } else {
                        so_part = so_part & ( ~lut_in( j ) );
                    }
                }

                add_o6_so = add_o6_so | so_part;
            }
        }

        lut_init lut_add3( add_o5_co, add_o6_so );
        //cerr << "Lut for {" << state << "," << state2 << "} is " << lut_add3.get_hex() << endl;
        //cerr << lut_add3.truth_table() << endl << endl;
        return lut_add3.get_hex();
    }

    void Xilinx_TernaryAdd_2State::insertCarryInit() {
        int state_ccount = 0, state2_ccount = 0;

        for( int i = 0; i < 3; i++ ) {
            if( state_ & ( 1 << i ) ) {
                state_ccount++;
            }

            if( state2_ & ( 1 << i ) ) {
                state2_ccount++;
            }
        }

        string carry_cct, bbus;

        if( abs( state_ccount - state2_ccount ) == 0 ) {
            switch( state2_ccount ) {
                case 0:
                    carry_cct = "'0'";
                    bbus = "'0'";
                    break;

                case 1:
                    carry_cct = "'0'";
                    bbus = "'1'";
                    break;

                case 2:
                    carry_cct = "'1'";
                    bbus = "'1'";
                    break;
            }
        } else if( abs( state_ccount - state2_ccount ) == 1 ) {
            if( state_ccount == 1 ) {
                if( state2_ccount == 2 ) {
                    carry_cct = "'1'";
                    bbus = "sel_i";
                } else {
                    carry_cct = "'0'";
                    bbus = "not sel_i";
                }
            } else if( state_ccount == 0 ) {
                carry_cct = "'0'";
                bbus = "sel_i";
            } else {
                carry_cct = "'1'";
                bbus = "not sel_i";
            }
        } else if( abs( state_ccount - state2_ccount ) == 2 ) {
            if( state_ccount == 0 ) {
                carry_cct = "sel_i";
                bbus = "sel_i";
            } else {
                carry_cct = "not sel_i";
                bbus = "not sel_i";
            }
        }

        if( carry_cct.empty() || bbus.empty() ) {
            throw "No carry init found";
        }

        vhdl << "\tcarry_cct <= " << carry_cct <<  ";" << endl;
        vhdl << "\tbbus(0) <= "  <<  bbus  <<  ";"  << endl;
    }

    void Xilinx_TernaryAdd_2State::computeState() {
        if( state_ == state2_ ) {
            state2_ = -1;
        }

        short states[] = {state_, state2_};
        bool lock[] = {false, false, false};
        stringstream debug;
        int k_max = 2;

        if( state2_ == -1 ) {
            k_max = 1;
        }

        mapping[0] = 0;
        mapping[1] = 1;
        mapping[2] = 2;

        if( states[0] & 0x4 ) {
            lock[2] = true;
        }

        for( int k = 0; k < k_max; k++ ) {
            for( int i = 2; i >= 0; --i ) {
                if( states[k] & ( 1 << mapping[i] ) ) {
                    for( int j = i; j < 2; ++j ) {
                        if( !lock[j + 1] ) {
                            short t      = mapping[j + 1];
                            mapping[j + 1] = mapping[j];
                            mapping[j]   = t;
                        } else {
                            break;
                        }
                    }
                }
            }

            for( int i = 2; i >= 0; --i ) {
                if( states[k] & ( 1 << mapping[i] ) ) {
                    lock[i] = true;
                }
            }
        }

        debug << "Ternary_2State with states {" << state_ << "," << state2_ << "}" << endl;
        short states_postmap[] = {0, 0};

        for( int k = 0; k < k_max; k++ ) {
            for( int i = 0; i < 3; i++ ) {
                if( states[k] & ( 1 << mapping[i] ) ) {
                    states_postmap[k] |= ( 1 << i );
                }
            }
        }

        if( state2_ == -1 ) {
            states_postmap[1] = -1;
        } else {
            if( ( states_postmap[0] & 0x6 ) == 6 && ( states_postmap[1] & 0x6 ) == 2 ) {
                short t      = mapping[1];
                mapping[1] = mapping[2];
                mapping[2]   = t;
                states_postmap[1] = ( states_postmap[1] & 0x1 ) + 4;
            }
        }

        for( int i = 0; i < 3; i++ ) {
            debug << "map " << i << " to " << mapping[i] << endl;
        }

        debug << "postmap states {" << states_postmap[0] << "," << states_postmap[1] << "}" << endl;
        state_type = 0;

        if( states_postmap[0] == 0 ) {
            if      ( states_postmap[1] == 4 ) {
                state_type = 3;
            } else if ( states_postmap[1] == 6 ) {
                state_type = 5;
            } else if ( states_postmap[1] == -1 ) {
                state_type = 0;
            }
        } else if ( states_postmap[0] == 4 ) {
            if      ( states_postmap[1] == 2 ) {
                state_type = 4;
            } else if ( states_postmap[1] == 6 ) {
                state_type = 6;
            } else if ( states_postmap[1] == 3 ) {
                state_type = 7;
            } else if ( states_postmap[1] == 0 ) {
                state_type = 3;
            } else if ( states_postmap[1] == -1 ) {
                state_type = 1;
            }
        } else if ( states_postmap[0] == 6 ) {
            if      ( states_postmap[1] == 5 ) {
                state_type = 8;
            } else if ( states_postmap[1] == 1 ) {
                state_type = 7;
            } else if ( states_postmap[1] == 0 ) {
                state_type = 5;
            } else if ( states_postmap[1] == 4 ) {
                state_type = 6;
            } else if ( states_postmap[1] == -1 ) {
                state_type = 2;
            }
        }

        debug << "found type " << state_type << endl;
        //cerr << debug.str();
        REPORT( DEBUG, debug.str() );

        if( state_ > 0 && state_type == 0 ) {
            throw "2State type not valid";
        }
    }
}//namespace
