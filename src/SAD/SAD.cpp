#include <iostream>
#include <sstream>
#include "gmp.h"
#include "mpfr.h"


#include "SAD.hpp"
#include "PrimitiveComponents/Xilinx/Xilinx_GenericAddSub.hpp"
#include "PrimitiveComponents/Xilinx/Xilinx_Comparator.hpp"
#include "PrimitiveComponents/Xilinx/Xilinx_GenericMux.hpp"
#include "BitHeap/BitHeap.hpp"
#include "IntAddSubCmp/IntAdderTree.hpp"

using namespace std;
namespace flopoco {
    SAD::SAD( Target *target, const int &wIn, const int &dimension_x, const int &dimension_y, SAD_mode mode, const string &addertree_type )
        : Operator( target ),
          wordsize( wIn ) {
        setCopyrightString( UniKs::getAuthorsString( UniKs::AUTHOR_MKUMM | UniKs::AUTHOR_MKLEINLEIN ) );

        total_dimension = dimension_x * ( dimension_y > 0 ? dimension_y : dimension_x );
        srcFileName = "SAD";
        ostringstream name;
        name << "SAD";

        switch( mode ) {
            case flopoco::SAD::MODE_STD:
                break;

            case flopoco::SAD::MODE_KUMM:
                name << "_KUMM";
                break;

            case flopoco::SAD::MODE_STRAIGHT_FORW:
                name << "_STF";
                break;

            case flopoco::SAD::MODE_PARALLEL:
                name << "_PAR";
                break;
        }

        name << "_" << dimension_x << "x" << ( dimension_y > 0 ? dimension_y : dimension_x ) << "_" << wIn;
        setName( name.str() );
        useNumericStd();

        for( int i = 0; i < total_dimension; i++ ) {
            addInput( join( "iA", "_", i ), wIn );
            addInput( join( "iB", "_", i ), wIn );

            if( mode == MODE_KUMM || mode == MODE_PARALLEL || mode == MODE_STRAIGHT_FORW ) {
                vhdl << declare( join( "iAe", "_", i ), wIn + 1 ) << " <= '0' & " << join( "iA", "_", i ) << ";" << endl;
                vhdl << declare( join( "iBe", "_", i ), wIn + 1 ) << " <= '0' & " << join( "iB", "_", i ) << ";" << endl;
            }
        }

        int pack4_count = total_dimension / 2;
        //int wOut = wIn + ( sizeof( int ) * 8 - __builtin_clz( total_dimension ) ) - 1;
        Xilinx_GenericAddSub *addsub;

        if( mode == MODE_KUMM ) {
            int treedim = 0;

            for( int i = 0; i < pack4_count; i++ ) {
                if( UserInterface::useTargetSpecificOptimization ) {
                    addsub = new Xilinx_GenericAddSub( target, wIn + 1, 2 );
                    addSubComponent( addsub );
                    inPortMap( addsub, "x_i", join( "iAe", "_", i * 2 ) );
                    inPortMap( addsub, "y_i", join( "iBe", "_", i * 2 ) );
                    outPortMap( addsub, "sum_o", join( "oL", "_", i ) );
                    vhdl << instance( addsub, join( "subL_", i ) );
                    addsub = new Xilinx_GenericAddSub( target, wIn + 1, 2 );
                    addSubComponent( addsub );
                    inPortMap( addsub, "x_i", join( "iAe", "_", i * 2 + 1 ) );
                    inPortMap( addsub, "y_i", join( "iBe", "_", i * 2 + 1 ) );
                    outPortMap( addsub, "sum_o", join( "oR", "_", i ) );
                    vhdl << instance( addsub, join( "subR_", i ) );
                    addsub = new Xilinx_GenericAddSub( target, wIn + 1, true );
                    addSubComponent( addsub );
                    inPortMap( addsub, "x_i", join( "oL", "_", i ) );
                    inPortMap( addsub, "y_i", join( "oR", "_", i ) );
                    inPortMapCst( addsub, "neg_x_i", join( "oL_", i ) + of( wIn ) );
                    inPortMapCst( addsub, "neg_y_i", join( "oR_", i ) + of( wIn ) );
                    outPortMap( addsub, "sum_o", join( "o", "_", i ) );
                    vhdl << instance( addsub, join( "sub_", i ) );
                } else {
                    vhdl << declare( join( "oL", "_", i ), wIn + 1 ) << " <= std_logic_vector(signed(" << join( "iAe", "_", i * 2 ) << ") - signed(" <<  join( "iBe", "_", i * 2 ) << "));" << std::endl;
                    vhdl << declare( join( "oR", "_", i ), wIn + 1 ) << " <= std_logic_vector(signed(" << join( "iAe", "_", i * 2 + 1 ) << ") - signed(" <<  join( "iBe", "_", i * 2 + 1 ) << "));" << std::endl;
                    vhdl << declare( join( "o", "_", i ), wIn + 1 ) << " <= std_logic_vector(abs(signed(" << join( "oL", "_", i ) << ")) + abs(signed(" << join( "oR", "_", i ) << ")));" << std::endl;
                }

                treedim++;
            }

            if( ( total_dimension % 2 ) == 1 ) {
                if( UserInterface::useTargetSpecificOptimization ) {
                    Xilinx_Comparator *comp = new Xilinx_Comparator( target, wIn + 1, Xilinx_Comparator::ComparatorType_gt );
                    addSubComponent( comp );
                    inPortMap( comp, "a", join( "iAe", "_", total_dimension - 1 ) );
                    inPortMap( comp, "b", join( "iBe", "_", total_dimension - 1 ) );
                    outPortMap( comp, "o", join( "comp", "_", pack4_count ) );
                    vhdl << instance( comp, "compRem" );
                    vhdl << declare( join( "comp", "_", pack4_count, "_inv" ) ) << " <= not " << join( "comp", "_", pack4_count ) << ";";
                    addsub = new Xilinx_GenericAddSub( target, wIn + 1, false );
                    addSubComponent( addsub );
                    inPortMap( addsub, "x_i", join( "iAe", "_", total_dimension - 1 ) );
                    inPortMap( addsub, "y_i", join( "iBe", "_", total_dimension - 1 ) );
                    inPortMapCst( addsub, "neg_x_i", join( "comp", "_", pack4_count, "_inv" ) );
                    inPortMapCst( addsub, "neg_y_i", join( "comp", "_", pack4_count ) );
                    outPortMap( addsub, "sum_o", join( "oRem", "_", pack4_count ) );
                    vhdl << instance( addsub, join( "subRem_", pack4_count ) );
                } else {
                    vhdl << declare( join( "oRem", "_", pack4_count ), wIn + 1 ) << " <= std_logic_vector(abs(signed(" << join( "iAe", "_", total_dimension - 1 ) << ") - signed(" <<  join( "iBe", "_", total_dimension - 1 ) << ")));" << std::endl;
                }

                treedim++;
            }

            IntAdderTree *atree =  new IntAdderTree( target, wIn + 1 , treedim , addertree_type );
            addSubComponent( atree );

            for( int i = 0; i < pack4_count; i++ ) {
                inPortMap( atree, join( "X", i + 1 ), join( "o", "_", i ) );
            }

            if( ( total_dimension % 2 ) == 1 ) {
                inPortMap( atree, join( "X", pack4_count + 1 ), join( "oRem", "_", pack4_count ) );
            }

            outPortMap( atree, "Y", "adder_tree_sum" );
            addOutput( "S" ,  getSignalByName( "adder_tree_sum" )->width() );
            vhdl << instance( atree, "addertree" ) << std::endl;
            setCycleFromSignal( "adder_tree_sum" );
            vhdl << "S <= adder_tree_sum;" << std::endl;
        } else if( mode == MODE_STRAIGHT_FORW ) {
            for( int i = 0; i < total_dimension; ++i ) {
                if( UserInterface::useTargetSpecificOptimization ) {
                    addsub = new Xilinx_GenericAddSub( target, wIn + 1, 2 ); //(a-b)
                    addSubComponent( addsub );
                    inPortMap( addsub, "x_i", join( "iAe", "_", i ) );
                    inPortMap( addsub, "y_i", join( "iBe", "_", i ) );
                    outPortMap( addsub, "sum_o", join( "oa", "_", i ) );
                    vhdl << instance( addsub, join( "suba_", i ) );
                } else {
                    vhdl << declare( join( "oa", "_", i ), wIn + 1 ) << " <= std_logic_vector(signed(" << join( "iAe", "_", i ) << ") - signed(" << join( "iBe", "_", i ) << "));" << std::endl;
                }

                vhdl << declare( join( "o", "_", i ), wIn + 1 ) << " <= std_logic_vector(abs(signed(" << join( "oa", "_", i ) << ")));" << std::endl;
            }

            IntAdderTree *atree =  new IntAdderTree( target, wIn + 1 , total_dimension , addertree_type );
            addSubComponent( atree );

            for( int i = 0; i < total_dimension; ++i ) {
                inPortMap( atree, join( "X", i + 1 ), join( "o", "_", i ) );
            }

            outPortMap( atree, "Y", "adder_tree_sum" );
            addOutput( "S" ,  getSignalByName( "adder_tree_sum" )->width() );
            vhdl << instance( atree, "addertree" ) << std::endl;
            setCycleFromSignal( "adder_tree_sum" );
            vhdl << "S <= adder_tree_sum;" << std::endl;
        } else if( mode == MODE_PARALLEL ) {
            for( int i = 0; i < total_dimension; ++i ) {
                if( UserInterface::useTargetSpecificOptimization ) {
                    addsub = new Xilinx_GenericAddSub( target, wIn + 1, 2 ); //(a-b)
                    addSubComponent( addsub );
                    inPortMap( addsub, "x_i", join( "iAe", "_", i ) );
                    inPortMap( addsub, "y_i", join( "iBe", "_", i ) );
                    outPortMap( addsub, "sum_o", join( "oa", "_", i ) );
                    vhdl << instance( addsub, join( "suba_", i ) );
                    addsub = new Xilinx_GenericAddSub( target, wIn + 1, 1 ); //(a-b)
                    addSubComponent( addsub );
                    inPortMap( addsub, "x_i", join( "iAe", "_", i ) );
                    inPortMap( addsub, "y_i", join( "iBe", "_", i ) );
                    outPortMap( addsub, "sum_o", join( "ob", "_", i ) );
                    vhdl << instance( addsub, join( "subb_", i ) );
                    Xilinx_GenericMux *mux = new Xilinx_GenericMux( target, 2, wIn + 1 );
                    addSubComponent( mux );
                    inPortMap( mux, "x0_in", join( "oa_", i ) );
                    inPortMap( mux, "x1_in", join( "ob_", i ) );
                    inPortMap( mux, "s_in", join( "oa_", i ) + of( wIn ) );
                    outPortMap( mux, "x_out", join( "o", "_", i ) );
                    vhdl << instance( mux, join( "mux_", i ) );
                } else {
                    vhdl << declare( join( "oa", "_", i ), wIn + 1 ) << " <= std_logic_vector(signed(" << join( "iAe", "_", i ) << ") - signed(" <<  join( "iBe", "_", i ) << "));" << std::endl;
                    vhdl << declare( join( "ob", "_", i ), wIn + 1 ) << " <= std_logic_vector(signed(" << join( "iBe", "_", i ) << ") - signed(" <<  join( "iAe", "_", i ) << "));" << std::endl;
                    vhdl << declare( join( "o", "_", i ), wIn + 1 ) << " <= " << join( "oa", "_", i ) << " when (" << join( "oa", "_", i ) << of( wIn ) << " = '0') else " << join( "ob", "_", i ) << ";" << std::endl;
                }
            }

            IntAdderTree *atree =  new IntAdderTree( target, wIn + 1 , total_dimension , addertree_type );
            addSubComponent( atree );

            for( int i = 0; i < total_dimension; ++i ) {
                inPortMap( atree, join( "X", i + 1 ), join( "o", "_", i ) );
            }

            outPortMap( atree, "Y", "adder_tree_sum" );
            addOutput( "S" ,  getSignalByName( "adder_tree_sum" )->width() );
            vhdl << instance( atree, "addertree" ) << std::endl;
            setCycleFromSignal( "adder_tree_sum" );
            vhdl << "S <= adder_tree_sum;" << std::endl;
        } else {
            for( int i = 0; i < total_dimension; ++i ) {
                if( UserInterface::useTargetSpecificOptimization ) {
                    Xilinx_Comparator *comp = new Xilinx_Comparator( target, wIn, Xilinx_Comparator::ComparatorType_gt );
                    addSubComponent( comp );
                    inPortMap( comp, "a", join( "iA", "_", i ) );
                    inPortMap( comp, "b", join( "iB", "_", i ) );
                    outPortMap( comp, "o", join( "comp", "_", i ) );
                    vhdl << instance( comp, join( "compi_", i ) );
                    vhdl << declare( join( "comp", "_", i, "_inv" ) ) << " <= not " << join( "comp", "_", i ) << ";";
                    addsub = new Xilinx_GenericAddSub( target, wIn, false );
                    addSubComponent( addsub );
                    inPortMap( addsub, "x_i", join( "iA", "_", i ) );
                    inPortMap( addsub, "y_i", join( "iB", "_", i ) );
                    inPortMapCst( addsub, "neg_x_i", join( "comp", "_", i, "_inv" ) );
                    inPortMapCst( addsub, "neg_y_i", join( "comp", "_", i ) );
                    outPortMap( addsub, "sum_o", join( "o", "_", i ) );
                    vhdl << instance( addsub, join( "sub_", i ) );
                } else {
                    vhdl << declare( join( "o", "_", i ), wIn ) << " <= std_logic_vector(unsigned(" << join( "iA", "_", i ) << ") - unsigned(" << join( "iB", "_", i ) << ")) when (unsigned(" << join( "iA", "_", i ) << ") > unsigned(" << join( "iB", "_", i ) << ")) else " << std::endl;
                    vhdl << tab << tab << "std_logic_vector(unsigned(" << join( "iB", "_", i ) << ") - unsigned(" << join( "iA", "_", i ) << "));" << std::endl;
                }
            }

            IntAdderTree *atree =  new IntAdderTree( target, wIn , total_dimension , addertree_type );
            addSubComponent( atree );

            for( int i = 0; i < total_dimension; ++i ) {
                inPortMap( atree, join( "X", i + 1 ), join( "o", "_", i ) );
            }

            outPortMap( atree, "Y", "adder_tree_sum" );
            addOutput( "S" ,  getSignalByName( "adder_tree_sum" )->width() );
            vhdl << instance( atree, "addertree" ) << std::endl;
            setCycleFromSignal( "adder_tree_sum" );
            vhdl << "S <= adder_tree_sum;" << std::endl;
        }
    }


    void SAD::emulate( TestCase *tc ) {
        //cerr << "Testcase" << endl;
        mpz_class expected_output = 0;

        for( int i = 0; i < total_dimension; ++i ) {
            mpz_class a = tc->getInputValue( join( "iA", "_", i ) );
            mpz_class b = tc->getInputValue( join( "iB", "_", i ) );

            //cerr << "(" << a << "," << b << ",";
            if( a > b ) {
                expected_output += ( a - b );
                //cerr << (a - b) << "," << expected_output << ")" << endl;
            } else {
                expected_output += ( b - a );
                //cerr << (b - a) << "," << expected_output << ")" << endl;
            }
        }

        tc->addExpectedOutput( "S", expected_output );
	}


    void SAD::buildStandardTestCases( TestCaseList *tcl ) {
		// please fill me with regression tests or corner case tests!
    }

    OperatorPtr SAD::parseArguments( Target *target, vector<string> &args ) {
        int wIn, x = 0, y = 0;
        std::string typestr;
        std::string addertree_type = "add2";
        UserInterface::parseInt( args, "wIn", &wIn );
        UserInterface::parseString( args, "type", &typestr );
        UserInterface::parseInt( args, "dimX", &x );

        if( args.size() > 0 ) {
            UserInterface::parseInt( args, "dimY", &y );
        }

        if( args.size() > 0 ) {
            UserInterface::parseString( args, "treetype", &addertree_type );
        }

        SAD::SAD_mode mode;

        if( typestr == "STD" ) {
            mode = MODE_STD;
        } else if( typestr == "KUMM" ) {
            mode = MODE_KUMM;
        } else if( typestr == "PAR" ) {
            mode = MODE_PARALLEL;
        } else if( typestr == "STF" ) {
            mode = MODE_STRAIGHT_FORW;
        } else {
            std::cerr << "No valid type given." << std::endl;
            exit( 0 );
        }

        return new SAD( target, wIn, x, y, mode, addertree_type );
    }

    void SAD::registerFactory() {
        UserInterface::add( "SAD", // name
                            "A component for (S)um of (A)bsolute (D)ifferences computation. Please use target specific optimizations for best results (Xilinx only).", // description, string
                            "BasicInteger", // category, from the list defined in UserInterface.cpp
                            "Proposed at FPL 2016", //seeAlso
                            // where parameterDescription is parameterName (parameterType)[=defaultValue]: parameterDescriptionString
                            "wIn (int): Wordsize of comparator inputs; \
                            type (string): Type of the precompution stage.\
                            STD - build using comparator adder pairs, \
                            KUMM - build using adders which support inversion on both inputs, \
                            PAR - compute both possible and decide using a mux, \
                            STF - compute one possible and take absolute; \
                            dimX (int): Dimension of the SAD in x; \
                            dimY (int)=0: Dimension of the SAD in y, If not given is set to dimX; \
                            treetype (string)=add2: Type of the adders in the addertree used after precompution",
                            "Nope.",
                            SAD::parseArguments
                          ) ;
    }
}//namespace
