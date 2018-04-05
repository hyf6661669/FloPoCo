PROJECT(FloPoCo C CXX)

MESSAGE("Trying to build FloPoCo on " ${CMAKE_SYSTEM_NAME} ", compiler is " ${CMAKE_C_COMPILER})

# CMake's way of handling various versions
cmake_minimum_required(VERSION 2.4)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()
#####################For packaging (source file distrib so far)#######
# Use  make package_source to build a release
# Currently it doesn't work at all

SET(CPACK_SOURCE_GENERATOR "TGZ")
SET(CPACK_PACKAGE_VERSION_MAJOR "0")
SET(CPACK_PACKAGE_VERSION_MINOR "9")
SET(CPACK_PACKAGE_VERSION_PATCH "0")

#SET(Boost_INCLUDE_DIR /usr/include/boost)
#SET(Boost_LIBRARY_DIR /usr/lib64/boost)
find_package( Boost 1.48 REQUIRED )

INCLUDE(CPack)
# 
#Compilation flags
#SET(CMAKE_CXX_FLAGS_DEBUG "-Wall")
# Uncomment for profiling
SET(CMAKE_CXX_FLAGS_DEBUG "-g -pg -Wall")
#SET(CMAKE_CXX_FLAGS_DEBUG "-pg -Wall")
# Do you want an executable of 30MB ?
SET(CMAKE_BUILD_TYPE DEBUG)
#SET(CMAKE_BUILD_TYPE release)

FIND_PATH(MPFR_H
	  mpfr.h
	  HINTS "C:/Program Files/Microsoft Visual Studio 9.0/VC/include"
	  DOC "Path of mpfr.h, the include file for GNU MPFR library"
)

FIND_LIBRARY(MPFR_LIB
	NAMES mpfr mpfr.lib
	HINTS "C:/Program Files/Microsoft Visual Studio 9.0/VC/lib"
	DOC "Directory of the MPFR library")

FIND_PATH(MPFI_H
	  mpfi.h
	  HINTS "C:/Program Files/Microsoft Visual Studio 9.0/VC/include"
	  DOC "Path of mpfi.h"
)

FIND_LIBRARY(MPFI_LIB
	NAMES mpfi mpfi.lib
	HINTS "C:/Program Files/Microsoft Visual Studio 9.0/VC/lib"
	DOC "Directory of the MPFI library")

# Detect Python
FIND_PATH(PYTHON_H
	  python2.7/Python.h
	  DOC "Path of Python.h"
)

FIND_LIBRARY(PYTHON_LIB
	NAMES libpython2.7.so libpython2.7.dylib
	HINTS "/usr/lib/x86_64-linux-gnu"
	DOC "Directory of the Python library"
)

# Detect Sollya
FIND_LIBRARY( SOLLYA_LIB
	NAMES libsollya.so libsollya.dylib
	HINTS "/usr/lib/x86_64-linux-gnu"
	)

FIND_PATH(SOLLYA_H
	  sollya.h
)

IF (SOLLYA_H AND SOLLYA_H)
	MESSAGE(STATUS "Found Sollya: ${SOLLYA_LIB},${SOLLYA_H}")
ELSE (SOLLYA_H AND SOLLYA_H)
	MESSAGE(FATAL_ERROR "Sollya not found (found: ${SOLLYA_LIB},${SOLLYA_H})")
ENDIF (SOLLYA_H AND SOLLYA_H)


IF ( PYTHON_H AND PYTHON_LIB )
    SET(PYTHON_FOUND TRUE)
ENDIF ( PYTHON_H AND PYTHON_LIB )

IF ( PYTHON_FOUND )
	MESSAGE ( STATUS "Found Python : ${PYTHON_H},${PYTHON_LIB}" )
ELSE ( PYTHON_FOUND )
	MESSAGE ( STATUS "Could not find (found: ${PYTHON_H},${PYTHON_LIB}), no flopoco_soaktest will be available")
ENDIF ( PYTHON_FOUND )

IF (MPFR_H AND MPFR_LIB)
    SET(MPFR_FOUND TRUE)
ENDIF (MPFR_H AND MPFR_LIB)

IF (MPFR_FOUND)
    MESSAGE(STATUS "Found mpfr: ${MPFR_H} ${MPFR_LIB}")
ELSE (MPFR_FOUND)
    MESSAGE(FATAL_ERROR "Could not find MPFR.\nOn Debian-based systems type: 'sudo apt-get install libmpfr-dev'")
ENDIF (MPFR_FOUND)


FIND_PATH(WCPG_H
	  wcpg.h
	  DOC "Path of wcpg.h, the include file for GNU MPFR library"
		)

FIND_LIBRARY(WCPG_LIB
	NAMES wcpg
	DOC "Directory of the WCPG library"
	)

IF (WCPG_H AND WCPG_LIB)
    SET(WCPG_FOUND TRUE)
ENDIF (WCPG_H AND WCPG_LIB)

IF ( WCPG_FOUND )
	MESSAGE ( STATUS "Found WCPG : ${WCPG_H} ${WCPG_LIB}" )
  ADD_DEFINITIONS(-DHAVE_WCPG)
ELSE ( WCPG_FOUND )
	MESSAGE ( STATUS "Could not find WCPG, this will affect the FixIIR operator" )
ENDIF ( WCPG_FOUND )

# Detect presence of SCIP
FIND_PATH(SCIP_H
	  scip.h
	  HINTS "/usr/local/include/scip"
	  DOC "Path of scip.h, the include file for the SCIP library"
)
FIND_LIBRARY(SCIP_LIB
	NAMES scip
	HINTS "/usr/local/lib64/"
	DOC "Directory of the SCIP library")


FIND_PATH(SCALP_H
          ScaLP/Solver.h
	  HINTS "/opt/scalp/include"
	  DOC "Path of main header of scalp"
)
FIND_LIBRARY(SCALP_LIB
        NAMES libScaLP.so
	HINTS "/usr/local/lib64/"
	DOC "Directory of the SCALP main library")

FIND_LIBRARY(SCALP_CPLEX_LIB
        NAMES libScaLP-CPLEX.so
	HINTS "/usr/local/lib64/"
	DOC "Directory of the SCALP library for cplex solver backend")

FIND_LIBRARY(SCALP_GUROBI_LIB
        NAMES libScaLP-Gurobi.so
	HINTS "/usr/local/lib64/"
	DOC "Directory of the SCALP library for Gurobi solver backend")

FIND_LIBRARY(SCALP_SCIP_LIB
        NAMES libScaLP-SCIP.so
	HINTS "/usr/local/lib64/"
	DOC "Directory of the SCALP library for SCIP solver backend")

FIND_LIBRARY(SCALP_LPSOLVE_LIB
        NAMES libScaLP-LPSOLVE.so
	HINTS "/usr/local/lib64/"
	DOC "Directory of the SCALP library for LPSOLVE solver backend")

#Detect presence of supporting libraries for SCIP
FIND_LIBRARY(READLINE_LIB NAMES readline)
FIND_LIBRARY(Z_LIB NAMES z)
FIND_LIBRARY(SOPLEX_LIB NAMES soplex)
FIND_LIBRARY(NLPI_LIB NAMES nlpi.cppad)
FIND_LIBRARY(LPISPX_LIB NAMES lpispx)
FIND_LIBRARY(ZIMPL_LIB NAMES zimpl.linux.x86_64.gnu.opt)

IF (SCIP_H AND SCIP_LIB AND READLINE_LIB AND Z_LIB AND NLPI_LIB AND LPISPX_LIB AND SOPLEX_LIB AND ZIMPL_LIB)
    SET(SCIP_FOUND TRUE)
ENDIF (SCIP_H AND SCIP_LIB AND READLINE_LIB AND Z_LIB AND NLPI_LIB AND LPISPX_LIB AND SOPLEX_LIB AND ZIMPL_LIB)

IF (SCALP_H AND SCALP_LIB AND (SCALP_CPLEX_LIB OR SCALP_GUROBI_LIB OR SCALP_SCIP_LIB OR SCALP_LPSOLVE_LIB))
    SET(SCALP_FOUND TRUE)
ENDIF (SCALP_H AND SCALP_LIB AND (SCALP_CPLEX_LIB OR SCALP_GUROBI_LIB OR SCALP_SCIP_LIB OR SCALP_LPSOLVE_LIB))


IF (READLINE_LIB)
	MESSAGE(STATUS "Found readline (supporting library of SCIP): ${READLINE_LIB}")
ELSE (READLINE_LIB)
	MESSAGE(STATUS "readline (supporting library of SCIP) not found, no optimal bit heap compression available")
ENDIF (READLINE_LIB)

IF (Z_LIB)
	MESSAGE(STATUS "Found z (supporting library of SCIP): ${Z_LIB}")
ELSE (Z_LIB)
	MESSAGE(STATUS "z (supporting library of SCIP) not found, no optimal bit heap compression available")
ENDIF (Z_LIB)


IF (SOPLEX_LIB)
	MESSAGE(STATUS "Found soplex (supporting library of SCIP): ${SOPLEX_LIB}")
ELSE (SOPLEX_LIB)
	MESSAGE(STATUS "soplex (supporting library of SCIP) not found, no optimal bit heap compression available (download and compile scipoptsuite from http://scip.zib.de)")
ENDIF (SOPLEX_LIB)

IF (ZIMPL_LIB)
	MESSAGE(STATUS "Found zimpl (supporting library of SCIP): ${ZIMPL_LIB}")
ELSE (ZIMPL_LIB)
	MESSAGE(STATUS "zimpl (supporting library of SCIP) not found, no optimal bit heap compression available (download and compile scipoptsuite from http://scip.zib.de)")
ENDIF (ZIMPL_LIB)

IF (NLPI_LIB)
	MESSAGE(STATUS "Found nlpi (supporting library of SCIP): ${NLPI_LIB}")
ELSE (NLPI_LIB)
	MESSAGE(STATUS "nlpi (supporting library of SCIP) not found, no optimal bit heap compression available (download and compile scipoptsuite from http://scip.zib.de)")
ENDIF (NLPI_LIB)

IF (LPISPX_LIB)
	MESSAGE(STATUS "Found lpispx (supporting library of SCIP): ${LPISPX_LIB}")
ELSE (LPISPX_LIB)
	MESSAGE(STATUS "lpispx (supporting library of SCIP) not found, no optimal bit heap compression available (download and compile scipoptsuite from http://scip.zib.de)")
ENDIF (LPISPX_LIB)

# Detect GMP / GMPXX
FIND_PATH(GMPXX_H
	  gmpxx.h
)

FIND_LIBRARY(GMP_LIB
	NAMES gmp gmp.lib
	DOC "Directory of the GMP library")

FIND_PATH(GMP_H
	  gmp.h
)

FIND_LIBRARY(GMPXX_LIB
	NAMES gmpxx gmpxx.lib
	DOC "Directory of the GMPXX library")

IF (GMPXX_H AND GMPXX_LIB AND GMP_H AND GMP_LIB)
    MESSAGE(STATUS "Found gmp and gmpxx: ${GMP_H},${GMP_LIB},${GMPXX_H},${GMPXX_LIB}")
ELSE (GMPXX_H AND GMPXX_LIB AND GMP_H AND GMP_LIB)
    MESSAGE(FATAL_ERROR "Could not find gmp and/or gmpxx (found: ${GMP_H},${GMP_LIB},${GMPXX_H},${GMPXX_LIB})")
ENDIF (GMPXX_H AND GMPXX_LIB AND GMP_H AND GMP_LIB)


IF (SCIP_FOUND)
	ADD_DEFINITIONS(-DHAVE_SCIP)
	MESSAGE(STATUS "Found SCIP: ${SCIP_H} ${SCIP_LIB}")
ELSE (SCIP_FOUND)
	MESSAGE(STATUS "SCIP MILP optimizer not found (or supporting libraries missing), no optimal bit heap compression available (download and compile scipoptsuite from http://scip.zib.de)")
ENDIF (SCIP_FOUND)

SET(SCALP_INCLUDE_DIRS ${SCALP_H})
IF (SCALP_FOUND)
	ADD_DEFINITIONS(-DHAVE_SCALP)
    set(SCALP_LIBS ${SCALP_LIBS} ${SCALP_LIB}) 
	IF (SCALP_CPLEX_LIB)
		set(SCALP_LIBS ${SCALP_LIBS} ${SCALP_CPLEX_LIB})
	ENDIF (SCALP_CPLEX_LIB)
	IF (SCALP_GUROBI_LIB)
		set(SCALP_LIBS ${SCALP_LIBS} ${SCALP_GUROBI_LIB})
	ENDIF (SCALP_GUROBI_LIB)
	
	IF (SCALP_SCIP_LIB)
		set(SCALP_LIBS ${SCALP_LIBS} ${SCALP_SCIP_LIB})
	ENDIF (SCALP_SCIP_LIB)
	
	IF (SCALP_LPSOLVE_LIB)
		set(SCALP_LIBS ${SCALP_LIBS} ${SCALP_LPSOLVE_LIB})
	ENDIF (SCALP_LPSOLVE_LIB)
	
	MESSAGE(STATUS "Found SCALP: ${SCALP_INCLUDE_DIRS} ${SCALP_LIB}")
	MESSAGE(STATUS "      Using the following ILP Solvers: ${SCALP_LIBS}")  
ELSE (SCALP_FOUND)
	MESSAGE(STATUS "Could not find SCALP, no ILP solver available (used for advanced bit heap compression and large multipliers), download and compile from https://digidev.digi.e-technik.uni-kassel.de/scalp/ (${SCALP_INCLUDE_DIRS} ${SCALP_LIB})")
ENDIF (SCALP_FOUND)


# Detect pag_suite

SET(PagSuite_DIR ${CMAKE_CURRENT_SOURCE_DIR})
FIND_PACKAGE(PagSuite)

IF( PAGSUITE_FOUND )
	ADD_DEFINITIONS(-DHAVE_PAGLIB)
ELSE(PAGSUITE_FOUND)
	Message( "Could not find PagSuite. All related Operators (ConstMultPAG, ...) will not be build." )
ENDIF( PAGSUITE_FOUND ) 
#Uni KS stop

ADD_DEFINITIONS(-DHAVE_LNS)

# Create custom command for flex++/lex (note the outputs)
FIND_PROGRAM(FLEXPP_EXECUTABLE
	NAMES flex
	PATHS /usr/bin
	DOC "Flex++ executable")

IF(FLEXPP_EXECUTABLE)
	MESSAGE(STATUS "Found lex: ${FLEXPP_EXECUTABLE}")
ELSE(FLEXPP_EXECUTABLE)
  	MESSAGE(FATAL_ERROR "Could not find LEX.\nOn Debian-based systems install flex++")
ENDIF(FLEXPP_EXECUTABLE)

# necessary to include generated files
#INCLUDE_DIRECTORIES (${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/src/FPExpressions ${PAGSUITE_INCLUDE_DIRS})
INCLUDE_DIRECTORIES (${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/src/FPExpressions  ${PAGSUITE_INCLUDE_DIRS} /usr/include/boost148 ${SCALP_INCLUDE_DIRS} ${GMP_H} ${GMPXX_H} ${MPFI_H} ${MPFR_H} ${SOLLYA_H})
ADD_CUSTOM_COMMAND(
	OUTPUT VHDLLexer.cpp
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/VHDLLexerRule.l
	COMMAND ${FLEXPP_EXECUTABLE} 
	ARGS -o ${CMAKE_CURRENT_BINARY_DIR}/VHDLLexer.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/VHDLLexerRule.l
		)

#SET_SOURCE_FILES_PROPERTIES(VHDLLexer GENERATED) 


# Create custom command for bison/yacc (note the outputs)
FIND_PROGRAM(BISON_EXECUTABLE
	NAMES bison
	PATHS /usr/bin
	DOC "Bison executable")

IF(BISON_EXECUTABLE)
	MESSAGE(STATUS "Found bison: ${BISON_EXECUTABLE}")
ELSE(BISON_EXECUTABLE)
  	MESSAGE(FATAL_ERROR "Could not find Bison.\nOn Debian-based systems install bison")
ENDIF(BISON_EXECUTABLE)

ADD_CUSTOM_COMMAND(
	OUTPUT ExpressionParser.cpp
	COMMAND ${BISON_EXECUTABLE}
	ARGS -y -d --defines=${CMAKE_CURRENT_BINARY_DIR}/ExpressionParser.h --output=${CMAKE_CURRENT_BINARY_DIR}/ExpressionParser.cpp --name-prefix=FlopocoExpression
    ${CMAKE_CURRENT_SOURCE_DIR}/src/FPExpressions/ExpressionParser.y
)

SET_SOURCE_FILES_PROPERTIES(ExpressionParser GENERATED) 

# Create custom command for flex/lex (note the outputs)
FIND_PROGRAM(FLEX_EXECUTABLE
	NAMES flex
	PATHS /usr/bin
	DOC "Flex executable")

IF(FLEX_EXECUTABLE)
	MESSAGE(STATUS "Found flex: ${FLEX_EXECUTABLE}")
ELSE(FLEX_EXECUTABLE)
  	MESSAGE(FATAL_ERROR "Could not find LEX.\nOn Debian-based systems install flex++")
ENDIF(FLEX_EXECUTABLE)

ADD_CUSTOM_COMMAND(
	OUTPUT ExpressionScanner.cpp
	COMMAND ${FLEX_EXECUTABLE} 
	ARGS --prefix=FlopocoExpression --header-file=${CMAKE_CURRENT_BINARY_DIR}/FlopocoScanner.h -o ${CMAKE_CURRENT_BINARY_DIR}/ExpressionScanner.cpp 
    ${CMAKE_SOURCE_DIR}/src/FPExpressions/ExpressionScanner.l
	)

SET_SOURCE_FILES_PROPERTIES(ExpressionScanner GENERATED)

IF( PAGSUITE_FOUND )
	SET( PAGSUITE_DEP_SOURCES 
		src/ConstMultPAG/ConstMultPAG
		src/ConstMultPAG/ConstMultPAG_types
		src/ConstMultPAG/IntConstMultOpt
		src/ConstMultPAG/IntConstMultOptTernary
		src/ConstMultPAG/scm_solutions
		src/ConstMultPAG/tscm_solutions
		src/ConstMultPAG/pagexponents
		src/ConstMultPAG/FullyParallelFFT
	)
ELSE( PAGSUITE_FOUND )
	SET( PAGSUITE_DEP_SOURCES 
	)
ENDIF( PAGSUITE_FOUND )



ADD_LIBRARY(FloPoCoLib STATIC

#--------------------------------------------------------------------
#               FRAMEWORK
#--------------------------------------------------------------------

# src/FloPoCo.hpp
 src/Operator
 src/UserInterface
 src/Signal
 src/utils
 src/FlopocoStream
 src/TestState
 VHDLLexer # generated from .l
 
 # Resource Estimation related ---------------------------------------
 src/Tools/ResourceEstimationHelper
 
 # Floorplanning related ---------------------------------------------
 src/Tools/FloorplanningHelper

 # Targets -----------------------------------------------------------
 src/Targets/DSP
 src/Target 
 src/TargetModel
 
 src/Targets/Spartan3 
 src/Targets/Virtex4 
 src/Targets/Virtex5 
 src/Targets/Virtex6 
 
 src/Targets/StratixII 
 src/Targets/StratixIII 
 src/Targets/StratixIV
 src/Targets/StratixV
 
 src/Targets/CycloneII 
 src/Targets/CycloneIII 
 src/Targets/CycloneIV 
 src/Targets/CycloneV 

# TestBench-related -------------------------------------------------
 src/TestBenches/TestCase
 src/TestBenches/FPNumber
 src/TestBenches/IEEENumber
 src/TestBenches/Wrapper
 src/TestBenches/TestBench 


# Bit heap  ------------------------------------------------
 src/BitHeap/BitHeap
 src/BitHeap/WeightedBit
 src/BitHeap/Plotter
 src/BitHeap/BitHeapHeuristicCompression
 src/BitHeap/BitHeapILPCompression

 
#--------------------------------------------------------------------
#               OPERATORS          
#--------------------------------------------------------------------
#src/OperatorPipeline/OperatorPipeline
#src/FPExpressions/ExpressionParserData


# Conversions ----------------------------------------------------
src/Conversions/InputIEEE
src/Conversions/OutputIEEE
src/Conversions/FP2Fix
src/Conversions/Fix2FP

# Shifters/LZOC ----------------------------------------------------
 src/ShiftersEtc/LZOC 
 src/ShiftersEtc/LZOCShifterSticky
 src/ShiftersEtc/Shifters

# FixFilters ----------------------------------------------------------
 src/ShiftReg
 src/FixFilters/FixSOPC
 src/FixFilters/FixFIR
 src/FixFilters/FixHalfSine
 src/FixFilters/FixIIR


# pipelined Integer addition -------------------------------------------------
 src/IntAddSubCmp/IntAdder  
 src/IntAddSubCmp/IntAdderClassical 
 src/IntAddSubCmp/IntAdderAlternative 
 src/IntAddSubCmp/IntAdderShortLatency
 src/IntAddSubCmp/IntComparator
 src/IntAddSubCmp/IntDualSub
 src/IntAddSubCmp/BasicCompressor
 src/IntAddSubCmp/TargetOptCompressor
 src/IntAddSubCmp/VariableColumnCompressor
 src/IntAddSubCmp/FourToTwoCompressor
 src/IntAddSubCmp/IntAdderTree
 
# ---------- First generation of fast large adders ------------------------
 src/IntAddSubCmp/LongIntAdderAddAddMuxGen1 
 src/IntAddSubCmp/LongIntAdderCmpCmpAddGen1 
 src/IntAddSubCmp/LongIntAdderCmpAddIncGen1

# ---------- Second Generation of fast large adders -----------------------
 src/IntAddSubCmp/IntAdderSpecific 
 src/IntAddSubCmp/IntComparatorSpecific 
 src/IntAddSubCmp/CarryGenerationCircuit 
 src/IntAddSubCmp/LongIntAdderAddAddMuxGen2 
 src/IntAddSubCmp/LongIntAdderCmpCmpAddGen2 
 src/IntAddSubCmp/LongIntAdderCmpAddIncGen2
 src/IntAddSubCmp/LongIntAdderMuxNetwork

# Constant multiplications and divisions
 src/ConstMult/ShiftAddOp 
 src/ConstMult/ShiftAddDag 
 src/ConstMult/IntConstMult
 src/ConstMult/IntConstMCM  
 src/ConstMult/KCMTable
 src/ConstMult/IntIntKCM 
 src/ConstMult/FPConstMult 
 src/ConstMult/CRFPConstMult 
 src/ConstMult/FixRealKCM
 src/ConstMult/FPRealKCM
 src/ConstMult/IntConstDiv 
 src/ConstMult/FPConstDiv 

 
# Multiplication related --------------------------------------------
 src/IntMult/IntMultiplier
 src/IntMult/FixMultAdd
 src/IntMult/MultiplierBlock
 src/IntMult/IntSquarer
 src/IntMult/GenericBinaryPolynomial
 src/IntMult/IntPower
 src/IntMult/FormalBinaryProduct
#uni_ks
 src/IntMult/MultiplierSolutionParser
 src/IntMult/BaseMultiplierCollection
 src/IntMult/BaseMultiplier
 src/IntMult/BaseMultiplierLUT
 src/IntMult/BaseMultiplier2xk
 src/IntMult/BaseMultiplierDSP
 src/IntMult/BaseMultiplierDSPSuperTilesXilinx

# Complex numbers ---------------------------------------------------
src/Complex/FixComplexKCM
# TODO for Matei
#src/Complex/FixedComplexAdder

# Floating-point adders and subtracters -----------------------------
 src/FPAddSub/FPAdd
 src/FPAddSub/FPAddDualPath
 src/FPAddSub/FPAddSinglePath
 src/FPAddSub/FPAdd3Input
 src/FPAddSub/FPAddSub

# Floating-point multipliers and squarers
 src/FPMultSquare/FPMult
#  src/FPMultSquare/FPMultKaratsuba # See TODO
 src/FPMultSquare/FPSquare

# Floating-point dividers and square roots
src/FPDivSqrt/SelFunctionTable
src/FPDivSqrt/FPDiv
src/FPDivSqrt/FPSqrt
#src/FPDivSqrt/Tools/NbBitsMin
# src/FPDivSqrt/FPSqrtPoly # Resurrect some day?


# Composite operators
src/FPComposite/FPLargeAcc 
src/FPComposite/LargeAccToFP 
src/FPComposite/FPDotProduct 
src/FPComposite/FP2DNorm


# Generic fixed-point function evaluators ---------------------------
 src/Table 
 src/DualTable 

 src/FixConstant

 src/FixFunctions/FixFunction
 src/FixFunctions/FixFunctionByTable
 src/FixFunctions/BasicPolyApprox
 src/FixFunctions/PiecewisePolyApprox
 src/FixFunctions/FixHornerEvaluator
 src/FixFunctions/FixFunctionBySimplePoly
 src/FixFunctions/FixFunctionByPiecewisePoly

 src/FixFunctions/GenericTable
 src/FixFunctions/BipartiteTable
 src/FixFunctions/Multipartite
 src/FixFunctions/FixFunctionByMultipartiteTable
# src/FixFunctions/QuineMcCluskey

# Trigonometric functions ------------------------------------------
src/Trigs/ConstDiv3ForSinPoly
src/Trigs/FixSinPoly
src/Trigs/FixSinCos
# src/Trigs/CordicSinOrCos
#commented out due to strange error in mpfr (object expression of non-scalar type '__mpfr_struct [1]' cannot be used in a pseudo-destructor expression)
#src/Trigs/CordicSinCos
#src/Trigs/FixAtan2
#src/Trigs/FixAtan2ByCORDIC
src/Trigs/FixAtan2ByRecipMultAtan
src/Trigs/FixAtan2ByBivariateApprox
src/Trigs/Atan2Table

src/Trigs/Fix2DNorm


src/Tools/Plane
src/Tools/Point

# Exp, Log and Power
src/ExpLog/FPExp
src/ExpLog/IterativeLog
src/ExpLog/FPPow

# Primitive Components ----------------------------------------------
src/PrimitiveComponents/bool_eq
src/PrimitiveComponents/GenericLut
src/PrimitiveComponents/Primitive
src/PrimitiveComponents/GenericAddSub
src/PrimitiveComponents/GenericMux
src/PrimitiveComponents/Tests/GenericLut_Test

src/PrimitiveComponents/Xilinx/Xilinx_Primitive
src/PrimitiveComponents/Xilinx/Xilinx_LUT_compute
src/PrimitiveComponents/Xilinx/Xilinx_CARRY4
src/PrimitiveComponents/Xilinx/Xilinx_LUT5
src/PrimitiveComponents/Xilinx/Xilinx_LUT6
src/PrimitiveComponents/Xilinx/Xilinx_MUXF7
src/PrimitiveComponents/Xilinx/Xilinx_MUXF8
src/PrimitiveComponents/Xilinx/Xilinx_FDCE
src/PrimitiveComponents/Xilinx/Xilinx_n2m_Decoder
src/PrimitiveComponents/Xilinx/Xilinx_Comparator
src/PrimitiveComponents/Xilinx/Xilinx_GenericAddSub
src/PrimitiveComponents/Xilinx/Xilinx_GenericAddSub_slice
src/PrimitiveComponents/Xilinx/Xilinx_GenericMux
src/PrimitiveComponents/Xilinx/Xilinx_GenericMux_slice
src/PrimitiveComponents/Xilinx/Xilinx_TernaryAdd_2State
src/PrimitiveComponents/Xilinx/Xilinx_TernaryAdd_2State_slice

src/PrimitiveComponents/Altera/Altera_LCELL

# Uni_ks Components ------------------------------------------------------
src/SAD/SAD
${PAGSUITE_DEP_SOURCES}

# If you want to add your operator, feel free
src/UserDefinedOperator
src/NeuralNetworks/ConfigurationGenerator
src/NeuralNetworks/ConfigurationController
src/NeuralNetworks/ConvolutionalCore
src/NeuralNetworks/SOP_KCM
src/NeuralNetworks/CoreTE
src/NeuralNetworks/LUTMult
src/PrimitiveComponents/Xilinx/Xilinx_CFGLUT5
src/PrimitiveComponents/Xilinx/Xilinx_CFGLUTShadow

# activation functions
src/NeuralNetworks/ActivationFunctions/ReLU

# layers
src/NeuralNetworks/Layers/Layer
#   convolutional
src/NeuralNetworks/Layers/Convolutional/ConvolutionalLayer
src/NeuralNetworks/Layers/Convolutional/ConvolutionalCoreSimple
#   full connected
src/NeuralNetworks/Layers/FullConnected/FullConnectedLayer
#   pooling
src/NeuralNetworks/Layers/Pooling/PoolingLayer
src/NeuralNetworks/Layers/Pooling/PoolingCore
src/NeuralNetworks/Layers/Pooling/Max2Input
#   in/out
src/NeuralNetworks/Layers/Input/InputLayer
src/NeuralNetworks/Layers/Output/OutputLayer


# memory management
src/NeuralNetworks/MemoryManagement/BlockRam
src/NeuralNetworks/MemoryManagement/ControlledMemory
src/NeuralNetworks/MemoryManagement/MemoryController

# utility
src/NeuralNetworks/Utility/WindowShiftRegister
src/NeuralNetworks/Utility/GlobalController
src/NeuralNetworks/Utility/Register
src/NeuralNetworks/Utility/ModuloCounter
src/NeuralNetworks/Utility/PaddingGenerator

# general
src/NeuralNetworks/H5ToNeuralNetwork
src/NeuralNetworks/NeuralNetwork
src/NeuralNetworks/LayerArguments





)

set_target_properties(FloPoCoLib PROPERTIES OUTPUT_NAME "FloPoCo")

TARGET_LINK_LIBRARIES(
  FloPoCoLib 
  ${GMP_LIB} ${GMPXX_LIB} ${MPFI_LIB} ${MPFR_LIB} #mpfr gmp gmpxx xml2 mpfi
  )

IF (SOLLYA_LIB)
TARGET_LINK_LIBRARIES(
  FloPoCoLib
  ${SOLLYA_LIB}
	)
ENDIF (SOLLYA_LIB)

IF (PYTHON_LIB)
TARGET_LINK_LIBRARIES(
  FloPoCoLib
  ${PYTHON_LIB}
  )
ENDIF (PYTHON_LIB)

#Uni KS start
IF (SCIP_FOUND)
TARGET_LINK_LIBRARIES(
  FloPoCoLib
  ${SCIP_LIB} ${READLINE_LIB} ${Z_LIB} ${NLPI_LIB} ${LPISPX_LIB} ${SOPLEX_LIB} ${GMP_LIB} ${ZIMPL_LIB}
  )
ENDIF (SCIP_FOUND)

IF (SCALP_FOUND)
TARGET_LINK_LIBRARIES(
  FloPoCoLib
  ${SCALP_LIBS}
  )
ENDIF (SCALP_FOUND)
#Uni KS stop


IF (WCPG_FOUND)
TARGET_LINK_LIBRARIES(
  FloPoCoLib
  ${WCPG_LIB}
  )
ENDIF (WCPG_FOUND)

IF(PAGSUITE_FOUND)
TARGET_LINK_LIBRARIES(
	FloPoCoLib
	${PAGSUITE_LIBRARIES}
)
ENDIF(PAGSUITE_FOUND)


ADD_EXECUTABLE(flopoco
 src/main 
)
 
TARGET_LINK_LIBRARIES(
  flopoco  
	FloPoCoLib 
  )

# unplugged for now to speed up compilation, to revive some day
#ADD_EXECUTABLE(flopoco_soaktest
# src/main_soaktest
#)
#TARGET_LINK_LIBRARIES(
#  flopoco_soaktest  
#	${FLOPOCO_LIB}
#)

 
ADD_EXECUTABLE(fp2bin src/Tools/fp2bin  src/utils)
TARGET_LINK_LIBRARIES(fp2bin ${MPFR_LIB} ${GMP_LIB} ${GMPXX_LIB})

ADD_EXECUTABLE(bin2fp src/Tools/bin2fp  src/utils)
TARGET_LINK_LIBRARIES(bin2fp ${MPFR_LIB} ${GMP_LIB} ${GMPXX_LIB})

ADD_EXECUTABLE(longacc2fp src/Tools/longacc2fp  src/utils)
TARGET_LINK_LIBRARIES(longacc2fp ${MPFR_LIB} ${GMP_LIB} ${GMPXX_LIB})



#add_subdirectory(src/random)
#target_link_libraries(FloPoCoLib flopoco_random_lib)
