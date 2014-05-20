set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})

message(${CMAKE_MODULE_PATH})

FIND_PACKAGE(Boost COMPONENTS test_exec_monitor)
#FIND_PACKAGE(GSL)

FIND_LIBRARY(GSL_LIBRARIES NAMES gsl )

FIND_PATH(NTL_H NTL/RR.h )

FIND_LIBRARY(NTL_LIB NAMES ntl ntl.lib )

IF (NTL_H AND NTL_LIB)
    SET(NTL_FOUND TRUE)
	ADD_DEFINITIONS(-DHAVE_NTL)
	
	FIND_LIBRARY(GMP_LIB NAMES gmp gmp.lib )

	INCLUDE_DIRECTORIES(${NTL_H})
	LINK_DIRECTORIES(${NTL_LIB}  )
ENDIF (NTL_H AND NTL_LIB)

FIND_PATH(MPREAL_H mpreal.h)

IF(MPREAL_H)
	INCLUDE_DIRECTORIES(${MPREAL_H})
ELSE(MPREAL_H)
	INCLUDE_DIRECTORIES(utils/mpreal)
ENDIF(MPREAL_H)

#ADD_LIBRARY(flopoco_random_asa_lib STATIC
#	utils/asa/asa.c
#	utils/asa/asa_usr.c
#)
#set_target_properties(flopoco_random_asa_lib	PROPERTIES
#	COMPILE_DEFINITIONS "ASA_LIB=TRUE"
#	COMPILE_DEFINITIONS "USER_COST_FUNCTION=asa_user_cost_function"
#	COMPILE_DEFINITIONS "OPTIONAL_DATA_PTR=TRUE"
#	COMPILE_DEFINITIONS "OPTIONAL_PTR_TYPE=\"void *\""
#	COMPILE_DEFINITIONS "INT_LONG=TRUE"
#	COMPILE_DEFINITIONS "INT_ALLOC=FALSE"
#)

ADD_LIBRARY(flopoco_random_lib STATIC
	random_main.cpp
	utils/operator_factory.cpp	
	utils/comparable_float_type.cpp
	utils/polynomial_evaluator_factory.cpp
	utils/output_combiner.cpp
	
	poly/fixed_point_polynomial_evaluator.cpp
	poly/full_precision_polynomial_evaluator.cpp
	poly/fixed_format_t.cpp
	
	urng/lut_sr_rng.cpp
	urng/mwc_rng.cpp
	urng/mwc_transform.cpp
	
	transforms/blocks/HadamardTransform.cpp
	transforms/blocks/CLTTransform.cpp
	transforms/blocks/TableTransform.cpp
	transforms/blocks/GRNGTableTransform.cpp
	transforms/blocks/CLTCorrectedTransform.cpp
	transforms/blocks/OutputShuffle.cpp
	transforms/blocks/OutputShuffleTransform.cpp
	transforms/TransformStats.cpp
	
	float_approx/static_quantiser.cpp
	float_approx/range.cpp
	float_approx/range_polys.cpp
	float_approx/FloatApprox.cpp
)
#target_link_libraries(flopoco_random_lib ${GMP_LIB} ${FLOPOCO_LIB}  ${GSL_LIBRARIES})
target_link_libraries(flopoco_random_lib ${GMP_LIB} ${GSL_LIBRARIES} ${MPFR_LIB})


INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR} ${GSL_INCLUDE_DIR})
LINK_DIRECTORIES ( ${Boost_LIBRARY_DIRS} ${GSL_LIBRARY_DIRS})

#add_executable(test_fixed_point_exp fixed_point_exp/test_fixed_point_exp.cpp)
#target_link_libraries(test_fixed_point_exp ${Boost_LIBRARIES} ${GSL_LIBRARIES} ${NTL_LIB} ${GMP_LIB} FloPoCoLib)

#add_executable(test_accumulation utils/test_accumulation.cpp)
#target_link_libraries(test_accumulation ${Boost_LIBRARIES} ${GSL_LIBRARIES} ${NTL_LIB} ${GMP_LIB})

#add_executable(test_distributions distributions/test_distributions.cpp)
#target_link_libraries(test_distributions ${Boost_LIBRARIES} ${GSL_LIBRARIES} ${MPFR_LIB} ${GMP_LIB})

#add_executable(test_optimisation utils/test_optimisation.cpp)
#target_link_libraries(test_optimisation ${Boost_LIBRARIES} ${GSL_LIBRARIES})

add_executable(test_moment_correction moment_correction/test_moment_correction.cpp)
target_link_libraries(test_moment_correction ${Boost_LIBRARIES} ${GSL_LIBRARIES} ${NTL_LIB} ${GMP_LIB} ${MPFR_LIB})

#add_executable(test_moment_correction_weighted moment_correction/test_moment_correction_weighted.cpp)
#target_link_libraries(test_moment_correction_weighted ${Boost_LIBRARIES} ${GSL_LIBRARIES} ${NTL_LIB} ${GMP_LIB}  ${MPFR_LIB})


#add_executable(gaussian_moment_correction_tables moment_correction/gaussian_moment_correction_tables.cpp)
#target_link_libraries(gaussian_moment_correction_tables ${Boost_LIBRARIES} ${GSL_LIBRARIES} ${NTL_LIB} ${GMP_LIB} ${MPFR_LIB})

#add_executable(find_hadamard_table moment_correction/find_hadamard_table.cpp)
#target_link_libraries(find_hadamard_table ${Boost_LIBRARIES} ${GSL_LIBRARIES} ${NTL_LIB} ${GMP_LIB})


#add_executable(test_dag_parser dag/test_dag_parser.cpp dag/dag_parser.cpp dag/dag_printer.cpp)
#target_link_libraries(test_dag_parser ${Boost_LIBRARIES})


#add_executable(test_comparable_float_type utils/test_comparable_float_type.cpp)
#target_link_libraries(test_comparable_float_type ${Boost_LIBRARIES} ${GMP_LIB} ${MPFR_LIB})

#add_executable(test_range float_approx/test_range.cpp float_approx/range.cpp)
#target_link_libraries(test_range flopoco_random_lib ${FLOPOCO_LIB})

#add_executable(test_fft utils/fft/test_fft.cpp)
#target_link_libraries(test_fft ${MPFR_LIB} ${GSL_LIBRARIES} ${GMP_LIB})

#add_executable(test_poly_min_max utils/test_poly_min_max.cpp)
#target_link_libraries(test_poly_min_max FloPoCoLib ${Boost_LIBRARIES})

#add_executable(test_decompose_period urng/test_decompose_period.cpp)
#target_link_libraries(test_decompose_period ${GMP_LIB})



