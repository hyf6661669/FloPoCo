/*
 * MultiFunctApprox.hpp
 *
 *  Author: Matei Istoan
 *  2019
 */

#ifndef _MULTIFUNCTAPPROX_HPP_
#define _MULTIFUNCTAPPROX_HPP_

#include <string>
#include <iostream>
#include <sstream>
#include <limits.h>
#include <float.h>

#include <sollya.h>
#include <gmpxx.h>

#include "Operator.hpp" // mostly for reporting
#include "FixFunction.hpp"
#include "FixConstant.hpp"
#include "BasicPolyApprox.hpp"
#include "PiecewisePolyApprox.hpp"

using namespace std;

/* Stylistic convention here: all the sollya_obj_t have names that end with a capital S */
namespace flopoco {

	class MultiFunctApprox {
	public:

		/**
		 * A minimal constructor
		 */
		MultiFunctApprox(vector<FixFunction*> functs, double targetAccuracy, int degree);

		/**
		 * A minimal constructor that parses a sollya string
		 */
		MultiFunctApprox(vector<string> sollyaStringList, double targetAccuracy, int degree);

		/**
		 * class destructor
		 */
		virtual ~MultiFunctApprox();

		/**
		 * parse the command line arguments
		 */
		static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args);

		/**
		 * register the operator for use on the command line
		 */
		static void registerFactory();

		/**
		 * function that regroups most of the constructor code
		 */
		void build();

		int degree;                               /**< max degree of the polynomial approximations */
		int alpha;                                /**< the input domain [0,1] will be split into at most 2^alpha subdomains */
		vector<BasicPolyApprox*> functApprox;     /**< The vector of polynomial approximations, each with its format */
		int LSB;                                  /**< min weight of the LSBs of the polynomial approximations */
		vector<int> MSB;                          /**< vector of the max MSB weights for each coefficient */
		double approxErrorBound;                  /**< guaranteed upper bound on the approx error of each approximation provided. Should be smaller than targetAccuracy */

	private:

		vector<FixFunction*> functs;              /**< The function to be approximated */
		double targetAccuracy;                    /**< please build an approximation at least as accurate as that */

		string srcFileName;                       /**< useful only to enable same kind of reporting as for FloPoCo operators. */
		string uniqueName_;                       /**< useful only to enable same kind of reporting as for FloPoCo operators. */
		bool needToFreeF;                         /**< in an ideal world, this should not be needed */

		fstream cacheFile;                        /**< file storing the cached parameters for the polynomials */
		vector<int> nbIntervals;                  /**< the total number of intervals the domain is split into */
		vector<string> cacheFileNames;            /**< the name of the file used for storing the cached polynomial parameters */
	};

} /* namespace flopoco */

#endif /* _MULTIFUNCTAPPROX_HPP_ */
