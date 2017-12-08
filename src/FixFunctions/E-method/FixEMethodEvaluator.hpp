/*

  A hardware implementation of the E-method for the evaluation of polynomials and rational polynomials.

  Author : Matei Istoan

*/

#ifndef FIXEMETHODEVALUATOR_HPP_
#define FIXEMETHODEVALUATOR_HPP_

#include <vector>
#include <sstream>
#include <string>

#include <sollya.h>
#include <gmpxx.h>

#include "Operator.hpp"
#include "Table.hpp"
#include "IntMult/FixMultAdd.hpp"
#include "BitHeap/BitHeap.hpp"
#include "FixFunction.hpp"
#include "FixConstant.hpp"

namespace flopoco {

	/**
	 * A hardware implementation of the E-method for the evaluation of polynomials and rational polynomials.
	 * Computing:
	 * 		y_0 = R(x) = P(x)/Q(x)
	 * , where:
	 * 		P(x) = p_n*x^n + p_{n-1}*x^{n-1} + ... + p_1*x^1 + p_0
	 * 		Q(x) = q_m*x^m + q_{m-1}*x^{m-1} + ... + q_1*x^1 + 1
	 * 		so q_0 = 1
	 * Currently only supporting radix 2.
	 * 		Planned support for higher radixes (4, 8 etc.).
	 * Format of the input X is given as a parameter.
	 * Format of the coefficients is given as a parameter.
	 * Checks can, and should, be made to insure that all parameters are in ranges where the method converges.
	*/

  class FixEMethodEvaluator : public Operator
  {
  public:
    /**
     * A constructor that exposes all options.
     * @param   n              degree of the polynomial P
     * @param   m              degree of the polynomial Q
     * @param   msbIn          MSB of the input
     * @param   lsbIn          LSB of the input
     * @param   msbOut         MSB of the output
     * @param   lsbOut         LSB of the output
     * @param   coeffsP        vector holding the coefficients of polynomial P
     * @param   coeffsQ        vector holding the coefficients of polynomial Q
     */
	FixEMethodEvaluator(Target* target,
			  int n,
			  int m,
			  int msbIn,
			  int lsbIn,
			  int msbOut,
			  int lsbOut,
			  vector<mpfr_t> coeffsP,
			  vector<mpfr_t> coeffsQ,
			  map<string, double> inputDelays = emptyDelayMap);

	/**
	 * Class destructor
	 */
    ~FixEMethodEvaluator();


  private:
    int n;                            /**< degree of the polynomial P */
    int m;                            /**< degree of the polynomial Q */
    int msbIn;                        /**< MSB of the input */
    int lsbIn;                        /**< LSB of the input */
    int msbOut;                       /**< MSB of the output  */
    int lsbOut;                       /**< LSB of the output */
    vector<mpfr_t> coeffsP;           /**< vector of the coefficients of P */
    vector<mpfr_t> coeffsQ;           /**< vector of the coefficients of Q */

    int maxDegree;                    /**< the maximum between the degrees of the polynomials P and Q */
    int nbIter;                       /**< the number of iterations */
    int g;                            /**< number of guard bits */
  };

} /* namespace flopoco */

#endif /* _FIXEMETHODEVALUATOR_HPP_ */
