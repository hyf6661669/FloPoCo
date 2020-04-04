/*
  
  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Author : Florent de Dinechin, Florent.de-Dinechin@insa-lyon.fr

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL, INSA,
  2008-2014.
  All rights reserved.

*/
#ifndef __FIXHORNEREVALUATOR_HPP
#define __FIXHORNEREVALUATOR_HPP
#include <vector>
#include <sstream>

#include "Operator.hpp"
#include "FixPolyEval.hpp"


namespace flopoco{

	/** An Horner polynomial evaluator computing just right.
	 It assumes the input X is an signed number in [-1, 1[ so msbX=-wX.
	*/

  class FixHornerEvaluator : public FixPolyEval
  {
  public:



		/** The constructor.
     * @param    lsbIn input lsb weight, 
			 @param    msbOut  output MSB weight, used to determine wOut
			 @param    lsbOut  output LSB weight
			 @param    poly the vector of polynomials that this evaluator should accomodate
			 @param   finalRounding: if false, the operator outputs its guard bits as well, saving the half-ulp rounding error. 
			                 This makes sense in situations that further process the result with further guard bits.
     */

																				
    FixHornerEvaluator(OperatorPtr parentOp, Target* target,
											 int lsbIn,
											 int msbOut,
											 int lsbOut,
											 vector<BasicPolyApprox*> p,
											 bool finalRounding=true);


#if 0


		/** The constructor with manual control of all options.
     * @param    lsbIn input lsb weight, 
			 @param    msbOut  output MSB weight, used to determine wOut
			 @param    lsbOut  output LSB weight
			 @param    degree  degree of the polynomial
			 @param    msbCoeff vector (of size degree+1) holding the MSB of the polynomial coefficients
			 @param    lsbCoeff vector (of size degree+1) holding the LSB of the polynomial coefficients
			 @param    roundingErrorBudget The rounding error budget, excluding final rounding. If -1, will be set to 2^(lsbOut-2) 
			 @param    signedXandCoeffs  true if the coefficients are signed numbers (usually true)
			 @param   finalRounding: if false, the operator outputs its guard bits as well, saving the half-ulp rounding error. 
			                 This makes sense in situations that further process the result with further guard bits.
     */

																				
    FixHornerEvaluator(OperatorPtr parentOp, Target* target, 
											 int lsbIn,
											 int msbOut,
											 int lsbOut,
											 int degree, 
											 vector<int> msbCoeff, 
											 int lsbCoeff,
											 double roundingErrorBudget=-1,
											 bool signedXandCoeffs=true, 
											 bool finalRounding=true);

		
		/** An optimized constructor if the caller has been able to compute the signs and MSBs of the sigma terms */
    FixHornerEvaluator(OperatorPtr parentOp, Target* target, 
											 int lsbIn,
											 int msbOut,
											 int lsbOut,
											 int degree, 
											 vector<int> msbCoeff, 
											 int lsbCoeff,
											 vector<int> sigmaSign, vector<int> sigmaMSB,
											 double roundingErrorBudget=-1,
											 bool signedXandCoeffs=true, 
											 bool finalRounding=true);

#endif

		
    ~FixHornerEvaluator();
		

  private: // we also inherit attribute of FixPolyEval
		// internal architectural parameters; 
		vector<double> wcMaxAbsSum; /**< from 0 to degree */
		vector <int> wcSumSign; /**< 1: always positive; -1: always negative; 0: can be both  */  
		vector<int> wcSumMSB; /**< from 0 to degree */
		vector<int> wcSumLSB; /**< from 0 to degree */
		vector<int> wcYLSB; /**< from 0 to degree */
		// vector <int>  wcProductMSB; /**< from 0 to degree */
		// vector <int>  wcProductLSB; /**< from 0 to degree-1 */

		void initialize(); /**< initialization factored out between various constructors */ 
		void computeArchitecturalParameters(); /**< error analysis that ensures the rounding budget is met */ 
		void generateVHDL(); /**< generation of the VHDL once all the parameters have been computed */ 
  };

}
#endif
