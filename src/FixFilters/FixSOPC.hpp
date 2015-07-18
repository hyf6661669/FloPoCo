#ifndef FixSOPC_HPP
#define FixSOPC_HPP

#include "Operator.hpp"
#include "utils.hpp"

#include "BitHeap/BitHeap.hpp"

/*  All flopoco operators and utility functions are declared within
  the flopoco namespace.
    You have to use flopoco:: or using namespace flopoco in order to access these
  functions.
*/

namespace flopoco{

	class FixSOPC : public Operator {
	public:

		/** simplest constructor for inputs in the fixed-point format (0, lsbIn), computing msbOut out of the coeffs, computing the internal format.
		 This constructor is all we need for a FIR */
		FixSOPC(Target* target, int lsbIn, int lsbOut, std::vector<std::string> coeff);


		/** Generic constructor for inputs in various formats and/or for splitting a SOPC into several ones, etc. 
				msbOut must be provided. 
				If g=-1, the number of needed guard bits will be computed for a faithful result, and a final round bit added in position lsbOut-1. 
				If g=0, the architecture will have no guard bit, no final round bit will be added. The architecture will not be faithful. 
				If g>0, the provided number of guard bits will be used and a final round bit added in position lsbOut-1.
 */
		FixSOPC(Target* target, std::vector<int> &msbIn, std::vector<int> &lsbIn, int msbOut, int lsbOut, std::vector<std::string> coeff_, int g=-1);

		/** destructor */
		~FixSOPC();

		/** The method that does most of operator construction for the two constructors */
		void initialize();

		/** Overloading the method of Operator */
		void emulate(TestCase * tc);

		/** Overloading the method of Operator */
		void buildStandardTestCases(TestCaseList* tcl);


		/** Returns an interval in which to find the result, given a std::vector (x) of inputs.
		  returns a result in very large precision ("6400 bits should be enough for anybody")
		  */
		int computeSOPCForEmulate(std::vector<mpz_class> inputs, mpfr_t &s);

		/** This method does most of the work for emulate(), because we want to call it also from the emulate() of FixFIR
			Returns an interval in which to find the result, given a std::vector (x) of inputs.
		  */
		std::pair<mpz_class,mpz_class> computeSOPCForEmulate(std::vector<mpz_class> x);

	protected:
		int n;							        /**< number of products, also size of the std::vectors coeff, msbIn and lsbIn */
		std::vector<int> msbIn;			    /**< MSB weights of the inputs */
		std::vector<int> lsbIn;			    /**< LSB weights of the inputs */
	public: // readable by FIR etc
		int msbOut;							    /**< MSB weight of the output, may be computed out of the constants (depending on the constructor used) */
		int lsbOut;							    /**< LSB weight of the output */
	protected:
		std::vector<std::string> coeff;			  /**< the coefficients as strings */
		mpfr_t mpcoeff[10000];			/**< the coefficients as MPFR numbers -- 10000 should be enough for anybody */
		int g;                      /**< Number of guard bits; the internal format will have LSB at lsbOut-g  */


	private:
		bool computeMSBOut;     /** <*/
		bool computeGuardBits;     /** <*/
		bool addFinalRoundBit;     /** <*/
		BitHeap* bitHeap;    			 /**< The heap of weighted bits that will be used to do the additions */
	};


}

#endif
