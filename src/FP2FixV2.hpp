#ifndef FP2FixV2_HPP
#define FP2FixV2_HPP
#include <vector>
#include <sstream>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>

#include "Operator.hpp"
#include "FPNumber.hpp"
#include "IntAdder.hpp"


namespace flopoco{

   /** The FP2Fix class */
   class FP2FixV2 : public Operator
   {
   public:
      /**
		 * The  constructor
		 * @param[in]		target		the target device
                 * @param[in]		MSB			the MSB of the output number for the convertion result
                 * @param[in]		LSB			the LSB of the output number for the convertion result
                 * @param[in]		wER			the with of the exponent in input
                 * @param[in]		wFR			the with of the fraction in input
                 * @param[in]		trunc_p			the output is not rounded when trunc_p is true
                 */
      FP2FixV2(Target* target, int LSBO, int MSBO, int Signed,int wER, int wFR, bool trunc_p);

      /**
		 *  destructor
		 */
      ~FP2FixV2();


      void emulate(TestCase * tc);
      
      /* Make sure we exercise numbers in the right range, as well as pure random*/
      TestCase* buildRandomTestCase(int i);
      
      void buildStandardTestCases(TestCaseList* tcl);

   private:

      /** The width of the exponent for the input */
      int wEI;
      /** The width of the fraction for the input */
      int wFI;
      /** are all numbers positive or not */
      int Signed;
      /** The LSB for the output */
      int LSBO;
      /** The MSB for the output */
      int MSBO; 
      /** when true the output is not rounded */
      bool trunc_p;
   };
}
#endif
