/*
*/

#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include <string.h>

#include <gmp.h>

#include <gmpxx.h>
#include "utils.hpp"
#include "Operator.hpp"

#include "FP2FixV2.hpp"

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <list>
#include <map>
#include <vector>
#include <math.h>
#include <locale>

#include <stdio.h>
#include <mpfr.h>

#include "Shifters.hpp"

using namespace std;

namespace flopoco{

#define DEBUGVHDL 0


   FP2FixV2::FP2FixV2(Target* target, int _LSBO, int _MSBO, int _Signed,int _wEI, int _wFI, bool _trunc_p) :
         Operator(target), wEI(_wEI), wFI(_wFI), Signed(_Signed), LSBO(_LSBO), MSBO(_MSBO), trunc_p(_trunc_p) {

      int MSB=MSBO;
      int LSB=LSBO;
      int Width=MSB-LSB+1;
            
      int wE=wEI, wF=wFI;
      
      ostringstream name;
      name<<"FP2FixV2_" << wE << "_" << wF <<"_"<<(LSB<0?"M":"")<<std::abs(LSB)<<"_"<<(MSB<0?"M":"")<<std::abs(MSB)<<"_"<<(Signed==1?"S":"US") << "_" << (trunc_p==1?"T":"NT")<<"_uid"<<getNewUId();
      setName(name.str());

      setCopyrightString("David B. Thomas & Fabrizio Ferrandi (2013)");
            
      addFPInput("I", _wEI, _wFI);
      addOutput("O", Width);

      if(_Signed){
         cerr<< " FP2FixV2 : Cannot handle signed numbers yet."<<endl;
         exit(EXIT_FAILURE);
      }
      if ((MSB < LSB)){
         cerr << " FP2FixV2: Input constraint LSB <= MSB not met."<<endl;
         exit (EXIT_FAILURE);
      }
      
      int bias=(1<<(wE-1))-1;
      // The smallest shift we care about is when input = 2^(LSB-1) * 1.xxxxxxx
      int minUsefulExpnt=(LSB-1) + bias;
      // And the biggest shift is when input = 2^MSB * 1.xxxxxxx
      int maxUsefulExpnt=MSB + bias;
      
      if(minUsefulExpnt>((1<<wE)-1)){
         std::cerr<<"For these parameters all inputs will overflow.\n";
         exit(1);
      }
      if(maxUsefulExpnt<0){
         std::cerr<<"For these parameters all inputs will underflow.\n";
         exit(1);
      }
      
      vhdl<<tab<<declare("I_flags",2)<<"<=I"<<range(2+wE+wF,1+wE+wF)<<";\n";
      vhdl<<tab<<declare("I_sign")<<"<=I("<<wE+wF<<");\n";
      vhdl<<tab<<declare("I_expnt",wE)<<"<=I"<<range(wE+wF-1, wF)<<";\n";
      vhdl<<tab<<declare("I_frac",wF)<<"<=I"<<range(wF-1, 0)<<";\n";
      
      ////////////////////////////////////////////////////
      // From here until further notice, we only care about
      // normal values in the valid input range. We'll patch
      // up abnormals at the end
      
      int ShiftAmountWidth=intlog2(Width);
      if(ShiftAmountWidth==wE){
         vhdl<<tab<<declare("shift_amount", wE)<<" <= I_expnt - ("<<minUsefulExpnt<<");\n";
      }else if(ShiftAmountWidth<wE){
         // Common case: fixed-point width is small compared to exponent range
         vhdl<<tab<<declare("shift_amount_raw", wE)<<" <= I_expnt - ("<<minUsefulExpnt<<");\n";
         vhdl<<tab<<declare("shift_amount", ShiftAmountWidth)<<" <= shift_amount_raw"<<range(ShiftAmountWidth-1,0)<<";\n";
      }else{
         // This means the output width is huge, or the exponent is very narrow. Likely to be inefficient.
         // It would be better handled by adding one to the exponent further up if necessary - at the very
         // least the user should think about it
         THROWERROR("The output range for FP2Fix is greater than the exponent range. If this is really what you want, you should increase the exponent width.");
         
      }
      syncCycleFromSignal("shift_amount");

      int ShifterLSB=LSB-wF-1;
      int ShifterMSB=MSB;
      int ShifterWidth=ShifterMSB-ShifterLSB+1;
      
      //           MSB      LSB LSB-2     LSB-wF-1
      //            |         | |           |
      //            |         | |           |
      // Set up:    00000...0001XXXX....XXXXX
      vhdl<<tab<<declare("shift_input", wF+1)<<" <= '1' & I_frac;\n";
      
      Operator *shifter=new Shifter(target, wF+1, Width, Shifter::Left);
      oplist.push_back(shifter);
      
      inPortMap(shifter, "S", "shift_amount");
      inPortMap(shifter, "X", "shift_input");
      outPortMap(shifter, "R", "shifted");
      vhdl<<tab<<instance(shifter, "shifter")<<"\n";
      syncCycleFromSignal("shifted");
      
      /////////////////////////////////////////////////////////////////
      // Mantissa has been put in the right place, now worry about rounding
      
      if(trunc_p){
         // No rounding needed, nice and easy
         vhdl<<tab<<declare("rounded", Width)<<" <= shifted"<<range(ShifterWidth-1,wF+1)<<";\n"; 
      }else{
         // We now need to round, and so there are two cases:
         // wF <= MSB-LSB+1 : It is possible for rounding to cause overflow
         // wF > MSB-LSB+1  : No overflow due to rounding is possible
         bool canOverflowOnRound=(wF<=Width);
         // Now round the values
         //              MSB     LSB        LSB-wF-1
         //              |         |            |
         // For overflow possible: |            |
         //              |         |            |
         // Min shift:  000000...0001xxxx....xxxx
         // Max shift:  01xxxx...xxxxx000....0000
         // +            |         |            |
         // Rounding:   000000...00001111....1111
         //              |         |            |
         // For overflow not possible:          |
         //              |         |            |
         // Min shift:   00000...0001xxxx....xxxx
         // Max shift:   1xxxx...xxxx0000....0000
         // +            |         |            |
         // Rounding:    00000...00001111....1111
         if(canOverflowOnRound){
            vhdl<<tab<<"-- Overflow possible while rounding.\n";
            vhdl<<tab<<declare("rounding_offset", ShifterWidth+1)<<" <= "<<zg(Width+2)<<"&"<<og(wF)<<";\n";
            vhdl<<tab<<declare("rounded_raw", ShifterWidth+1)<<" <= shifted + rounding_offset;\n";
            vhdl<<tab<<declare("rounded", Width)<<" <= "<<og(Width)<<" when rounded_raw("<<ShifterWidth<<")='1' "
                                   <<" else rounded_raw"<<range(ShifterWidth-1,wF+1)<<";\n";
         }else{
            vhdl<<tab<<"-- Overflow not possible while rounding.\n";
            vhdl<<tab<<declare("rounding_offset", ShifterWidth)<<" <= "<<zg(Width+1)<<"&"<<og(wF)<<";\n";
            vhdl<<tab<<declare("rounded_raw", ShifterWidth)<<" <= shifted + rounding_offset;\n";
            vhdl<<tab<<declare("rounded", Width)<<" <= rounded_raw"<<range(ShifterWidth-1,wF+1)<<";\n"; 
         }
      }
      
      //////////////////////////////////////////////////////
      // Ok, in all normal cases we have sorted out the output,
      // so there are only overflows and underflows to deal with
      
      vhdl<<tab<<declare("max",Width)<<"<="<<og(Width)<<";\n";
      vhdl<<tab<<declare("zero",Width)<<"<="<<zg(Width)<<";\n";
      
      vhdl<<tab<<declare("res", Width)<< "<= ";
      vhdl<<tab<<tab<< "zero when (I_flags=\"00\") else\n";        // Number is zero
      vhdl<<tab<<tab<< "zero when (I_flags=\"11\") else\n";        // Number is NaN
      vhdl<<tab<<tab<< "zero when (I_sign='1') else\n" ;           // Number is negative (whether normal or infinity)
      // Must be positive and (infinite or normal)
      vhdl<<tab<<tab<< "max when (I_flags=\"10\") else\n";         // Number is positive infinity
      // Must be positive and normal
      vhdl<<tab<<tab<< "zero when (I_expnt<"<<minUsefulExpnt<<") else\n";  // Number has underflowed
      vhdl<<tab<<tab<< "max when (I_expnt>"<<maxUsefulExpnt<<") else\n";   // Number has overflowed
      vhdl<<tab<<tab<< "rounded;\n";              

      vhdl<<tab<<"O <= res;\n";
   }


   FP2FixV2::~FP2FixV2() {
   }


   void FP2FixV2::emulate(TestCase * tc)
   {
     // NOTE : There is a question about what FP2Fix should do on overflow. For underflow it seems
     // obvious to go to zero. I'm going to enforce a stronger condition on things than I think it
     // had, which is that for overflow it saturates at the largest possible value, and for underflow
     // it saturates at the most negative value. This adds a touch more hardware, but makes behaviour
     // absolutely unambiguous for all inputs except NaN. For NaN I'll arbitrarily define it to be 0 output.
      
      int Width=MSBO-LSBO+1;
      
      FPNumber  fpi(wEI, wFI, tc->getInputValue("I"));
      
      mpfr_t input, output;
      mpfr_init2(input, 1+wFI);
      fpi.getMPFR(input);
      
      mpfr_init2(output, 2+wFI+Width);
      // Lossless transform
      mpfr_mul_2si(output, input, -LSBO, MPFR_RNDN);   
      // Here is where it is rounded
      if(trunc_p){
         mpfr_floor(output, output);
      }else{
         mpfr_t tmp, tmp2;
         mpfr_init2(tmp, 2+wFI+Width);
         mpfr_init2(tmp2, 2+wFI+Width);
         
         mpfr_set(tmp, output, MPFR_RNDN);
         mpfr_round(output, output);
         
         // We're always going to round down with ties...
         // Convergent rounding be damned
         // TODO : Umm, maybe it should be convergent.
         
         mpfr_sub(tmp2, tmp, output, MPFR_RNDN);
         mpfr_abs(tmp2, tmp2, MPFR_RNDN);
         if(mpfr_cmp_d(tmp2, 0.5)==0){
            REPORT(DEBUG, "Breaking tie down");
            mpfr_floor(output, tmp);
         }
            
         mpfr_clear(tmp);
         mpfr_clear(tmp2);
      }
      
      // Might as well pull it out here
      mpz_class rawOut;
      mpfr_get_z(rawOut.get_mpz_t(), output, MPFR_RNDN);    
      
      // Then put it back in the right position
      mpfr_mul_2si(output, output, LSBO, MPFR_RNDN);      
      
      // Clamp to correct values. The extra mpz_class(.) wrappers are because of expression templates
      if(Signed){
         rawOut=std::max(mpz_class(mpz_class(-1)<<(Width-1)), std::min( mpz_class((mpz_class(1)<<(Width-1))-1), rawOut));
      }else{
         rawOut=std::max(mpz_class(0), std::min( mpz_class((mpz_class(1)<<Width)-1), rawOut));
      }
      
      REPORT(FULL, "Input : "<<input<<", Output : "<<output<<", raw="<<rawOut);
      
      tc->addExpectedOutput("O", rawOut);
      
      mpfr_clears(input, output, NULL);
  }
  
  TestCase* FP2FixV2::buildRandomTestCase(int i)
  {
     TestCase *tc;
     mpz_class a;
     tc = new TestCase(this); 
    
     mpfr_t fp;
     mpfr_init2(fp, 1+wFI);
     FPNumber fpr(wEI, wFI);
     
      // DT10 : 1/8th of the time generate a floating-point number in the input range,
      // rest of the time generate a fixed-point number with two more MSBs and four more LSBs
      if(getLargeRandom(3)==7){
         // floating-point number
         mpz_class mantissa=getLargeRandom(wFI);
         mpz_class exponent=getLargeRandom(wEI);
         mpz_class sign = Signed ? getLargeRandom(1) : 0;
         mpz_class flags = getLargeRandom(4)==0 ? 0 : 1;   // 1/16 chance of zero,  but still leaving rubbish in exponent and fraction
         mpz_class raw= (flags << (1+wFI+wEI)) +  (sign<<(wFI+wEI)) + (exponent<<wFI) + mantissa;
         fpr=FPNumber(wEI, wFI, raw);
      }else{
         // Fixed point value with MSB=MSB+1 and LSB=LSB-4
         mpz_class raw=getLargeRandom((MSBO+1)-(LSBO-4));
         
         mpfr_set_z(fp, raw.get_mpz_t(), MPFR_RNDN);
         mpfr_mul_2si(fp, fp, LSBO-4, MPFR_RNDN);
         
        fpr=FPNumber(wEI, wFI, fp);
      }  
      tc->addFPInput("I", &fpr);
     
      REPORT(FULL, "Emulate FP2FixV2 : input = "<<fp);
      /* Get correct outputs */
      emulate(tc);
      
      mpfr_clear(fp);
      
     return tc;		
   }
   
   void FP2FixV2::buildStandardTestCases(TestCaseList* tcl)
   {
      int Width=MSBO-LSBO+1;  
      int wE=wEI, wF=wFI;
	   
      TestCase *tc=new TestCase(this);
      tc->addInput("I", mpz_class(0));  
      emulate(tc);
      tcl->add(tc);
      
      tc=new TestCase(this);
      tc->addInput("I", mpz_class(1)<<(wE+wF));  
      emulate(tc);
      tcl->add(tc);
      
      tc=new TestCase(this);
      tc->addInput("I", mpz_class(3)<<(wE+wF));  
      emulate(tc);
      tcl->add(tc);
      
      tc=new TestCase(this);
      tc->addInput("I", mpz_class(2)<<(wE+wF));  
      emulate(tc);
      tcl->add(tc);
      
      mpfr_t x, delta;
      mpfr_init2(x, 1+wFI+Width);
      mpfr_init2(delta, 1+wFI+Width);
      
      // Walk over the 64 smallest numbers, some of which will underflow, some won't.
      // This isn't very smart, for some floating-point types they will always
      // be the same
      mpfr_set_d(x, pow(2.0, -LSBO-4), MPFR_RNDN);
      mpfr_set(x, delta, MPFR_RNDN);
      for(unsigned i=0;i<64;i++){
         FPNumber fpr(wE, wF, x);
         tc=new TestCase(this);
         tc->addFPInput("I", &fpr);  
         emulate(tc);
         tcl->add(tc);
         
         mpfr_add(x, x, delta, MPFR_RNDN);
      }
      
      mpfr_set_d(x, pow(2.0, MSBO), MPFR_RNDN);
      mpfr_set_d(delta, pow(2.0, LSBO-2), MPFR_RNDN);
      for(unsigned i=0;i<64;i++){
         mpfr_sub(x, x, delta, MPFR_RNDN);
      }
      for(unsigned i=0;i<128;i++){
         FPNumber fpr(wE, wF, x);
         tc=new TestCase(this);
         tc->addFPInput("I", &fpr);  
         emulate(tc);
         tcl->add(tc);
         
         mpfr_add(x, x, delta, MPFR_RNDN);
      }
      
      
      mpfr_clear(x);
      mpfr_clear(delta);
   }
}
