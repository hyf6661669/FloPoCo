#ifndef __lut_sr_urng_HPP
#define __lut_sr_urng_HPP
#include <vector>
#include <sstream>
#include <set>

#include "Operator.hpp"
#include "../transforms/RngTransformOperator.hpp"

/* This file contains a lot of useful functions to manipulate vhdl */
#include "utils.hpp"

#include "../transforms/RngTransformOperator.hpp"

namespace flopoco{
namespace random{


// new operator class declaration
class MWCRng : public Operator {
  public:
    static string operatorInfo;
    
    struct MWCRngParams{
      unsigned w;  //  w The underlying width of an underlying digit
      mpz_class M;  // M The multiplier, with requirement that 0 < M < 2^w
      bool useCSAInternally;    // Use CSA for the internals, rather than standard multipliers
      bool useCSAExternally;    // Use CSA for the outputs, resulting in oRngCarry and oRngSum. Requires useCSAInternally
      
      MWCRngParams()
        : w(0)
        , M(0)
        , useCSAInternally(false)
        , useCSAExternally(false)
      {}
    };
      
    
  private:
    MWCRngParams m_params;
    unsigned m_w;
    mpz_class m_M;
    unsigned m_wM;
    unsigned m_wO;
    bool m_useCSAInternally;
    bool m_useCSAExternally;
  
    mpz_class m_state;
  public:

   /*! \note Output value has width w+log2ceil(M)*/
    MWCRng(Target* target, const MWCRngParams &p);

    // destructor
    ~MWCRng();


    // Below all the functions needed to test the operator
    /* the emulate function is used to simulate in software the operator
      in order to compare this result with those outputed by the vhdl opertator */
    void emulate(TestCase * tc);

    /* function used to create Standard testCase defined by the developper */
  //  void buildStandardTestCases(TestCaseList* tcl);


	void buildStandardTestCases(TestCaseList* tcl);


  /** Take the given transform, and supply it with uniform bits. The inputs of the
    * resulting operator will be the LUT-SR inputs, while the outputs will be those
    * of the transform operator
   */
  //static Operator *DriveTransform(std::string name, RngTransformOperator *base);
};

}; // random
}; // flopoco
#endif
