#ifndef __mwc_rng_HPP
#define __mwc_rng_HPP
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
    MWCRngParams m_p;
    unsigned m_w;
    mpz_class m_M;
    unsigned m_wM;
    unsigned m_wO;
  
    mpz_class m_state, m_state2;

	void BuildCSA();
    std::string Invert(std::string prefix, std::string x, bool dryRun);
    std::pair<std::string,std::string> MakeHalfAdder(std::string prefix, std::string A, std::string B, bool dryRun);
    std::pair<std::string,std::string> MakeFullAdder(std::string prefix, std::string A, std::string B, std::string C, bool dryRun);
  
    std::vector<int> DecomposeM();  
  
    std::vector<std::vector<std::string> > OptimiseOnesCSA(const std::vector<std::vector<std::string> > &curr);
    std::vector<std::vector<std::string> > NegateCSA(std::string prefix, std::vector<std::vector<std::string> > curr, bool dryRun);  
    std::vector<std::vector<std::string> > ProcessCSALevel(std::string prefix, const std::vector<std::vector<std::string> > &curr, bool dryRun);  
	std::vector<std::vector<std::string> > MakeMWCFromCSA(const std::vector<std::vector<std::string> > &acc, bool dryRun);
  
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
