#ifndef __mwc_transform_HPP
#define __mwc_transform_HPP
#include <vector>
#include <sstream>
#include <set>

#include "Operator.hpp"
#include "utils.hpp"

#include "../transforms/RngTransformOperator.hpp"


namespace flopoco{
namespace random{


// new operator class declaration
class MWCTransform : public Operator {
  public:
    static string operatorInfo;
    
    struct MWCTransformParams{
      unsigned w;  //  w The underlying width of an underlying digit
      mpz_class M;  // M The multiplier, with requirement that 0 < M < 2^w
      bool useCSA;  // Both input and output will be in CSA format
      
      MWCTransformParams()
        : w(0)
        , M(0)
        , useCSA(false)
      {}
    };
      
    
  private:
    MWCTransformParams m_p;
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
    MWCTransform(Target* target, const MWCTransformParams &p);

    // destructor
    ~MWCTransform();

    void emulate(TestCase * tc);
    void buildStandardTestCases(TestCaseList* tcl);
};

}; // random
}; // flopoco
#endif
