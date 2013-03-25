// general c++ library for manipulating streams
#include <iostream>
#include <sstream>
#include <math.h>	// for NaN
#include <fstream>

/* header of libraries to manipulate multiprecision numbers
  There will be used in the emulate function to manipulate arbitraly large
  entries */
#include "gmp.h"
//#include "mpfr.h"
//#include "FPNumber.hpp"

// include the header of the Operator
#include "mwc_rng.hpp"
#include "../utils/chain_operator.hpp"
#include "../utils/operator_factory.hpp"

#include <boost/math/common_factor.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;

namespace flopoco
{
	extern std::vector<Operator *> oplist;
	
namespace random
{


// personalized parameter
string MWCRng::operatorInfo = "MWCRng <NotAValue/FixMe>";


MWCRng::MWCRng(Target* target, const MWCRngParams &p)
: Operator(target)
, m_params(p)
, m_w(p.w)
, m_M(p.M)
, m_wM(mpz_sizeinbase(m_M.get_mpz_t(), 2))
, m_wO(m_w+m_wM)
{
  if(p.useCSAInternally)
    throw std::string("MWCRng - CSA mode not implemented yet.");
  if(p.useCSAExternally)
    throw std::string("MWCRng - CSA mode not implemented yet.");

  setCopyrightString("David Thomas 2013");
  
  std::stringstream name;
  name<<"MWCRng_w"<<m_w<<"_M"<<m_M;
  setName(name.str());

  // We have recurrences of one cycle, anything more is an error
  setHasDelay1Feedbacks();

  addInput("m");			//mode-> m=1 load; m=0 RNG
  addInput("Sin");		//serial load input in load mode
  addOutput("RNG", m_wO, 1, true);
  addOutput("Sout");

  declare("state_n",m_wO,false, Signal::registeredWithZeroInitialiser);
  
  nextCycle();
  declare("state", m_wO);
  vhdl<<tab<<"state <= state_n;\n\n";

  setCycleFromSignal("m");
  
  vhdl<<declare("multiplier",m_wO)<<" <= state"<<range(m_w-1,0)<<" * \""<<unsignedBinary(m_M, m_wM)<<"\";\n";
  addAttribute("mult_style",  "string",  "multiplier:signal", "lut" );
  
  vhdl<<tab<<"state_n <= state"<<range(m_wO-2,0)<<" & Sin when m='1' else\n";
  vhdl<<tab<<tab<<"multiplier + ("<<zg(m_w)<<" & state"<<range(m_wO-1,m_w)<<");\n\n";
  
  vhdl << tab << "RNG <= state;\n";
  vhdl << tab << "Sout <= state("<<m_wO-1<<");" << endl;
}

/*
std::pair<std::string,std::string> MWCRng::MakeCompressor(std::vector<std::string> srcs)
{
  if(srcs.size()==0){
    return "'0'";
  }else if(srcs.size()==1){
    return srcs[0];
  }else if(srcs.size()==2){
    return srcs[1];
  }
}

void MWCRng::MakeCSAMult(std::string dstC, std::string dstS, std::string srcC, std::string srcS, unsigned i, const std::vector<unsigned> &shifts)
{
  
}
*/
//============================================================================================================================

MWCRng::~MWCRng(){	
}

//============================================================================================================================
void MWCRng::emulate(TestCase * tc) {
	//std::cerr<<"LutSrRng::emulate\n";
  
  mpz_class smode = tc->getInputValue("m");
  mpz_class ssin= tc->getInputValue("Sin");

  if(smode==1){
    m_state= (m_state<<1) + ssin;
    mpz_fdiv_r_2exp (m_state.get_mpz_t(), m_state.get_mpz_t(), m_wO);
  }else{
    mpz_class carry=m_state>>m_w;
    mpz_fdiv_r_2exp (m_state.get_mpz_t(), m_state.get_mpz_t(), m_w);
    m_state=(m_state*m_M) + carry;
  }
  
  mpz_class s_out=m_state>>(m_wO-1);

  tc->addExpectedOutput("RNG",m_state);
  tc->addExpectedOutput("Sout",s_out);
}

//============================================================================================================================
void MWCRng::buildStandardTestCases(TestCaseList* tcl)
{
    TestCase *tc;
  
    // Loading in the state
	for(unsigned i=0;i<m_wO;i++)
	{
		
      tc = new TestCase(this);
      tc->setSetupCycle(true);
      
      tc->addInput("Sin", getLargeRandom(1));
      tc->addInput("m", mpz_class(1));

      emulate(tc);
      tcl->add(tc);
    }

  // Running the RNG
  for(unsigned i=0;i<10*m_wO;i++)
	{
      tc = new TestCase(this);
      tc->addInput("Sin", mpz_class(0));
      tc->addInput("m", mpz_class(0));

      emulate(tc);
      tcl->add(tc);
	}
	
    // Check readback
    for(unsigned i=0;i<m_wO;i++)
	{

      tc = new TestCase(this); 
      tc->addInput("Sin", mpz_class(1));
      tc->addInput("m", mpz_class(1));
      emulate(tc);
      tcl->add(tc);
	}
}

static bool IsPrime(const mpz_class &x)
{
  return mpz_probab_prime_p(x.get_mpz_t(),32)!=0;
}

static bool IsSafePrime(const mpz_class &x)
{
  return IsPrime(x) && IsPrime((x-1)/2);
}

static bool IsSafeMultiplier(unsigned w, const mpz_class &M)
{
  mpz_class base=mpz_class(1)<<w;
  mpz_class x=base*M-1;
  return IsSafePrime(x);
}


static void MWCFactoryUsage(std::ostream &dst)
{
	OperatorFactory::classic_OP(dst, "MWCRng", "w m", false);
	dst << "    Generates a multiply-with-carry RNG\n";
	dst << "	      w - Width of base digit.\n";
	dst << "	      M - multiplier, must have 0<M<2^w. If M==0, then it will try to auto-select\n";
}

static Operator *MWCFactoryParser(Target *target ,const std::vector<std::string> &args,int &consumed)
{
	unsigned nargs = 2;
	if (args.size()<nargs)
		throw std::string("MWCFactory - Not enough arguments, check usage.");
	consumed += nargs;
    
    MWCRng::MWCRngParams params;
	
	params.w = atoi(args[0].c_str());
    params.M = mpz_class(args[1]);
    
	if(params.w<1)
		throw std::string("MWCFactory - w must be a positive integer.");
	if((params.M < 0) || (params.M >= (mpz_class(1)<<params.w)))
		throw std::string("MWCFactory - Must have 0<M<2^2");
    if(params.M==0){
      if(::flopoco::verbose >=DETAILED)
        std::cerr<<"Selecting a multiplier.\n";
      
      params.M=(mpz_class(1)<<params.w)-1;
      while(params.M>0){
        std::cerr<<params.M<<"\n";
        if(IsSafeMultiplier(params.w, params.M))
          break;
        params.M-=2;    // Move down over odd numbers
      }
      if(params.M<=0)
        throw std::string("MWCFactory - Couldn't find a multiplier.");
    }
	
	return new MWCRng(target, params);
}

void MWCRng_registerFactory()
{
	DefaultOperatorFactory::Register(
		"MWCRng",
		"operator;urng",
		flopoco::random::MWCFactoryUsage,
		flopoco::random::MWCFactoryParser,
		DefaultOperatorFactory::ParameterList(
			DefaultOperatorFactory::Parameters("32", "0")
		)
	);
}

}; // random
}; // flopoco
