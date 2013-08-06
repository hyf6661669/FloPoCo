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
, m_p(p)
, m_w(p.w)
, m_M(p.M)
, m_wM(mpz_sizeinbase(m_M.get_mpz_t(), 2))
, m_wO(m_w+m_wM)
{
  setCopyrightString("David Thomas 2013");
  
  std::stringstream name;
  name<<"MWCRng_w"<<m_w<<"_M"<<m_M;
  setName(name.str());

  // We have recurrences of one cycle, anything more is an error
  setHasDelay1Feedbacks();

  if(p.useCSAInternally){
    BuildCSA();
  }else{
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
}

void MWCRng::BuildCSA()
{
  addInput("m");			//mode-> m=1 load; m=0 RNG
  addInput("Sin");		//serial load input in load mode
  if(m_p.useCSAExternally){
    addOutput("RNG", 2*m_wO, 1, true);
  }else{
    addOutput("RNG", m_wO, 1, true);
  }
  //addOutput("Sout");   // Too much hassle supporting proper shift out.

  declare("state_n_c",m_wO-1,false, Signal::registeredWithZeroInitialiser);
  declare("state_n_s",m_wO,false, Signal::registeredWithZeroInitialiser); // This has one less bit, because we don't care about overflow
  
  nextCycle();
  declare("state_c", m_wO-1);
  declare("state_s", m_wO);
  vhdl<<tab<<"state_c <= state_n_c;\n\n";
  vhdl<<tab<<"state_s <= state_n_s;\n\n";
  
    setCycleFromSignal("m");    
  
  // This is the current state of the accumulator, with each bit represented by one or two bits
  std::vector<std::vector<std::string> > state(m_wO);
  for(unsigned i=0;i<m_wO;i++){
    state[i].push_back(join("state_s(",i,")"));
    if(i>0)
      state[i].push_back(join("state_c(",i-1,")"));
  }
  
  /*
  std::vector<std::vector<std::string> > state_i(state), state_o;
  while(1){
    state_o=MakeMWCFromCSA(state_i, true);
    if(state_o==state_i)
      break;
    std::cerr<<"Repeat\n";
    state_i=state_o;
  }
  
  for(unsigned i=0;i<m_wO;i++){
    std::cerr<<i<<" : "<<state_o[i].size()<<"\n";
  }
  exit(1);
  */
  
  // Build up the recurrent multiply and add
  std::vector<std::vector<std::string> > state_n=MakeMWCFromCSA(state, false);
  for(unsigned i=0;i<m_wO;i++){
    assert(state_n[i].size()>=1);
    assert(state_n[i].size()<=2);
    vhdl<<tab<<"state_n_s("<<i<<") <= "<<state_n[i][0]<<" when m='0' else ";
    if(i==0){
      vhdl<<"Sin;\n";
    }else{
      vhdl<<"state_s("<<i-1<<");\n";
    }
    if((i>0) && (i<m_wO-1) && (state_n[i-1].size()==2)){
      vhdl<<tab<<"state_n_c("<<i<<") <= "<<state_n[i][1]<<" when m='0' else '0';\n";
    }else if(i==0){
      vhdl<<tab<<"state_n_c("<<i<<") <= '0';\n";
    }
  }
  
  
  if(m_p.useCSAExternally){
    vhdl<<tab<<"RNG <= state_c & state_s;\n";
  }else{
    // Do the final CSA to binary conversion, still requires a wO bit ripple carry add
    vhdl<<tab<<"RNG <= state_s + (state_c&'0');\n";
    nextCycle();
  }
}


std::pair<std::string,std::string> MWCRng::MakeHalfAdder(std::string prefix, std::string A, std::string B, bool dryRun)
{
  std::string S=join(prefix,"_S");
  std::string C=join(prefix,"_C");
  
  if(!dryRun){
    declare(S);
    declare(C);
  
    vhdl<<tab<<S<<" <= ("<<A<<" xor "<<B<<");\n";
    vhdl<<tab<<C<<" <= ("<<A<<" and "<<B<<");\n";
  }
  
  return std::make_pair(C,S);
}

std::pair<std::string,std::string> MWCRng::MakeFullAdder(std::string prefix, std::string A, std::string B, std::string C, bool dryRun)
{
  std::string dS=join(prefix,"_S");
  std::string dC=join(prefix,"_C");
  
  if(!dryRun){
    declare(dS);
    declare(dC);
  
    vhdl<<tab<<dS<<" <= ("<<A<<" xor "<<B<<" xor "<<C<<");\n";
    vhdl<<tab<<dC<<" <= ("<<A<<" and "<<B<<")  or ("<<C<<" and ("<<A<<" xor "<<B<<"));\n";
  }
  
  return std::make_pair(dC,dS);
}

std::string MWCRng::Invert(std::string prefix, std::string x, bool dryRun)
{
  if(x=="'0'")
    return "'1'";
  if(x=="'1'")
    return "'0'";
  std::string res=prefix+"_inv";
  if(!dryRun){
    declare(res);
    vhdl<<tab<<res<<" <= not "<<x<<";\n";
  }
  return res;
}

/*! Given a CSA number, take out any zeros, and try to redistribute ones such that there is at most a single one in any column */
std::vector<std::vector<std::string> > MWCRng::OptimiseOnesCSA(const std::vector<std::vector<std::string> > &curr)
{
	std::vector<std::vector<std::string> > res(curr.size());
	
	// Remove any zeros, and count all the ones
	std::vector<int> ones(curr.size());
	for(unsigned i=0;i<res.size();i++){
		for(unsigned j=0;j<curr[i].size();j++){
			if(curr[i][j]=="'0'"){
				// Do nothing
			}else if(curr[i][j]=="'1'"){
				ones[i]++;
			}else{
				res[i].push_back(curr[i][j]);
			}
		}
	}
	
	// Propagate ones up the chain, so there is at most a single one at any level
	for(unsigned i=0;i<res.size();i++){
		while(ones[i]>1){
			if(i+1<res.size()){
				ones[i+1]++;
			}
			ones[i]-=2;
		}
		assert(ones[i]<=1);
	}
	
	// Distribute signals to minimise height
	for(unsigned i=0;i<res.size();i++){
		// Given a one, we can choose to push a single "real" signal up to the next level,
		// or keep it here,  as   sum(1+X)==1 and carry(1+X)==X
		if(ones[i]){
			if(i==res.size()-1){
				// If this is the MSB, always promote, as that loses us a signal
				res[i].pop_back();
				res[i].push_back("'1'");
			}else{
				// If the next level has lower height, then move it up
				if(res[i+1].size() < res[i].size()){
					res[i+1].push_back(res[i].back());
					res[i].pop_back();
					res[i].push_back("'1'");
				}
			}
		}
	}
	
	return res;
}

std::vector<std::vector<std::string> > MWCRng::NegateCSA(std::string prefix, std::vector<std::vector<std::string> > curr, bool dryRun)
{
	curr=OptimiseOnesCSA(curr);
	
  // Work out the maximum height of bits
  unsigned count=0;
  for(unsigned i=0;i<curr.size();i++){
    count=std::max(count, curr[i].size());
  }
  
  // Invert everything
  std::vector<std::vector<std::string> > res(curr.size());
  for(unsigned i=0;i<curr.size();i++){
    for(unsigned j=0;j<curr[i].size();j++){
      res[i].push_back(Invert(str(boost::format("%1%_%2%_%3%")%prefix%i%j), curr[i][j], dryRun));
    }
    for(unsigned j=curr[i].size();j<count;j++){
      res[i].push_back("'1'");
    }
  }
  
  // And add one for each bit of height at the MSB
  unsigned pos=0;
  while(count>0){
    if(count&1){
      if(pos<res.size())
		res[pos].push_back("'1'");
    }
    pos++;
    count=count>>1;
  }
  
  return OptimiseOnesCSA(res);
}

std::vector<std::vector<std::string> > MWCRng::ProcessCSALevel(std::string prefix, const std::vector<std::vector<std::string> > &curr, bool dryRun)
{
  std::vector<std::vector<std::string> > next;
  
  bool newRound=false;
  
  for(unsigned i=0;i<curr.size();i++){
    std::vector<std::string> todo=curr[i];
    std::vector<std::string> sum, carry;
    
    while(todo.size()>0){
      if(todo.size()==1){
        sum.push_back(todo[0]);
        todo.clear();
      }else if(todo.size()==2){
        std::string name=str(boost::format("%1%_%2%_%3%")%prefix%i%sum.size());
        std::pair<std::string,std::string> ha=MakeHalfAdder(name, todo[0],todo[1],dryRun);
        todo.clear();
        carry.push_back(ha.first);
        sum.push_back(ha.second);
      }else{
        newRound=true;
        std::string name=str(boost::format("%1%_%2%_%3%")%prefix%i%sum.size());
        std::pair<std::string,std::string> fa=MakeFullAdder(name, todo[todo.size()-3],todo[todo.size()-2],todo[todo.size()-1],dryRun);
        todo.pop_back(); todo.pop_back(); todo.pop_back();
        carry.push_back(fa.first);
        sum.push_back(fa.second);
      }
    }
    
    if(next.size()==i){
      next.push_back(sum);
    }else{
      next[i].insert(next[i].begin(), sum.begin(), sum.end());
    }
    if(carry.size()>0){
      next.push_back(carry);
    }
  }
  
  if(!newRound)
    next.clear();
  
  return next;
}

std::vector<int> MWCRng::DecomposeM()
{
	// TODO : Booth encode stuff...
	std::vector<int> res;
	for(unsigned i=0;i<m_wM;i++){
		if(mpz_tstbit(m_M.get_mpz_t(), i)){
			res.push_back(1);
		}else{
			res.push_back(0);
		}
	}
	return res;
}

std::vector<std::vector<std::string> > MWCRng::MakeMWCFromCSA(
  const std::vector<std::vector<std::string> > &acc, // This is the current state of the accumulator
  bool dryRun   // Don't create any signals
)
{
	REPORT(DEBUG, "MakeMWCFromCSA - Begin");

	std::vector<std::vector<std::string> > parts(m_wO); // We don't care about overflows

	// Add the +C part
	for(unsigned i=0;i<m_wM;i++){
		parts[i]=acc.at(i+m_w);
	}
	// Add all the *M parts
	std::vector<int> M_parts=DecomposeM();
	  
	for(unsigned i=0;i<M_parts.size();i++){
		std::vector<std::vector<std::string> > x(m_wO);
		if(M_parts[i]==1){
			REPORT(DEBUG, "    Adding +x<<"<<i);
			for(unsigned j=0;j<m_w;j++){
				assert(i+j<m_wO);
				x[i+j]=acc[j];
			}
		}else if(M_parts[i]==-1){
			REPORT(DEBUG, "    Adding -x<<"<<i);
			for(unsigned j=0;j<m_w;j++){
				assert(i+j<m_wO);
				x[i+j]=acc[j];
			}
			x=NegateCSA(str(boost::format("init_neg_%1%")%i), x, dryRun);
		}else if(M_parts[i]==0){
			// Do nothing
		}else{
			throw std::string("Unexpected weight in M_parts.");
		}
		for(unsigned j=0;j<m_w;j++){
			assert(i+j<m_wO);
			parts[i+j].insert(parts[i+j].end(), acc[j].begin(), acc[j].end());
		}
	}
  
	std::vector<std::vector<std::string> > curr=parts;
	unsigned stage=0;
	while(1){
		unsigned maxHeight=0;
		for(unsigned i=0;i<curr.size();i++){
			maxHeight=std::max(maxHeight, curr[i].size());
		}
		REPORT(DEBUG, "    Stage="<<stage<<", CurrentHeight="<<maxHeight);

		std::vector<std::vector<std::string> > next=ProcessCSALevel(join("csa_",stage), curr, dryRun);
		if(next.size()==0)
			break;
		curr=next;
		++stage;
	}

	REPORT(DEBUG, "MakeMWCFromCSA - Done");

	return curr;
}

//============================================================================================================================

MWCRng::~MWCRng(){	
}

//============================================================================================================================
void MWCRng::emulate(TestCase * tc) {
	//std::cerr<<"LutSrRng::emulate\n";
	
	m_state2=m_state;
  
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

  if(!m_p.useCSAInternally){
		tc->addExpectedOutput("RNG",m_state);
		tc->addExpectedOutput("Sout",s_out);
  }else{
	  if(!m_p.useCSAExternally){
		  tc->addExpectedOutput("RNG",m_state);
	  }
  }
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
	OperatorFactory::classic_OP(dst, "MWCRng", "[-csa] w m", false);
	dst << "    Generates a multiply-with-carry RNG\n";
	dst << "	      w - Width of base digit.\n";
	dst << "	      M - multiplier, must have 0<M<2^w. If M==0, then it will try to auto-select\n";
    dst <<"          -csa - Use CSA internally.\n";
}

static Operator *MWCFactoryParser(Target *target ,const std::vector<std::string> &args,int &consumed)
{
    MWCRng::MWCRngParams params;
  
    unsigned src=0;
    if(args.size()-src > 0){
      if(args[src]=="-csa"){
        params.useCSAInternally=true;
        ++src;
      }
    }
  
	unsigned nargs = 2;
	if (args.size()-src<nargs)
		throw std::string("MWCFactory - Not enough arguments, check usage.");
	consumed = nargs+src;
	
	params.w = atoi(args[src].c_str());
    params.M = mpz_class(args[src+1]);
    
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
