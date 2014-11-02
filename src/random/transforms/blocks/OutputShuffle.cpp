// general c++ library for manipulating streams
#include <iostream>
#include <sstream>
#include <math.h>	// for NaN

#include <random/utils/mpreal/boost_math_mpreal.hpp>



/* header of libraries to manipulate multiprecision numbers
  There will be used in the emulate function to manipulate arbitraly large
  entries */
#include "gmp.h"
#include "mpfr.h"
#include "FPNumber.hpp"

// include the header of the Operator
#include "OutputShuffle.hpp"
#include "random/utils/operator_factory.hpp"


using namespace std;

namespace flopoco
{

extern vector<Operator *> oplist;	
	
namespace random
{

OutputShuffle::OutputShuffle(Target* target, unsigned k, unsigned w)
	: Operator(target)
	, m_k(k)
	, m_w(w)
{
	REPORT(DETAILED, "  OutputShuffle() begin");	
	
	setSequential();
	
	std::stringstream acc;
	acc<<"OutputShuffle_"<<k<<"_"<<w<<"_uid"<<getNewUId();
	setName(acc.str());
	
	REPORT(DETAILED, "    Adding inputs.");
	
	addInput("iData", w);
	addInput("iIncrement", k-1);
	addOutput("oData", w);
	
	nextCycle();
	
	REPORT(DETAILED, "  OutputShuffle() complete.");
}

OutputShuffle::~OutputShuffle()
{}
	
void OutputShuffle::outputVHDL(std::ostream& o, std::string name)
{
	if(!isSequential()){
		throw std::string("OutputShuffle must be sequential.");
	}
	
	licence(o);
		o << "library ieee; " << endl;
		o << "use ieee.std_logic_1164.all;" << endl;
		o << "use ieee.numeric_std.all;" << endl;
		o << "library work;" << endl;
		outputVHDLEntity(o);
	
	newArchitecture(o,name);
	
	o << tab << "type storage_t is array(0 to "<<(1<<m_k)-1<<") of  std_logic_vector("<< m_w-1<<" downto 0);" << endl;
	o << tab << "signal theRam : storage_t := (others=>(others=>'0'));"<<std::endl;
	o << tab << "signal offset : unsigned("<<m_k-1<<" downto 0) := (others=>'0');"<<std::endl;
	o << tab << "signal prev_offset : unsigned("<<m_k-1<<" downto 0) := (0=>'1', others=>'0');"<<std::endl;
	o << tab << "signal ram_buffer : std_logic_vector("<<m_w-1<<" downto 0) := (others=>'0');"<<std::endl;
	
	outputVHDLSignalDeclarations(o);
	beginArchitecture(o);
	
	o << "	process(clk)" << endl;
	o << tab << "begin" << endl;
	o << tab << "if(rising_edge(clk)) then" << endl;
	o << tab << tab << "offset <= offset + unsigned('0'&iIncrement)+1;"<<std::endl;
	o << tab << tab << "prev_offset <= offset;"<<std::endl;
	o << tab << tab << "ram_buffer <= theRam(TO_INTEGER(offset));"<<std::endl;
	o << tab << tab << "theRam(TO_INTEGER(prev_offset)) <= iData;"<<std::endl;
	o << tab << "end if;" << endl;
	o << tab << "end process;" << endl;
	o << tab << "oData <= ram_buffer;"<<endl;
	endArchitecture(o);
}
	

void OutputShuffle::emulate(TestCase * tc)
{
	throw std::string("OutputShuffle::emulate - Not implemented.");
}
	
TestCase* OutputShuffle::buildRandomTestCase(int i)
{
	throw std::string("OutputShuffle::buildRandomTestCase - Not implemented.");
}


static void OutputShuffleFactoryUsage(std::ostream &dst)
{
	OperatorFactory::classic_OP(dst, "OutputShuffle", "k w", false);
	dst << "    Applies a shuffle according toa random index\n";
	dst << "	      k - Number of table address bits, iIndex will have k-1 bits.\n";
	dst << "	      w - Width of each element\n";
}

static Operator *OutputShuffleFactoryParser(Target *target ,const std::vector<std::string> &args,int &consumed)
{
	unsigned nargs = 2;
	if (args.size()<nargs)
		throw std::string("TableFactory - Not enough arguments, check usage.");
	consumed += nargs;
	
	int k = atoi(args[0].c_str());
	int w = atoi(args[1].c_str());
	
	if(::flopoco::verbose>DETAILED)
		std::cerr<<"OutputShuffleParser : k="<<k<<", w="<<w<<"\n";

	if(k<1)
		throw std::string("TableFactory - k must be at least 3.");
	if(w<1)
		throw std::string("TableFactory - w must be a positive integer.");
	
	return new OutputShuffle(target, k ,w);
}

void OutputShuffle::registerFactory()
{
	DefaultOperatorFactory::Register(
		"OutputShuffle",
		"operator",
		flopoco::random::OutputShuffleFactoryUsage,
		flopoco::random::OutputShuffleFactoryParser,
		DefaultOperatorFactory::ParameterList(
			DefaultOperatorFactory::Parameters("4", "8"),
			DefaultOperatorFactory::Parameters("10","16")
		)
	);
}

};
};


