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
#include "TableTransform.hpp"
#include "random/utils/operator_factory.hpp"

// For twos complement stuff
#include "CLTTransform.hpp"

#include "FixedPointFunctions/Function.hpp"
#include "Table.hpp"

#include "random/utils/make_table.hpp"


using namespace std;

namespace flopoco
{

extern vector<Operator *> oplist;	
	
namespace random
{
	
TableTransform::TableTransform(Target* target, int wElts, const std::vector<mpz_class> &elements, bool addRandomSign, int wF)
	: RngTransformOperator(target)
	, m_wElts(wElts)
	, m_addRandomSign(addRandomSign)
	, m_elements(elements)
	, m_log2n((int)round(log(elements.size())/log(2.0)))
	, m_wF(wF)
{
	REPORT(DETAILED, "  TableTransform() begin");	
	
	if((1<<m_log2n) != elements.size())
		throw std::string("TableTransform - Table must have binary power size.");
	
	mpz_class mm=mpz_class(1)<<wElts;
	for(unsigned i=0;i<elements.size();i++){
		if(elements[i] <0)
			throw std::string("TableTransform - Elements should not be negative (should be encoded already).");
		if(elements[i] >=mm){
			std::cerr<<"  Element "<<i<<" = "<<elements[i]<<", max="<<mm<<"\n";
			throw std::string("TableTransform - Element exceeds specified width of table.");
		}
	}
	
	std::stringstream acc;
	acc<<"TableTransform_wElts"<<wElts<<"_nElts"<<elements.size()<<"_uid"<<getNewUId();
	setName(acc.str());
	
	REPORT(DETAILED, "    Adding inputs.");
	
	addInput(uniformInputName(), uniformInputBits());
	addOutput(nonUniformOutputName(0), nonUniformOutputWidth(0));
	
	REPORT(DETAILED, "    Extracting uniform bits.");
	
	if(m_addRandomSign){
		vhdl << declare("index",m_log2n) << "<= iUniformBits"<<range(m_log2n-1,0)<<";\n";
		vhdl << declare("sign_bit") << "<= iUniformBits("<<m_log2n<<");\n";
	}else{
		vhdl << declare("index",m_log2n) << "<= iUniformBits;\n";
	}
	
	REPORT(DETAILED, "    Building single port table.");
	
	bool hardRam=elements.size() > (1<<(1+target->lutInputs()));
	Operator *table=MakeSinglePortTable(target, acc.str()+"_Contents", wElts, m_elements, hardRam);
	oplist.push_back(table);
	
	inPortMap(table, "X", "index");
	outPortMap(table,"Y", "elt");
	syncCycleFromSignal("elt");
	vhdl << instance(table, "elements");
	if(hardRam){
		useHardRAM(table);
		nextCycle(); // Force pipeline register in
	}else{
		useSoftRAM(table);
		nextCycle();	// Wwant FF for soft RAM
	}
	
	REPORT(DETAILED, "    Sorting out output and/or sign change.");
	
	if(m_addRandomSign){
		// This seems to cause weirdness in xst
		//vhdl<<declare("res",wElts+1) << " <= "<<zeroExtend("elt",wElts+1) << " when (sign_bit='0') else ("<<zg(wElts-1)<<" - "<<zeroExtend("elt",wElts+1)<<");\n";
		vhdl << declare("elt_ext", wElts+1) << "<=" << zeroExtend("elt", wElts+1)<<";\n";
		vhdl << declare("pos_mask", wElts+1) << "<= (others => not sign_bit);\n";
		vhdl << declare("neg_mask", wElts+1) << "<= (others => sign_bit);\n";
		vhdl<<declare("res",wElts+1) << " <= (pos_mask and elt_ext) - (neg_mask and elt_ext);\n";
		nextCycle();
		vhdl<<nonUniformOutputName(0)<<" <= res;\n";
	}else{
		vhdl<<nonUniformOutputName(0)<<" <= elt;\n";
	}
	
	REPORT(DETAILED, "  TableTransform() complete.");
}

TableTransform::~TableTransform()
{}
	

void TableTransform::emulate(TestCase * tc)
{
	mpz_class bits=tc->getInputValue(uniformInputName());
	
	if(m_addRandomSign){
		mpz_class sign_bit, index;
		mpz_tdiv_r_2exp(index.get_mpz_t(), bits.get_mpz_t(), m_log2n);
		sign_bit=bits>>m_log2n;
		
		mpz_class res=m_elements.at(index.get_ui());
		
		if(sign_bit!=0)
			res=-res;
		
		tc->addExpectedOutput(nonUniformOutputName(0), flopoco::random::toTwosComplement(res, m_wElts+1));
	}else{
		tc->addExpectedOutput(nonUniformOutputName(0), m_elements.at(bits.get_ui()));
	}
}
	
TestCase* TableTransform::buildRandomTestCase(int i)
{
	TestCase *tc=new TestCase(this);
	
	tc->addInput(uniformInputName(), getLargeRandom(uniformInputBits()));
	emulate(tc);
	
  	return tc;
}

typename Distribution<mpfr::mpreal>::TypePtr TableTransform::nonUniformOutputDistribution(int i, unsigned prec) const
{
	if(i!=0)
		throw std::string("Only one output distribution.");
	
	if(!m_distribution || (prec>m_distributionPrec)){
		std::vector<mpfr::mpreal> contents(m_elements.size() * (m_addRandomSign?2:0), mpfr::mpreal(0.0, prec));
		for(unsigned i=0;i<m_elements.size();i++){
			contents[i]=ldexp(mpfr::mpreal( m_elements[i].get_mpz_t(), prec), -m_wF);
			assert(contents[i].get_prec()>=(int)prec);
			if(m_addRandomSign){
				contents[2*m_elements.size()-i-1]=-contents[i];
			}
		}
		
		m_distribution= boost::make_shared<TableDistribution<mpfr::mpreal> >(&contents[0], &contents[contents.size()], m_wF);
		m_distributionPrec=prec;
	}
	assert(m_distribution->Cdf(0).get_prec()>=(int)prec);
	return m_distribution;
}


static void TableFactoryUsage(std::ostream &dst)
{
	OperatorFactory::classic_OP(dst, "table_transform", "k w func addSign", false);
	dst << "    Generates a table transform according to the given function\n";
	dst << "	      k - Number of input address bits, table will have 2^k elements.\n";
	dst << "	      w - Width of each element\n";
	dst << "        func - Function over (0,1) -> [0,1) defining entries.\n";
	dst << "        addSign - Whether to use the MSB of the input to attach a sign (resulting in output of w+1 bits, and an input of k+1 bits).\n";
	dst << "    The table will contain the entries:\n";
	dst << "      round(2^w * func[ (i+0.5)/(2^k) ] ) for i in [0..2^k)\n";
}

static Operator *TableFactoryParser(Target *target ,const std::vector<std::string> &args,int &consumed)
{
	unsigned nargs = 4;
	if (args.size()<nargs)
		throw std::string("TableFactory - Not enough arguments, check usage.");
	consumed += nargs;
	
	int k = atoi(args[0].c_str());
	int w = atoi(args[1].c_str());
	std::string funcStr = args[2];
	bool addSign=atoi(args[3].c_str())!=0;

	if(k<1)
		throw std::string("TableFactory - k must be a positive integer.");
	if(w<1)
		throw std::string("TableFactory - w must be a positive integer.");
	if(w>=48)
		throw std::string("TableFactory - w must be less than 48 (currently table is built in double-precision).");
	
	if(DETAILED<=::flopoco::verbose)
		std::cerr<<"  Parsing sollya string \""<<funcStr<<"\" ... ";
	Function func(funcStr);
	if(DETAILED<=::flopoco::verbose)
		std::cerr<<"done\n";
	
	std::vector<mpz_class> contents(1<<k);
	for(int i=0;i<(1<<k);i++){
		double x=(i+0.5)/(1<<k);

		double y=func.eval(x);
		double ry=round(ldexp(y,w));
		
		if(DEBUG<=::flopoco::verbose)
			std::cerr<<"    "<<i<<", x="<<x<<", y="<<y<<", ry="<<ry<<"\n";
		
		if((ry<0.0) || (ry>=ldexp(1.0,w))){
			std::stringstream acc;
			acc<<"TableFactoryParser : For index "<<i<<", the value "<<y<<"=f("<<x<<") leads to out of range value "<<ry;
			throw std::string(acc.str());
		}
		contents[i]=ry;
	}
	
	return new TableTransform(target, w, contents, addSign);
}

void TableTransform::registerFactory()
{
	DefaultOperatorFactory::Register(
		"table_transform",
		"operator;rng_transform",
		flopoco::random::TableFactoryUsage,
		flopoco::random::TableFactoryParser,
		DefaultOperatorFactory::ParameterList(
			DefaultOperatorFactory::Parameters("8", "8", "sin(x)", "0"),
			DefaultOperatorFactory::Parameters("8", "8", "sin(x)", "1")
		)
	);
}

};
};


