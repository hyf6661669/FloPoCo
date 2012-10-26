// general c++ library for manipulating streams
#include <iostream>
#include <sstream>
#include <math.h>	// for NaN



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

using namespace std;

namespace flopoco
{

extern vector<Operator *> oplist;	
	
namespace random
{
	
Operator *MakeSinglePortTable(Target *target, std::string name, int wElts, const std::vector<mpz_class> &contents, map<string, double> inputDelays = emptyDelayMap )
{
	class SinglePortTable
		: public Table
	{
	private:
		std::vector<mpz_class> m_elements;
	public:
		SinglePortTable(Target* target, std::string name, int wIn, int wOut, const std::vector<mpz_class> &elements, map<string, double> inputDelays = emptyDelayMap )
			: Table(target, wIn, wOut, /*minIn*/ 0, /*maxIn*/-1, /*logicTable*/0,  inputDelays)
			, m_elements(elements)
		{
			setName(name);
		}
		
		virtual mpz_class function(int x)
		{ return m_elements.at(x);	}
	};
	
	int wIn=(int)round(log(contents.size())/log(2.0));
	if((1<<wIn) != contents.size())
		throw std::string("MakeSinglePortTable - Number of elements must currently be a power of two.");
	
	mpz_class mm=mpz_class(1)<<wElts;
	for(unsigned i=0;i<contents.size();i++){
		if(contents[i] < 0)
			throw std::string("MakeSinglePortTable - Currently elements must be positive (need to be encoded for signed values).");
		if(contents[i] >=mm)
			throw std::string("MakeSinglePortTable - Element is not in range [0,2^wElts).");
	}
	
	return new SinglePortTable(target, name, wIn, wElts, contents, inputDelays);
}
	
TableTransform::TableTransform(Target* target, int wElts, const std::vector<mpz_class> &elements, bool addRandomSign)
	: RngTransformOperator(target)
	, m_wElts(wElts)
	, m_addRandomSign(addRandomSign)
	, m_elements(elements)
	, m_log2n((int)round(log(elements.size())/log(2.0)))
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
	
	Operator *table=MakeSinglePortTable(target, acc.str()+"_Contents", wElts, m_elements);	
	oplist.push_back(table);
		
	inPortMap(table, "X", "index");
	outPortMap(table,"Y", "elt");
	syncCycleFromSignal("elt");
	vhdl << instance(table, "elements");
	
	REPORT(DETAILED, "    Sorting out output and/or sign change.");
	
	if(m_addRandomSign){
		vhdl<<declare("res",wElts+1) << " <= "<<zeroExtend("elt",wElts+1) << " when (sign_bit='0') else ("<<zg(wElts-1)<<" - "<<zeroExtend("elt",wElts+1)<<");\n";
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


