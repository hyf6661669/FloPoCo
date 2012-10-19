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
	if((1<<m_log2n) + elements.size())
		throw std::string("TableTransform - Table must have binary power size.");
	
	mpz_class mm=mpz_class(1)<<m_log2n;
	for(unsigned i=0;i<elements.size();i++){
		if(elements[i] <0)
			throw std::string("TableTransform - Elements should not be negative (should be encoded already).");
		if(elements[i] >=mm)
			throw std::string("TableTransform - Element exceeds specified width of table.");
	}
	
	std::stringstream acc;
	acc<<"TableTransform_wElts"<<wElts<<"_nElts"<<elements.size()<<"_uid"<<getNewUId();
	setName(acc.str());
	
	addInput(uniformInputName(), uniformInputBits());
	addOutput(nonUniformOutputName(0), nonUniformOutputWidth(0));
	
	if(m_addRandomSign){
		vhdl << declare("index",m_log2n) << "<= iUniformBits"<<range(m_log2n-1,0)<<";\n";
		vhdl << declare("sign_bit",1) << "<= iUniformBits("<<m_log2n<<");\n";
	}else{
		vhdl << declare("index",m_log2n) << "<= iUniformBits;\n";
	}
	
	Operator *table=MakeSinglePortTable(target, acc.str()+"_Contents", wElts, m_log2n, m_elements);	
	inPortMap(table, "X", "index");
	outPortMap(table,"Y", "elt");
	syncCycleFromSignal("elt");
	instance(table);
	
	if(m_addRandomSign){
		vhdl<<declare("res",baseWidth+1) << " <= "<<zeroExtend("left",baseWidth+1) << " when (sign_bit='1') else (- "<<zeroExtend("right",baseWidth+1)<<");\n";
		nextCycle();
		vhdl<<nonUniformOutputName(0)<<" <= res;\n";
	}else{
		vhdl<<nonUniformOutputName(0)<<" <= elt;\n";
	}
}

TableTransform::~TableTransform()
{}
	

void CLTTransform::emulate(TestCase * tc)
{
	mpz_class bits=tc->getInputValue(uniformInputName());
	
	if(m_addRandomSign){
		mpz_class sign_bit, index;
		mpz_tdiv_r_2exp(index.get_mpz_t(), bits.get_mpz_t(), m_log2n);
		sign_bit=index>>m_log2n;
		
		mpz_class res=m_elements.at(index.get_uint());
		
		if(sign_bit!=0)
			res=-res;
		
		return toTwosComplement(res, m_baseWidth+1);
	}else{
		return m_elements.at(bits.get_uint());
	}
}
	
TestCase* CLTTransform::buildRandomTestCase(int i)
{
	TestCase *tc=new TestCase(this);
	
	tc->addInput(uniformInputName(), getLargeRandom(uniformInputBits()));
	emulate(tc);
	
  	return tc;
}

};
};
