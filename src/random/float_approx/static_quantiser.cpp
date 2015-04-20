#include "static_quantiser.hpp"
#include <assert.h>

#include "Table.hpp"

#include "random/utils/comparable_float_type.hpp"
#include "random/utils/chain_operator.hpp"
#include "random/utils/operator_factory.hpp"

#include "mpreal.h"
#include <sstream>

namespace flopoco
{
namespace random
{

/*! Given an input value x, and a set of n+1 bucket boundaries
	find 0<=oY<n where boundaries[oY] <=iX < boundaries[oY+1]
*/
StaticQuantiser::StaticQuantiser(Target *target, int wValue, const std::vector<mpz_class> &_boundaries, map<string, double> inputDelays)
	: Operator(target, inputDelays)
	, m_wValue(wValue)
	, boundaries(_boundaries)
{
	ostringstream name;
	name << "StaticQuantiser_w"<<wValue<<"_uid"<<getNewUId();
	setName(name.str());
	setCopyrightString("Imperial College 2012");
	
	assert(boundaries.size()>=2);
	int nBuckets=boundaries.size()-1;
	
	int wIndex=(int)ceil(log(nBuckets)/log(2.0));
	
	// Force a binary power by adding extra empty buckets.
	// TODO : This wastes RAM, potentially quite a lot. However, it makes indexing easier
	while( (1<<wIndex) != nBuckets){
		boundaries.push_back(boundaries.back());
		nBuckets++;
	}
	
	REPORT(DETAILED, "StaticQuantiser, orig(n)="<<_boundaries.size()-1<<", curr(n)="<<boundaries.size()-1<<"\n");
	
	addInput ("iX" , wValue);
	addOutput("oY", wIndex);
	
	setCriticalPath( getMaxInputDelays(inputDelays) );
	for(int i=1;i<=wIndex;i++){
		AddLevel(i);
	}
	
	vhdl << "oY <= "<<join("idx_",wIndex)<<";\n";
	outDelayMap["oY"] = getCriticalPath();
}

StaticQuantiser::~StaticQuantiser()
{}
	
void StaticQuantiser::emulate(TestCase * tc)
{
	
	mpz_class iX=tc->getInputValue("iX");
	
	unsigned i=1;
	while(true){
		if(boundaries[i]>iX){
			break;
		}
		i++;
		if(i==boundaries.size()-1)
			break;
	}
	
	mpz_class oY(i-1);
	
	//std::cerr<<std::hex<<" StaticQuantiser::emulate("<<iX<<") -> "<<oY<<"\n"<<std::dec;
	
	tc->addExpectedOutput("oY", oY);
}

void StaticQuantiser::AddLevel(int level)
{
	class StaticQuantiserTable
	: public Table
	{
	private:	
		std::vector<mpz_class> m_values;
	public:
		StaticQuantiserTable(Target* target, int _wIn, int _wOut, const std::vector<mpz_class> &values, int logicTable = 0,  map<string, double> inputDelays = emptyDelayMap )
			: Table(target, _wIn, _wOut, 0, values.size()-1, logicTable, inputDelays)
			, m_values(values)
		{
			ostringstream name;
			name << "StaticQuantiser_Table"<<"_uid"<<getNewUId();
			setName(name.str());
		}

		virtual mpz_class function(int x)
		{
			return m_values.at(x);
		}
	};
	
	if(level==1){
		vhdl << declare(join("bit_",level),1) << " <= \"0\" when (iX < \""<<unsignedBinary(GetBoundary(1, 0),m_wValue)<<"\" ) else \"1\";\n";
		vhdl << declare(join("idx_",level),level) << " <= "<< join("bit_",level)<<";\n";
	}else{
		std::vector<mpz_class> boundaries(1<<(level-1));
		for(int i=0;i<(1<<(level-1));i++){
			boundaries.at(i)=GetBoundary(level, i);
		}
		
		bool useBRAM=level > 9;	// level==9 implies 2^8 elements, so will use LUTs for up to 256 elements, then switch to BRAM for 512
		
		StaticQuantiserTable *table=new StaticQuantiserTable(target_, level-1, m_wValue, boundaries, useBRAM?1:0, inDelayMap( "X", target_->localWireDelay() + getCriticalPath() ));
		oplist.push_back(table);
		inPortMap(table, "X", join("idx_",level-1));
		outPortMap(table, "Y", join("thresh_",level));
		vhdl << instance(table, join("boundaries_",level));
		syncCycleFromSignal(join("thresh_",level), table->getOutputDelay("Y") );
		if(useBRAM){
			nextCycle();	// Force the output register
		}

		manageCriticalPath(target_->adderDelay(m_wValue));
		vhdl << declare(join("bit_",level),1) << " <= \"0\" when (iX < "<<join("thresh_",level)<<" ) else \"1\";\n";
		
		// No delay, just joining
		vhdl << declare(join("idx_",level),level) << " <= " << join("idx_",level-1) <<" & "<<join("bit_",level)<<";\n";
	}
}

mpz_class StaticQuantiser::GetBoundary(int level, int index)
{
	int nb=boundaries.size();	// Should be a binary power
	// level == 1 -> index==0					address=nb/2
	// level == 2 -> index==0,1				address=nb/4,3*nb/4
	// level == 3 -> index==0,1,2,3			address=nb/8,3*nb/8,5*nb/8,7*nb/8
	// level == 4 -> index==0,1,2,3,4,5,6,7	address=nb/16,3*nb/16, ... 15*nb/16
	//												address= (index*2+1) * nb/(2^level)
	int address=(index*2+1) * nb / (1<<level);
	return boundaries.at(address);
}

  void StaticQuantiser::emitHLSBody
  (
   HLSContext &ctxt,
   HLSScope &scope
   )const
  {
    auto iX=hls_get("iX");
    auto oY=hls_get("oY");
    oY=hls_og(oY.getNode()->getType()->getWidth());
  }

// Apply a quantisation of the form y=floor( f(x)*n )

mpfr::mpreal eval(const Function *f, mpfr::mpreal x)
{
	mpfr::mpreal tmp(0, getToolPrecision());
	f->eval(get_mpfr_ptr(tmp), get_mpfr_ptr(x));
	return tmp;
}

Operator *StaticQuantiser::BuildFloatQuantiser(Target *target, int wE, int wF, int n, const Function *f, map<string, double> inputDelays)
{
	std::vector<mpz_class> boundaries(n+1);
	
	mpfr::mpreal a(-pow(2.0,wE+1),getToolPrecision());
	mpfr::mpreal b=-a;	// Use interval over all representable values...
	
	mpfr::mpreal x(0, getToolPrecision());
	
	ComparableFloatType type(wE,wF);
	
	boundaries[0]=mpz_class(0);
	for(int i=1;i<n;i++){
		// Values in bucket i  have to lie in the range [ f(x)*(i-1)..f(x)*i ), so the
		// upper boundary of bucket (i-1) is at f(x)=i/n
		
		mpfr::mpreal y(i, getToolPrecision());
		y=y/n;
		f->eval_inverse(get_mpfr_ptr(x), get_mpfr_ptr(y), get_mpfr_ptr(a), get_mpfr_ptr(b));
		
		boundaries[i]=type.ToBits(get_mpfr_ptr(x), true);	// round and convert
		if(::flopoco::verbose>=DETAILED){
			std::cerr<<" Quantiser["<<i<<"] : f("<<x<<")="<<y<<", bits="<<std::hex<<boundaries[i]<<std::dec<<"\n";
		}
	}
	boundaries[n]=(mpz_class(1)<<(type.Width()))-1;
	
	for(int i=0;i<n;i++){
		if(boundaries[i] > boundaries[i+1]){

			throw std::string("BuildFloatQuantiser - Function does not result in monotonically increasing intervals.");
		}
	}
	
	Operator *encoder=type.MakeEncoder(target, inputDelays);
	
	map<string, double> interDelays;
	interDelays["iX"]=encoder->getOutputDelay("oY");
	Operator *quantiser=new StaticQuantiser(target, type.Width(), boundaries, interDelays);
	
	std::stringstream acc;
	acc<<"StaticQuantiser_wE"<<wE<<"_wF"<<wF<<"_n"<<(boundaries.size()-1)<<"_uid"<<Operator::getNewUId();
	ChainOperator::mapping_list_t mapping;
	mapping.push_back(ChainOperator::mapping_t("oY","iX"));
	
	ChainOperator *res=ChainOperator::Create(acc.str(), encoder, mapping, quantiser);
	
	FloPoCoRandomState::init(1);	// Eventually this will move to main?
	
	FPNumber fpRep(wE,wF);
	mpfr::mpreal vv=mpfr::create_zero(wF+1);
	for(unsigned i=1;i<boundaries.size();i++){
		mpz_class a=boundaries[i-1], b=boundaries[i];
		
		type.FromBits(get_mpfr_ptr(vv), a);
		fpRep=get_mpfr_ptr(vv);
		res->addStandardTestCaseInputs("iX", fpRep.getSignalValue());
	
		type.FromBits(get_mpfr_ptr(vv), b);
		fpRep=get_mpfr_ptr(vv);
		res->addStandardTestCaseInputs("iX", fpRep.getSignalValue());
		
		if(b>a){
			for(int i=0;i<10;i++){
				mpz_class v=getRandomBetween(a,b);
				type.FromBits(get_mpfr_ptr(vv), v);
				fpRep=get_mpfr_ptr(vv);
				res->addStandardTestCaseInputs("iX", fpRep.getSignalValue());
			}
		}
	}
	
	return res;
}

static void StaticQuantiserFactoryUsage(std::ostream &dst)
{
	OperatorFactory::classic_OP(dst, "static_quantiser", "wE wF n f", false);
	dst << "    Generates a non-linear quantiser that calculates y=floor(f(x)*n).\n";
	dst << "	      wE - Exponent bits of input.\n";
	dst << "	      wF - Fractional bits of each table element\n";
	dst << "	      n - How many levels to quantise into (storage requirements are O(n))\n";
	dst << "      The transform will take a (wE,wF) input float, and produced a log2ceil(n) bit unsigned output.\n";
}

static Operator *StaticQuantiserFactoryParser(Target *target ,const std::vector<std::string> &args,int &consumed)
{
	unsigned nargs = 4;
	if (args.size()<nargs)
		throw std::string("StaticQuantiserFactoryParser - Not enough arguments, check usage.");
	consumed += nargs;
	
	if(::flopoco::verbose>=INFO){
		std::cerr<<"StaticQuantiserFactoryParser(wE='"<<args[0]<<"',wF='"<<args[1]<<"',n='"<<args[2]<<"',f='"<<args[3]<<"\n";
	}
	
	int wE = atoi(args[0].c_str());
	int wF = atoi(args[1].c_str());
	int n=atoi(args[2].c_str());
	Function f(args[3]);

	if(wE<1)
		throw std::string("StaticQuantiserFactoryParser - wE must be a positive integer.");
	if(wF<1)
		throw std::string("StaticQuantiserFactoryParser - wF must be a positive integer.");
	if(n<1)
		throw std::string("StaticQuantiserFactoryParser - n must be a positive integer.");
	
	if((n>16384) && (::flopoco::verbose>=INFO)){
		std::cerr<<"Warning : Building a quantiser with "<<n<<" output levels, this will require huge amounts of RAM.\n";
	}
	
	return StaticQuantiser::BuildFloatQuantiser(target, wE, wF, n, &f);
}

void StaticQuantiser_registerFactory()
{
	DefaultOperatorFactory::Register(
		"static_quantiser",
		"operator",
		StaticQuantiserFactoryUsage,
		StaticQuantiserFactoryParser,
		DefaultOperatorFactory::ParameterList(
			DefaultOperatorFactory::Parameters("6", "12", "1.0", "auto", "auto")
		)
	);
}

}; 
};
