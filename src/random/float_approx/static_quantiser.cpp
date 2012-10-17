#include "static_quantiser.hpp"
#include <assert.h>

#include "Table.hpp"

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
	
	addInput ("iX" , wValue);
	addOutput("oY", wIndex);
	
	setCriticalPath( getMaxInputDelays(inputDelays) );
	for(int i=1;i<=wIndex;i++){
		AddLevel(i);
	}
	
	vhdl << "index <= "<<join("idx_",wIndex)<<";\n";
	outDelayMap["index"] = getCriticalPath();
}

StaticQuantiser::~StaticQuantiser()
{}


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
		vhdl << declare(join("bit_",level),1) << " <= '0' when (x < "<<unsignedBinary(GetBoundary(1, 0),m_wValue)<<" ) else '1';\n";
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
		vhdl << declare(join("bit_",level),1) << " <= '0' when (x < "<<join("thresh_",level)<<" ) else '1';\n";
		
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

/*
Operator *StaticQuantiser::BuildFloatQuantiser(Target *target, int wE, int WF, const MPFRVec &boundaries, map<string, double> inputDelays = emptyDelayMap)
{
	std::vector<mpz_class> fixBoundaries;
	
	ComparableFloatType type(wE, wF);
	for(unsigned i=1;i<boundaries.size()-1;i++){
		fixBoundaries.push_back(type.ToBits(boundaries[i]));
	}
	
	Operator *encoder=type.MakeEncoder(target, inputDelays);
	oplist.push_back(encoder);
	
	map<string, double> interDelays;
	interDelays["iX"]=encoder->getOutputDelay("oY");
	Operator *quantiser=new StaticQuantiser(target, type.GetWidth(), fixBoundaries, interDelays);
	oplist.push_back(quantiser);
	
	std::strstream acc;
	acc<<"StaticQuantiser_wE"<<wE<<"_wF"<<wF<<"_n"<<(boundaries.size()-1)<<"_uid"<<getNewUid();
	Operator *res=ChainOperator::Create(acc.str(), encoder, quantiser);
	oplist.push_back(res);
	
	return res;
}
*/

// Apply a quantisation of the form y=floor( f(x)*n )
/*
Operator *StaticQuantiser::BuildFloatQuantiser(Target *target, int wE, int WF, int n, Function *f, map<string, double> inputDelays = emptyDelayMap)
{
	MPFRVec vec(wF, n+1);
	
	mpfr_t tmp, a, b;
	mpfr_inits2(getToolPrecision(), a,b,tmp ,(mpfr_ptr)0);
	
	mpfr_set_d(a,0.0, MPFR_RNDN);
	
	for(int i=1;i<n;i++){
		// Values in bucket i  have to lie in the range [ f(x)*(i-1)..f(x)*i ), so the
		// upper boundary of bucket (i-1) is at f(x)=1/i
		mpfr_set_d(tmp, (double)i, MPFR_RNDN);
		mpfr_div_d(tmp, tmp, (double)n, MPFR_RNDN);
		
		f->eval_inverse(vec[i], tmp, 
	}
	
	Operator *encoder=type.MakeEncoder(target, inputDelays);
	oplist.push_back(encoder);
	
	map<string, double> interDelays;
	interDelays["iX"]=encoder->getOutputDelay("oY");
	Operator *quantiser=new StaticQuantiser(target, type.GetWidth(), fixBoundaries, interDelays);
	oplist.push_back(quantiser);
	
	std::strstream acc;
	acc<<"StaticQuantiser_wE"<<wE<<"_wF"<<wF<<"_n"<<(boundaries.size()-1)<<"_uid"<<getNewUid();
	Operator *res=ChainOperator::Create(acc.str(), encoder, quantiser);
	oplist.push_back(res);
	
	return res;
}
*/

}; 
};
