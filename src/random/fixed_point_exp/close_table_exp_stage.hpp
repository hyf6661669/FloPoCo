#ifndef flopoco_random_close_table_exp_stage_hpp
#define flopoco_random_close_table_exp_stage_hpp

#include "fixed_point_exp_stage.hpp"

#include "Table.hpp"

#include "find_close_value.hpp"

namespace flopoco
{
namespace random
{

/* Creates a table that maps (residual_hi@residual_lo) to
	(residual_lo+offset) and result_hi.
	The width of result_hi is determined by (rel_err_min,rel_err_max),
	and will be increased until it can be met.
*/
class CloseTableExpStage
	: public FixedPointExpStage
{
protected:
	residual_type<T> m_outputResidualType;
	result_type<T> m_outputResultType;
	residual_type<T> m_tableIndexType;
	residual_type<T> m_offsetResidualType;
	residual_type<T> m_baseResidualType;

	int m_tableInputBits, m_tableEntryBits;
	int m_tableMinIndex, m_tableMaxIndex;
	
	struct table_entry{
		table_entry()
			: inRange(false)
			, input(0)
			, trueResult(1)
			, roundedResult(1)
		{}
		
		bool inRange;	// not all entries are included in residual range
		T input;
		T residualOffset;	// What to pass on
		T trueResult;
		T roundedResult;
	};
	std::vector<table_entry> m_table;
	
	int TableInputBits() const { return m_tableInputBits; }
	int TableEntryBits() const { return m_tableEntryBits; }
	int TableMinIndex() const { return m_tableMinIndex; }
	int TableMaxIndex() const { return m_tableMaxIndex; }
	
	mpz_class TableEntry(int index) const
	{
		assert(m_table.at(index).inRange);
		uint64_t lo= m_outputResultType.ToBits(m_table[index].roundedResult);
		uint64_t hi=m_offsetResidualType.ToBits(m_table[index].residualOffset);
		return to_mpz_class(lo) + (to_mpz_class(hi)<<m_outputResultType.Width());
	}
	
	
	class TableOp : public flopoco::Table
	{
	private:	
		CloseTableExpStage *m_parent;
	public:
		TableOp(Target *target,
			CloseTableExpStage *parent, int logicTable = 0,  map<string, double> inputDelays = emptyDelayMap)
			: Table(target, 
				parent->TableInputBits(), parent->TableEntryBits(),
				parent->TableMinIndex(),parent->TableMaxIndex(),
				logicTable, inputDelays
			)
			, m_parent(parent)
		{
			setName(parent->getName()+"_table");
		}
		
		virtual mpz_class function(int x)
		{return m_parent->TableEntry(x);}
	};
public:
	CloseTableExpStage(
		Target *target,
		T expMu,
		T expSigma,
		residual_type<T> inputResidualType,
		int addressBits,
		T errorMin, T errorMax,
		map<string, double> inputDelays = emptyDelayMap
	)
		: FixedPointExpStage(target, inputDelays, expMu, expSigma, inputResidualType, result_type<T>(0))
		, m_outputResidualType(addressBits==inputResidualType.Width() ? inputResidualType.drop_msbs(addressBits) : inputResidualType.drop_msbs(addressBits-1))
		, m_outputResultType(0)
		, m_tableIndexType(inputResidualType.take_msbs(addressBits))
		, m_offsetResidualType(inputResidualType.drop_msbs(addressBits))
		, m_baseResidualType(inputResidualType.drop_msbs(addressBits))
	{
		if(addressBits<2)
			throw std::invalid_argument("CloseTableExpStage - Can't have less than two address bits.");
		
		m_table.resize(1<<addressBits);
		typename residual_type<T>::iterator curr=m_tableIndexType.begin_msbs(addressBits);
		typename residual_type<T>::iterator end=m_tableIndexType.end_msbs(addressBits);
		
		close_value_table_t<T> ct=FindCloseValues(
			expMu, expSigma,
			m_tableIndexType,
			inputResidualType.drop_msbs(addressBits),
			errorMin, errorMax
		);
		m_offsetResidualType=ct.offsetSpace;
		m_outputResultType=ct.resultSpace;
		
		m_outputResidualType.ResetRange();
		m_outputResidualType.Add(m_baseResidualType.RangeMin()+m_offsetResidualType.RangeMin());
		m_outputResidualType.Add(m_baseResidualType.RangeMax()+m_offsetResidualType.RangeMax());
		
		m_tableInputBits=m_tableIndexType.Width();
		m_tableEntryBits=m_outputResultType.Width()+m_offsetResidualType.Width();
		
		std::cerr<<"tableEntryBits="<<m_tableEntryBits<<"\n";
		std::cerr<<"tableIndexType="<<m_tableIndexType<<"\n";
		std::cerr<<"offsetType="<<m_offsetResidualType<<"\n";
		std::cerr<<"resultType="<<m_outputResultType<<"\n";
		
		m_tableMinIndex=INT_MAX;
		m_tableMaxIndex=0;
		while(curr!=end){
			T x=*curr;
			bool valid=inputResidualType.IsInRange(x);
			if(valid){
				unsigned index=m_tableIndexType.ToBits(x);
				m_tableMinIndex=std::min(m_tableMinIndex, (int)index);
				m_tableMaxIndex=std::max(m_tableMaxIndex, (int)index);
				m_table.at(index).trueResult=ct.table[x].exact;
				m_table.at(index).roundedResult=ct.table[x].approx;
				m_table.at(index).residualOffset=ct.table[x].offset;
				m_table.at(index).inRange=true;
			}
			++curr;
		}
		
		ostringstream name;
		name << "CloseTableExpStage_"<<m_inputResidualType.DescriptionId();
		setName(name.str());
		setCopyrightString("Imperial College 2012");

		addInput ("iResidual" , m_inputResidualType.Width(), true);
		//addInput ("iResult" , m_inputResultType.Width());
		if(m_outputResidualType.Width()>0){
			addOutput ("oResidual" , m_outputResidualType.Width(), 1, true);
		}
		addOutput ("oResult" , m_outputResultType.Width(), 1, true);
		
		TableOp *table=new TableOp(target, this, 0 /*logicTable*/, inputDelays);
		AddOperator(table); // add to oplist
		
		vhdl << tab  << declare("tableAddr", addressBits, true) << " <= iResidual" << range(m_inputResidualType.Width()-1, m_inputResidualType.Width()-addressBits) << ";" << endl;
		inPortMap(table, "X", "tableAddr");
		outPortMap(table, "Y", "tableData");
		vhdl << instance(table, "table");
		syncCycleFromSignal("tableData");
		
		if(m_offsetResidualType.LsbPower()!=m_baseResidualType.LsbPower())
			throw std::logic_error("CloseTableExpStage - LSBs of base and offset residuals don't match.");
		
		if(m_outputResidualType.Width()>0){	
			vhdl<<tab<<declare("residual_base", m_outputResidualType.Width(), true) << " <= '0'&iResidual("<<m_baseResidualType.Width()-1<<" downto 0);\n";
			vhdl<<tab<<declare("residual_offset", m_outputResidualType.Width(), true) << " <= '0'&tableData("<<m_tableEntryBits-1<<" downto "<<m_outputResultType.Width()<<");\n";
			vhdl<<tab<<declare("residual_next", m_outputResidualType.Width(), true) << " <= residual_base+residual_offset;\n";
		}
		
		nextCycle();
		
		// Technically this can go out a cycle earlier, if later components are careful
		vhdl << tab << "oResult <= tableData("<<m_outputResultType.Width()-1<<" downto 0);\n";
		if(m_outputResidualType.Width()>0){
			vhdl<<tab<<"oResidual <= residual_next;\n";
		}
	}
	
	virtual ~CloseTableExpStage()
	{}
	
	virtual residual_type<T> OutputResidualType() const
	{ return m_outputResidualType; }
	
	virtual result_type<T> OutputResultType() const
	{ return m_outputResultType; }
	
	virtual std::pair<T,T> Execute(const T &residual, const T &result) const
	{
		assert(result==1.0);
		uint64_t bits=m_inputResidualType.ToBits(residual);
		
		unsigned index=bits>>(m_inputResidualType.Width()-m_tableInputBits);
		T baseResidual=m_outputResidualType.FromBits(bits & ((1ull<<(m_inputResidualType.Width()-m_tableInputBits))-1));
		
		assert(m_tableIndexType.FromBits(index)+baseResidual==residual);
		
		return std::make_pair(baseResidual+m_table[index].residualOffset, m_table[index].roundedResult);
	}
};

}; // random
}; // flopoco

#endif
