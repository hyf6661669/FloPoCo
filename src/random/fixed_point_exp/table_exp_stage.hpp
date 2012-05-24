#ifndef flopoco_random_fixed_point_exp_stage_hpp
#define flopoco_random_fixed_point_exp_stage_hpp

#include "Operator.hpp"
#include "utils.hpp"

namespace flopoco
{
namespace random
{

class TableExpStage
	: public FixedPointExpStage
{
protected:
	residual_type m_outputResidualType;
	result_type m_outputResultType;

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
		return m_outputResultType.AsBits(m_table[roundedResult]);
	}
	
	
	class TableOp : public flopoco::Table
	{
	public:
		TableOp(Target *target, TableExpStage *parent, int logicTable = 0,  map<string, double> inputDelays = emptyDelayMap)
			: Table(target, 
				m_parent->TableInputBits(), m_parent->TableEntryBits(),
				m_parent->TableMinIndex(),m_parent->TableMaxIndex(),
				logicTable, inputDelays
		){}
		
		virtual mpz_class function(int x)
		{return m_parent->TableEntry(x);}
	};
public:
	TableExpStage(
		Target *target,
		residual_type inputResidualType,
		int addressBits,
		int resultFracBits,	// explicit for this version, we just round to it
		int logicTable = 0, 
		map<string, double> inputDelays = emptyDelayMap
	)
		: FixedPointExpStage(target, inputResidualType, result_type())
		, m_outputResidualType(m_inputResidualType.drop_msbs(addressBits))
		, m_outputResultType(resultFracBits)
	{
		m_table.resize(1<<addressBits);
		residual_type::iterator curr=m_inputResidualType.begin_msbs(addressBits);
		residual_type::iterator end=m_inputResidualType.end_msbs(addressBits);
		
		int dst=0;
		m_tableMinIndex=INT_MAX;
		m_tableMaxIndex=0;
		while(curr!=end){
			value_t x=*curr;
			bool valid=m_inputResidualType.IsInRange(x);
			if(valid){
				value_t fx=ReferenceExp(x);
				value_t rx=m_outputResultType.Round(fx);
				m_outputResultType.Add(fx, rx);
				m_tableMinIndex=std::min(m_tableMinIndex, dst);
				m_tableMaxIndex=std::max(m_tableMaxIndex, dst);
				m_table[dst].isValid=true;
				m_table[dst].trueValue=fx;
				m_table[dst].roundedValue=rx;
			}
			++dst;
			++curr;
		}
		
		m_tableInputBits=
		
		ostringstream name;
		name << "TableExpStage_"<<m_inputResidualType.DescriptionId();
		setName(name.str());
		setCopyrightString("Imperial College 2012");

		addInput ("iResidual" , m_inputResidualType.Width());
		addInput ("iResult" , m_inputResultType.Width()-1);
		addOutput ("oResidual" , m_outputResidualType.Width());
		addOutput ("oResult" , m_outputResultType.Width()-1);
		
		TableOp *table=new TableOp(target, this, logicTable, inputDelays);
		AddOperator(table); // add to oplist
		
		vhdl << tab  << declare("tableAddr", addressBits) << " <= iResidual" << range(m_inputResidualType.Width()-1, m_inputResidualType.Width()-addressBits) << ";" << endl;
		inPortMap(table, "R", "tableAddr");
		outPortMap(lshift, "S", "tableData");
		vhdl << instance(table, "table");
		syncCycleFromSignal("tableData");
		
		vhdl <<tab<< declare("unused_signal", m_inputResultType.Width()-1) << iResult<<"\n";
		
		vhdl << tab << "oResidual <= iResidual"<<range(m_inputResidualType.Width()-addressBits-1, downto 0)<<";\n";
		vhdl << tab << "oResult <= tableData;\n";
	}
	
	virtual residual_type OutputResidualType() const
	{ return m_outputResidualType; }
	
	virtual result_type OutputResultType() const
	{ return m_outputResultType; }
	
	virtual std:pair<value_t,value_t> Execute(const value_t &residual, const value_t &result) const
	{
		assert(result==1.0);
		unsigned index=residual.TakeMsbsAsBits(residual);
		value_t nextResidual=residual.DropBsbs(residual);
		assert(m_table[index].input+nextResidual==residual);
		return std::make_pair(nextResidual, m_table[index].roundedResult);
	}
};

}; // random
}; // flopoco

#endif
