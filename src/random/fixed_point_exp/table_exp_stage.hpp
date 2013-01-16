#ifndef flopoco_random_table_exp_stage_hpp
#define flopoco_random_table_exp_stage_hpp

#include "fixed_point_exp_stage.hpp"

#include "Table.hpp"

namespace flopoco
{
namespace random
{

class TableExpStage
	: public FixedPointExpStage
{
protected:
	residual_type<T> m_outputResidualType;
	result_type<T> m_outputResultType;
	residual_type<T> m_tableIndexType;

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
		uint64_t x= m_outputResultType.ToBits(m_table[index].roundedResult);
		mpz_class res(uint32_t(x>>32));
		return (res<<32) + (uint32_t)(x&0xFFFFFFFF);
	}
	
	
	class TableOp : public flopoco::Table
	{
	private:	
		TableExpStage *m_parent;
	public:
		TableOp(Target *target,
			TableExpStage *parent, int logicTable = 0,  map<string, double> inputDelays = emptyDelayMap)
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
	TableExpStage(
		Target *target,
		T expMu,
		T expSigma,
		residual_type<T> inputResidualType,
		int addressBits,
		int resultFracBits,	// explicit for this version, we just round to it
		int logicTable = 0, 
		map<string, double> inputDelays = emptyDelayMap
	)
		: FixedPointExpStage(target, inputDelays, expMu, expSigma, inputResidualType, result_type<T>(0))
		, m_outputResidualType(inputResidualType.drop_msbs(addressBits))
		, m_outputResultType(resultFracBits)
		, m_tableIndexType(inputResidualType.take_msbs(addressBits))
	{
		
		m_table.resize(1<<addressBits);
		residual_type<T>::iterator curr=m_tableIndexType.begin_msbs(addressBits);
		residual_type<T>::iterator end=m_tableIndexType.end_msbs(addressBits);
		
		m_tableMinIndex=INT_MAX;
		m_tableMaxIndex=0;
		while(curr!=end){
			T x=*curr;
			bool valid=inputResidualType.IsInRange(x);
			if(valid){
				unsigned index=m_tableIndexType.ToBits(x);
				T fx=ReferenceExp(x);
				T rx=m_outputResultType.Round(fx);
				assert(rx=m_outputResultType.Round(rx));
				m_outputResultType.Add(rx, fx);
				m_tableMinIndex=std::min(m_tableMinIndex, (int)index);
				m_tableMaxIndex=std::max(m_tableMaxIndex, (int)index);
				m_table.at(index).trueResult=fx;
				m_table.at(index).roundedResult=rx;
				m_table.at(index).inRange=true;
			}
			++curr;
		}
		
		m_tableInputBits=m_tableIndexType.Width();
		m_tableEntryBits=m_outputResultType.Width();
		
		std::cerr<<"tableEntryBits="<<m_tableEntryBits<<"\n";
		std::cerr<<"resultType="<<m_outputResultType<<"\n";
		
		ostringstream name;
		name << "TableExpStage_"<<m_inputResidualType.DescriptionId();
		setName(name.str());
		setCopyrightString("Imperial College 2012");

		addInput ("iResidual" , m_inputResidualType.Width(), true);
		//addInput ("iResult" , m_inputResultType.Width());
		if(m_outputResidualType.Width()>0){
			addOutput ("oResidual" , m_outputResidualType.Width(), 1, true);
		}
		addOutput ("oResult" , m_outputResultType.Width(), 1, true);
		
		TableOp *table=new TableOp(target, this, logicTable, inputDelays);
		AddOperator(table); // add to oplist
		
		vhdl << tab  << declare("tableAddr", addressBits, true) << " <= iResidual" << range(m_inputResidualType.Width()-1, m_inputResidualType.Width()-addressBits) << ";" << endl;
		inPortMap(table, "X", "tableAddr");
		outPortMap(table, "Y", "tableData");
		vhdl << instance(table, "table");
		
		nextCycle();
		
		if(m_outputResidualType.Width()>0){
			vhdl<<tab<<declare("residual_temp", m_inputResidualType.Width(), true) << " <= iResidual;\n";
			syncCycleFromSignal("residual_temp");
			vhdl << tab << "oResidual <= residual_temp"<<range(m_inputResidualType.Width()-addressBits-1, 0)<<";\n";
		}
		syncCycleFromSignal("tableData");
		vhdl << tab << "oResult <= tableData;\n";
	}
	
	virtual ~TableExpStage()
	{}
	
	virtual residual_type<T> OutputResidualType() const
	{ return m_outputResidualType; }
	
	virtual result_type<T> OutputResultType() const
	{ return m_outputResultType; }
	
	virtual std::pair<T,T> Execute(const T &residual, const T &result) const
	{
		assert(result==1.0);
		uint64_t bits=m_inputResidualType.ToBits(residual);
		
		unsigned index=bits>>m_outputResidualType.Width();
		T nextResidual=m_outputResidualType.FromBits(bits & ((1ull<<m_outputResidualType.Width())-1));
		
		assert(m_tableIndexType.FromBits(index)+nextResidual==residual);
		
		return std::make_pair(nextResidual, m_table[index].roundedResult);
	}
	
	virtual std::vector<std::pair<T,T> > GetAllResultValues() const
	{
		std::vector<std::pair<T,T> > res;
		res.reserve(m_table.size());
		
		for(int i=0;i<(int)m_table.size();i++){
			if(m_table[i].inRange){
				res.push_back(std::make_pair(m_table[i].trueResult, m_table[i].roundedResult));
			}
		}
		
		return res;
	}
};

}; // random
}; // flopoco

#endif
