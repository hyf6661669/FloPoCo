#ifndef flopoco_random_chained_exp_stage_hpp
#define flopoco_random_chained_exp_stage_hpp

#include "fixed_point_exp_stage.hpp"

#include "Table.hpp"

#include <boost/lexical_cast.hpp>

namespace flopoco
{
namespace random
{

class ChainedExpStage
	: public FixedPointExpStage
{
protected:
	std::vector<FixedPointExpStagePtr> m_stages;
	
	T CalcMu(const std::vector<FixedPointExpStagePtr> &stages)
	{
		T acc=0;
		for(int i=0;i<(int)stages.size();i++)
			acc+=stages[i]->Mu();
		return acc;
	}
public:
	ChainedExpStage(
		Target *target,
		std::vector<FixedPointExpStagePtr> stages,
		const std::map<std::basic_string<char>, double> &inputDelays=emptyDelayMap
	)
		: FixedPointExpStage(target, inputDelays, CalcMu(stages), stages[0]->Sigma(), stages[0]->InputResidualType(), stages[0]->InputResultType())
		, m_stages(stages)
	{
		ostringstream name;
		name << "ChainedExpStage";
		for(unsigned i=0;i<stages.size();i++){
			name<<"_"<<stages[i]->getName();
		}
		setName(shortenId(name.str()));
		setCopyrightString("Imperial College 2012");

		addInput ("iResidual" , InputResidualType().Width(), false);
		if(InputResultType().Width()!=0){
			addInput ("iResult" , InputResultType().Width(), false);
		}
		if(OutputResidualType().Width()!=0){
			addOutput ("oResidual" , OutputResidualType().Width(), 1, false);
		}
		addOutput ("oResult" , OutputResultType().Width(), 1, false);
		
		std::string residual_prev=declare("residual_0", InputResidualType().Width(), true);
		vhdl<<tab<<residual_prev<<" <= iResidual;\n";
		std::string result_prev;
		if(InputResultType().Width()!=0){
			result_prev=declare("result_0", InputResultType().Width(), true);
			vhdl<<tab<<result_prev<<" <= iResult;\n";
		}
		
		for(int i=0;i<(int)stages.size();i++){
			inPortMap(stages[i].get(), "iResidual", residual_prev);
			if(result_prev.size()){
				inPortMap(stages[i].get(), "iResult", result_prev);
			}
			result_prev="result_"+boost::lexical_cast<std::string>(i+1);
			outPortMap(stages[i].get(), "oResult", result_prev);
			if(stages[i]->OutputResidualType().Width()){
				residual_prev="residual_"+boost::lexical_cast<std::string>(i+1);
				outPortMap(stages[i].get(), "oResidual", residual_prev);
			}else{
				residual_prev.clear();
			}
			vhdl<<tab<<instance(stages[i].get(), "stage_"+boost::lexical_cast<std::string>(i));
			syncCycleFromSignal(result_prev);
			if(residual_prev.size()!=0){
				syncCycleFromSignal(residual_prev);
			}
		}
		
		vhdl<<tab<<"oResult <= "<<result_prev<<";\n";
		if(residual_prev.size()){
			vhdl<<tab<<"oResidual <= "<<residual_prev<<";\n";
		}
	}
	
	virtual ~ChainedExpStage()
	{}
	
	virtual residual_type<T> OutputResidualType() const
	{ return m_stages.back()->OutputResidualType(); }
	
	virtual result_type<T> OutputResultType() const
	{ return m_stages.back()->OutputResultType(); }
	
	virtual std::pair<T,T> Execute(const T &residual, const T &result) const
	{
		std::pair<T,T> curr(residual, result);
		for(unsigned i=0;i<m_stages.size();i++){
			curr=m_stages[i]->Execute(curr.first, curr.second);
		}
		return curr;
	}
};

}; // random
}; // flopoco

#endif
