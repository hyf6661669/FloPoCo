#ifndef flopoco_random_multiplier_v2_exp_stage_hpp
#define flopoco_random_multiplier_v2_exp_stage_hpp

#include "fixed_point_exp_stage.hpp"

#include "Table.hpp"

namespace flopoco
{
namespace random
{

class MultiplierV2ExpStage
	: public FixedPointExpStage
{
protected:
	FixedPointExpStagePtr m_impl;

	residual_type<T> m_outputResidualType;
	result_type<T> m_outputResultType;


public:
	MultiplierV2ExpStage(
		Target *target,
		FixedPointExpStagePtr impl,
		result_type<T> inputResultType,
		unsigned outputFracBits,
		map<string, double> inputDelays = emptyDelayMap
	)
		: FixedPointExpStage(target, inputDelays, impl->Mu(), impl->Sigma(), impl->InputResidualType(), inputResultType)
		, m_impl(impl)
		, m_outputResidualType(impl->OutputResidualType())
		, m_outputResultType(std::min(outputFracBits, inputResultType.FracWidth()+impl->OutputResultType().FracWidth()))
	{
		std::cerr<<"MultiplierV2 src type : "<<inputResultType<<"\n";
		std::cerr<<"Stage res type : "<<impl->OutputResultType()<<"\n";
		std::cerr<<"range = "<<inputResultType.RangeMin()*impl->OutputResultType().RangeMin()<<", "<<inputResultType.RangeMax()*impl->OutputResultType().RangeMax()<<"\n";
		m_outputResultType.Add(m_outputResultType.Floor(inputResultType.RangeMin()*impl->OutputResultType().RangeMin()));
		m_outputResultType.Add(m_outputResultType.Floor(inputResultType.RangeMax()*impl->OutputResultType().RangeMax()));
		std::cerr<<"Multiplier output type : "<<m_outputResultType<<"\n";
		
		ostringstream name;
		name << "MultiplierExpV2Stage_r"<<outputFracBits<<"_"<<impl->getName();
		setName(name.str());
		setCopyrightString("Imperial College 2012");
		
		int FW_A=InputResultType().FracWidth();
		int FW_B=impl->OutputResultType().FracWidth();
		int AW_B=FW_B - (impl->OutputResultType().FracMax()==0.5 ? 0 : ceil(-log(impl->OutputResultType().FracMax()-0.5)));
		int AW_A=FW_A - (InputResultType().FracMax()==0.5 ? 0 : ceil(-log(InputResultType().FracMax()-0.5)));
		
		addConstant("FW_A", "NATURAL",  FW_A);
		addConstant("FW_B", "NATURAL", FW_B);
		addConstant("AW_B", "NATURAL", AW_B);
		addConstant("AW_A", "NATURAL", AW_A);
		
		/* If the stage output always has a fixed exponent, we can treat the multiplication as (1+b), and
			calculate 2^a_e  * a_f*(1+b_f) == 2^a_e * (a_f + a_f * b_f)
			For us to care, we want the fraction to be of the form 0.10XXXXb, so need the fraction to
			be < 0.11b == 0.75.
		*/
		bool splitMultiplier = (AW_A+1 < FW_A) || (AW_B+1<FW_B);

		addInput ("iResidual" , InputResidualType().Width(), true);
		addInput ("iResult" , InputResultType().Width(), true);
		if(OutputResidualType().Width()!=0){
			addOutput ("oResidual" , OutputResidualType().Width(), 1, true);
		}
		addOutput ("oResult" , OutputResultType().Width(), 1, true);
		
		//declare("stage_result", impl->OutputResultType().Width());
		//declare("stage_residual", impl->OutputResidualType().Width());
		
		vhdl<<tab<<declare("input_result", InputResultType().Width(), true)<<" <= iResult;\n";
		
		inPortMap(impl.get(), "iResidual", "iResidual");
		assert(impl->InputResultType().Width()==0);
		// Don't have iResult, as it is zero width
		if(OutputResidualType().Width()!=0){
			outPortMap(impl.get(), "oResidual", "stage_residual");
		}
		outPortMap(impl.get(), "oResult", "stage_result");
		vhdl<<tab<<instance(impl.get(), "stage");
		
		syncCycleFromSignal("input_result");
		syncCycleFromSignal("stage_result");
		
		if(impl->OutputResultType().ExpWidth()>0){
			vhdl<<tab<<declare("stage_result_exp", impl->OutputResultType().ExpWidth(), true)<<" <= stage_result("<<impl->OutputResultType().Width()-1<<" downto "<<impl->OutputResultType().FracWidth()<<");\n";
		}
		vhdl<<tab<<declare("stage_result_frac", impl->OutputResultType().FracWidth()+1, true)<<" <= '1' & stage_result("<<impl->OutputResultType().FracWidth()-1<<" downto 0);\n";
		
		if(InputResultType().ExpWidth()>0){
			vhdl<<tab<<declare("input_result_exp", InputResultType().ExpWidth(), true)<<" <= input_result("<<InputResultType().Width()-1<<" downto "<<InputResultType().FracWidth()<<");\n";
		}
		vhdl<<tab<<declare("input_result_frac", InputResultType().FracWidth()+1, true)<<" <= '1' & input_result("<<InputResultType().FracWidth()-1<<" downto 0);\n";
		
		unsigned partialFracWidth=2+impl->OutputResultType().FracWidth()+InputResultType().FracWidth();
		if((impl->OutputResultType().ExpWidth()>0) && (InputResultType().ExpWidth()>0)){
			vhdl<<tab<<declare("result_exp", OutputResultType().ExpWidth(), true) <<" <= stage_result_exp + input_result_exp;\n";
		}else if(impl->OutputResultType().ExpWidth()>0){
			vhdl<<tab<<declare("result_exp", OutputResultType().ExpWidth(), true) <<" <= stage_result_exp;\n";
		}else if(InputResultType().ExpWidth()>0){
			vhdl<<tab<<declare("result_exp", OutputResultType().ExpWidth(), true) <<" <= input_result_exp;\n";
		}
		if(!splitMultiplier){
			vhdl<<tab<<declare("result_frac", partialFracWidth, true) <<" <= stage_result_frac * input_result_frac;\n";
			nextCycle();
		}else{
			/*	0.1aaaa * 0.100bb
				0.0100000000 + 0.00aaaa0000 + 0.0000bb0000 + 0.0000rrrrrr
			
			*/
			
			declare("result_frac_r", partialFracWidth, true);
			vhdl<<tab<<"result_frac_r(FW_A+FW_B+1) <= '0';\n";
			vhdl<<tab<<"result_frac_r(FW_A+FW_B) <= '1';\n";
			vhdl<<tab<<"result_frac_r(FW_A+FW_B-1 downto AW_A+AW_B) <= (others=>'0');\n";
			vhdl<<tab<<"result_frac_r(AW_A+AW_B-1 downto 0) <= input_result_frac(AW_A-1 downto 0) * stage_result_frac(AW_B-1 downto 0);\n";
			
			declare("result_frac_a", partialFracWidth, true);
			vhdl<<tab<<"result_frac_a(FW_A+FW_B+1 downto AW_A+FW_B) <= (others=>'0');\n";
			vhdl<<tab<<"result_frac_a(AW_A+FW_B-1 downto FW_B) <= input_result_frac(AW_A-1 downto 0);\n";
			vhdl<<tab<<"result_frac_a(FW_B-1 downto 0) <= (others=>'0');\n";
			
			declare("result_frac_b", partialFracWidth, true);
			vhdl<<tab<<"result_frac_b(FW_B+FW_A+1 downto AW_B+FW_A) <= (others=>'0');\n";
			vhdl<<tab<<"result_frac_b(AW_B+FW_A-1 downto FW_A) <= stage_result_frac(AW_B-1 downto 0);\n";
			vhdl<<tab<<"result_frac_b(FW_A-1 downto 0) <= (others=>'0');\n";
			
			declare("result_frac_ab", partialFracWidth, true);
			vhdl<<tab<<"result_frac_ab <= result_frac_a+result_frac_b;\n";
			nextCycle();
			nextCycle();
			vhdl<<tab<<declare("result_frac", partialFracWidth, true)<<" <= result_frac_r + result_frac_ab;\n";
			nextCycle();
		}
		
		vhdl<<tab<<declare("result_frac_norm", partialFracWidth, true)<<" <= result_frac when result_frac("<<partialFracWidth<<"-1)='1' else result_frac("<<partialFracWidth-2<<" downto 0)&'0';\n";
		if((impl->OutputResultType().ExpWidth()>0) || (InputResultType().ExpWidth()>0)){
			vhdl<<tab<<declare("result_exp_norm", OutputResultType().ExpWidth(), true)<<" <= result_exp+1 when result_frac("<<partialFracWidth<<"-1)='1' else result_exp;\n";
		}else{
			vhdl<<tab<<declare("result_exp_norm", OutputResultType().ExpWidth(), true)<<" <= \"1\" when result_frac("<<partialFracWidth<<"-1)='1' else \"0\";\n";
		}
		nextCycle();
		
		declare("result_frac_rounded", OutputResultType().FracWidth(), true);
		declare("result_exp_rounded", OutputResultType().ExpWidth(), true);
		if(OutputResultType().FracWidth()>partialFracWidth-1)
			throw std::invalid_argument("MultiplierExpStage won't pad.");
		
		// Bleh, truncation FTW!
		vhdl<<tab<<"result_frac_rounded <= result_frac_norm("<<partialFracWidth-2<<" downto "<<partialFracWidth-1-OutputResultType().FracWidth()<<");\n";
		vhdl<<tab<<"result_exp_rounded <= result_exp_norm;\n";
		
		nextCycle();
		
		syncCycleFromSignal("result_exp_rounded");
		syncCycleFromSignal("result_frac_rounded");
		if(OutputResidualType().Width()!=0){
			syncCycleFromSignal("stage_residual");
			vhdl<<tab<<"oResidual <= stage_residual;\n";	
		}
		vhdl<<tab<<"oResult <= result_exp_rounded & result_frac_rounded;\n";
	}
	
	virtual ~MultiplierV2ExpStage()
	{}
	
	virtual residual_type<T> OutputResidualType() const
	{ return m_outputResidualType; }
	
	virtual result_type<T> OutputResultType() const
	{ return m_outputResultType; }
	
	virtual std::pair<T,T> Execute(const T &residual, const T &result) const
	{
		std::pair<T,T> stage=m_impl->Execute(residual, 1.0);
		
		return std::make_pair(stage.first, m_outputResultType.Floor(result*stage.second));
	}
};

}; // random
}; // flopoco

#endif
