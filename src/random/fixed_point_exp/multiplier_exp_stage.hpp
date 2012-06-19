#ifndef flopoco_random_multiplier_exp_stage_hpp
#define flopoco_random_multiplier_exp_stage_hpp

#include "fixed_point_exp_stage.hpp"

#include "Table.hpp"

/*	(1+ea)*(1+eb)
	// Output of result
	1+ea+eb+ea*eb
	// Then round with g guard bits, which can cause an error of 0.5 ulp up or down
	1+ea+eb+ea*eb + 2^(-wr-g-1)
	
	-2^(-wr-g-1) - min(ea+eb+ea*eb) .. max(ea+eb+ea*eb) + 2^(-wr-g-1)
	
	-2^(-wr-1) < -2^(-wr-g-1) - min(ea+eb+ea*eb)
	-2^(-wr-1) + 2^(-wr-g-1) < -min(ea+eb+ea*eb)
	2^(-wr-1) -2^(-wr-g-1) > min(ea+eb+ea*eb)
	
	g==0 -> 0 > min(ea+eb+ea*eb)
	g==1 -> 2^(-wr-1) - 2^(-wr-2) = 2^(-wr-2) > ...
	g==2 -> 2^(-
	
	0*2^(-wr-1)
	0.5*2^(-wr-1)
	0.75*2^(-wr-1)
*/

namespace flopoco
{
namespace random
{

class MultiplierExpStage
	: public FixedPointExpStage
{
protected:
	FixedPointExpStagePtr m_impl;

	residual_type<T> m_outputResidualType;
	result_type<T> m_outputResultType;

	std::string zeros(int n)
	{
		std::stringstream ss;
		ss<<"\"";
		for(int i=0;i<n;i++){
			ss<<"0";
		}
		ss<<"\"";
		return ss.str();
	}
	
public:
	MultiplierExpStage(
		Target *target,
		FixedPointExpStagePtr impl,
		result_type<T> inputResultType,
		unsigned outputFracBits,
		map<string, double> inputDelays = emptyDelayMap
	)
		: FixedPointExpStage(target, inputDelays, impl->Mu(), impl->Sigma(), impl->InputResidualType(), inputResultType)
		, m_impl(impl)
		, m_outputResidualType(impl->OutputResidualType())
		, m_outputResultType(outputFracBits)
	{
		std::cerr<<"Multiplier src type : "<<inputResultType<<"\n";
		std::cerr<<"Stage res type : "<<impl->OutputResultType()<<"\n";
		std::cerr<<"range = "<<inputResultType.RangeMin()*impl->OutputResultType().RangeMin()<<", "<<inputResultType.RangeMax()*impl->OutputResultType().RangeMax()<<"\n";
		m_outputResultType.Add(m_outputResultType.Floor(inputResultType.RangeMin()*impl->OutputResultType().RangeMin()));
		m_outputResultType.Add(m_outputResultType.Floor(inputResultType.RangeMax()*impl->OutputResultType().RangeMax()));
		std::cerr<<"Multiplier output type : "<<m_outputResultType<<"\n";
		
		ostringstream name;
		name << "MultiplierExpStage_r"<<outputFracBits<<"_"<<impl->getName();
		setName(name.str());
		setCopyrightString("Imperial College 2012");

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
		declare("stage_result_frac", impl->OutputResultType().FracWidth()+1, true);
		if(impl->OutputResultType().FracWidthNonZero()+1<impl->OutputResultType().FracWidth()){
			// This didn't convince XST to do anything exciting :(
			vhdl<<tab<<"stage_result_frac <= \"1";
			for(int i=0;i<impl->OutputResultType().FracWidth()-impl->OutputResultType().FracWidthNonZero();i++){
				vhdl<<"0";
			}
			vhdl<<"\" & stage_result("<<impl->OutputResultType().FracWidthNonZero()-1<<" downto 0);\n";
		}else{
			vhdl<<tab<<"stage_result_frac <= '1' & stage_result("<<impl->OutputResultType().FracWidth()-1<<" downto 0);\n";
		}
		
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
		
		if(impl->OutputResultType().FracWidthNonZero()+1<impl->OutputResultType().FracWidth()){
			int WA=impl->OutputResultType().FracWidth(), WB=InputResultType().FracWidth();
			int NZA=impl->OutputResultType().FracWidthNonZero();
			int WR=WA+WB+2;
			// 100AA * 1BBBB ->  0000RRRRRR + 00BBBB0000 + 0000AA0000 + 0100000000
			vhdl<<tab<<declare("result_frac_prod_amb", NZA+WB, true) << "<= stage_result_frac("<<NZA-1<<" downto 0) * input_result_frac("<<WB-1<<" downto 0);\n";
			vhdl<<tab<<declare("result_frac_prod_apb", WR, true) << " <= \n"
				<< tab<<tab<<"("<<zeros(WR-NZA-WB)<< "& stage_result_frac("<<NZA-1<<" downto 0) & "<<zeros(WB)<<")"
				<< tab<<tab<<" + \n"
				<< tab<<tab<<"("<<zeros(WR-WB-WA)<< "& input_result_frac("<<WB-1<<" downto 0) & "<<zeros(WA)<<");\n";
			vhdl<<tab<<declare("result_frac_prod", WR, true) << " <= \n"
				<< tab<<tab<<"(\"01\"&"<<zeros(WR-NZA-WB-2)<<"&result_frac_prod_amb)\n"
				<< tab<<tab<<" + \n" 
				<< tab<<tab<<"result_frac_prod_apb;\n";
		}else{
			vhdl<<tab<<declare("result_frac_prod", partialFracWidth, true) <<" <= stage_result_frac * input_result_frac;\n";
		}
		nextCycle();
		
		vhdl<<tab<<declare("result_frac_norm", partialFracWidth, true)<<" <= result_frac_prod when result_frac_prod("<<partialFracWidth<<"-1)='1' else result_frac_prod("<<partialFracWidth-2<<" downto 0)&'0';\n";
		if((impl->OutputResultType().ExpWidth()>0) || (InputResultType().ExpWidth()>0)){
			vhdl<<tab<<declare("result_exp_norm", OutputResultType().ExpWidth(), true)<<" <= result_exp+1 when result_frac_prod("<<partialFracWidth<<"-1)='1' else result_exp;\n";
		}else{
			vhdl<<tab<<declare("result_exp_norm", OutputResultType().ExpWidth(), true)<<" <= \"1\" when result_frac_prod("<<partialFracWidth<<"-1)='1' else \"0\";\n";
		}
		nextCycle();
		
		declare("result_frac_rounded", OutputResultType().FracWidth(), true);
		declare("result_exp_rounded", OutputResultType().ExpWidth(), true);
		if(OutputResultType().FracWidth()>partialFracWidth-1)
			throw std::invalid_argument("MultiplierExpStage won't pad.");
		
		if(OutputResultType().FracWidth()==partialFracWidth-1){
			vhdl<<tab<<"result_frac_rounded <= result_frac_norm("<<partialFracWidth-2<<" downto "<<partialFracWidth-1-OutputResultType().FracWidth()<<");\n";
			vhdl<<tab<<"result_exp_rounded <= result_exp_norm;\n";
		}else{
			vhdl<<tab<<declare("result_frac_rounded_pre", partialFracWidth, true)<<" <= (result_frac_norm + \"";
			for(int i=0;i<OutputResultType().FracWidth()+2;i++)
				vhdl<<"0";
			for(int i=OutputResultType().FracWidth()+2;i<partialFracWidth;i++)
				vhdl<<"1";
			vhdl<<"\" );\n";
			vhdl<<tab<<"result_frac_rounded <= result_frac_rounded_pre("<<partialFracWidth-2<<" downto "<<partialFracWidth-1-OutputResultType().FracWidth()<<");\n";
			vhdl<<tab<<"result_exp_rounded <= result_exp_norm;\n";
		}
		
		nextCycle();
		
		syncCycleFromSignal("result_exp_rounded");
		syncCycleFromSignal("result_frac_rounded");
		if(OutputResidualType().Width()!=0){
			syncCycleFromSignal("stage_residual");
			vhdl<<tab<<"oResidual <= stage_residual;\n";	
		}
		vhdl<<tab<<"oResult <= result_exp_rounded & result_frac_rounded;\n";
	}
	
	virtual ~MultiplierExpStage()
	{}
	
	virtual residual_type<T> OutputResidualType() const
	{ return m_outputResidualType; }
	
	virtual result_type<T> OutputResultType() const
	{ return m_outputResultType; }
	
	virtual std::pair<T,T> Execute(const T &residual, const T &result) const
	{
		std::pair<T,T> stage=m_impl->Execute(residual, 1.0);
		
		return std::make_pair(stage.first, m_outputResultType.RoundHalfDown(result*stage.second));
	}
};

}; // random
}; // flopoco

#endif
