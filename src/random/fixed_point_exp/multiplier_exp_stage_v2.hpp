#ifndef flopoco_random_multiplier_exp_stage_v2_hpp
#define flopoco_random_multiplier_exp_stage_v2_hpp

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

class MultiplierExpStageV2
	: public FixedPointExpStage
{
protected:
	FixedPointExpStagePtr m_hi, m_lo;

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
	MultiplierExpStageV2(
		Target *target,
		FixedPointExpStagePtr hi,
		FixedPointExpStagePtr lo,
		unsigned outputFracBits,
		map<string, double> inputDelays = emptyDelayMap
	)
		: FixedPointExpStage(target, inputDelays, hi->Mu()+lo->Mu(), hi->Sigma(), hi->InputResidualType(), hi->InputResultType())
		, m_hi(hi)
		, m_lo(lo)
		, m_outputResidualType(impl->OutputResidualType())
		, m_outputResultType(outputFracBits)
	{
		std::cerr<<"MultiplierV2 src type : "<<inputResultType<<"\n";
		std::cerr<<"  hi res type : "<<hi->OutputResultType()<<"\n";
		std::cerr<<"  lo res type : "<<lo->OutputResultType()<<"\n";
		std::cerr<<"  range = "<<hi->OutputResultType.RangeMin()*lo->OutputResultType().RangeMin()<<", "<<lo->OutputResultType.RangeMax()*lo->OutputResultType().RangeMax()<<"\n";
		m_outputResultType.Add(m_outputResultType.Floor(hi->OutputResultType.RangeMin()*lo->OutputResultType().RangeMin()));
		m_outputResultType.Add(m_outputResultType.Floor(hi->OutputResultType.RangeMax()*lo->OutputResultType().RangeMax()));
		std::cerr<<"Multiplier output type : "<<m_outputResultType<<"\n";
		
		bool needsPostNorm=true;
		
		ostringstream name;
		name << "MultiplierExpStage_r"<<outputFracBits<<"_"<<hi->getName()<<"_"<<lo->getName();
		setName(shortenId(name.str()));
		setCopyrightString("Imperial College 2012");

		addInput ("iResidual" , hi->InputResidualType().Width(), true);
		if(OutputResidualType().Width()!=0){
			addOutput ("oResidual" , OutputResidualType().Width(), 1, true);
		}
		addOutput ("oResult" , OutputResultType().Width(), 1, true);
		
		
		inPortMap(hi.get(), "iResidual", "iResidual");
		assert(hi->InputResultType().Width()==0);
		outPortMap(hi.get(), "oResidual", "hi_residual");
		outPortMap(hi.get(), "oResult", "hi_result");
		vhdl<<tab<<instance(hi.get(), "hi_stage");
		
		inPortMap(hi.get(), "iResidual", "hi_residual");
		assert(lo->InputResultType().Width()==0);
		if(lo->OutputResidualType().Width()!=0){
			outPortMap(hi.get(), "oResidual", "lo_residual");
		}
		outPortMap(hi.get(), "oResult", "lo_result");
		vhdl<<tab<<instance(lo.get(), "lo_stage");
		
		if(lo->OutputResidualType().Width()!=0){
			syncCycleFromSignal("lo_residual");
			vhdl<<tab<<"oResidual <= lo_residual;\n";
		}

		syncCycleFromSignal("hi_result");
		syncCycleFromSignal("lo_result");
		
		if(hi->OutputResultType().ExpWidth()>0){
			vhdl<<tab<<declare("hi_result_exp", hi->OutputResultType().ExpWidth(), true)<<" <= hi_result("<<hi->OutputResultType().Width()-1<<" downto "<<hi->OutputResultType().FracWidth()<<");\n";
		}
		declare("hi_result_frac", hi->OutputResultType().FracWidth()+1, true);
		vhdl<<tab<<"hi_result_frac <= '1' & stage_result("<<hi->OutputResultType().FracWidth()-1<<" downto 0);\n";
		
		if(lo->InputResultType().ExpWidth()>0){
			vhdl<<tab<<declare("lo_result_exp", lo->OutputResultType().ExpWidth(), true)<<" <= lo_result("<<lo->OutputResultType().Width()-1<<" downto "<<lo->OutputResultType().FracWidth()<<");\n";
		}
		declare("lo_result_frac", lo->OutputResultType().FracWidth()+1, true)
		vhdl<<tab<<<<"lo_result_frac <= '1' & input_result("<<lo->OutputResultType().FracWidth()-1<<" downto 0);\n";
		
		unsigned partialFracWidth=2+hi->OutputResultType().FracWidth()+lo->OutputResultType().FracWidth();
		if((hi->OutputResultType().ExpWidth()>0) && (lo->OutputResultType().ExpWidth()>0)){
			vhdl<<tab<<declare("result_exp", OutputResultType().ExpWidth(), true) <<" <= hi_result_exp + lo_result_exp;\n";
		}else if(hi->OutputResultType().ExpWidth()>0){
			vhdl<<tab<<declare("result_exp", OutputResultType().ExpWidth(), true) <<" <= hi_result_exp;\n";
		}else if(lo->OutputResultType().ExpWidth()>0){
			vhdl<<tab<<declare("result_exp", OutputResultType().ExpWidth(), true) <<" <= lo_result_exp;\n";
		}
		
		bool sparse_multiplier = lo->OutputResultType().FracWidthNonZero()+1 < lo->OutputResultType().FracWidth()
				|| hi->OutputResultType().FracWidthNonZero()+1 < hi->OutputResultType().FracWidth();
		
		if(sparse_multiplier){
			int WLO=lo->OutputResultType().FracWidth(), WHI=lo->OutputResultType().FracWidth();
			int NZLO=lo->OutputResultType().FracWidthNonZero(), NZHI=hi->OutputResultType().FracWidthNonZero();
			int WR=WHI+WLO+2;
			// 100AA * 10BBB ->  00000RRRRR + 00BBBB0000 + 0000AA0000 + 0100000000
			vhdl<<tab<<declare("result_frac_prod_amb", NZHI+NZLO, true) << "<= lo_result_frac("<<NZLO-1<<" downto 0) * hi_result_frac("<<NZHI-1<<" downto 0);\n";
			vhdl<<tab<<declare("result_frac_prod_apb", WR, true) << " <= \n"
				<< tab<<tab<<"("<<zeros(WR-NZLO-WHI)<< "& lo_result_frac("<<NZLO-1<<" downto 0) & "<<zeros(WHI)<<")"
				<< tab<<tab<<" + \n"
				<< tab<<tab<<"("<<zeros(WR-NZHI-WLO)<< "& hi_result_frac("<<NZHI-1<<" downto 0) & "<<zeros(WLO)<<");\n";
			nextCycle();
			vhdl<<tab<<declare("result_frac_prod", WR, true) << " <= \n"
				<< tab<<tab<<"(\"01\"&"<<zeros(WR-NZHI-NZLO-2)<<"&result_frac_prod_amb)\n"
				<< tab<<tab<<" + \n" 
				<< tab<<tab<<"result_frac_prod_apb;\n";
		}else{
			vhdl<<tab<<declare("result_frac_prod", partialFracWidth, true) <<" <= hi_result_frac * lo_result_frac;\n";
		}
		nextCycle();
		
		if(postNormUp && postNormDown){
			vhdl<<tab<<declare("result_frac_norm", partialFracWidth, true)<<" <= result_frac_prod when result_frac_prod("<<partialFracWidth<<"-1)='1' else result_frac_prod("<<partialFracWidth-2<<" downto 0)&'0';\n";
			if((impl->OutputResultType().ExpWidth()>0) || (InputResultType().ExpWidth()>0)){
				vhdl<<tab<<declare("result_exp_norm", OutputResultType().ExpWidth(), true)<<" <= result_exp+1 when result_frac_prod("<<partialFracWidth<<"-1)='1' else result_exp;\n";
			}else{
				vhdl<<tab<<declare("result_exp_norm", OutputResultType().ExpWidth(), true)<<" <= \"1\" when result_frac_prod("<<partialFracWidth<<"-1)='1' else \"0\";\n";
			}
			nextCycle();
		}else if(postNormUp){
			vhdl<<tab<<declare("result_frac_norm", partialFracWidth, true)<<" <= result_frac_prod;\n";
			vhdl<<tab<<declare("result_exp_norm", OutputResultType().ExpWidth(), true)<<" <= result_exp+1;\n";
			nextCycle(); // This should really be handled by fudging the exponent range of the output
		}else{
			vhdl<<tab<<declare("result_frac_norm", partialFracWidth, true)<<" <= result_frac_prod("<<partialFracWidth-2<<" downto 0)&'0';\n";
			vhdl<<tab<<declare("result_exp_norm", OutputResultType().ExpWidth(), true)<<" <= result_exp;\n";
		}
		
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
