#ifndef flopoco_random_func_approx_exp_stage_hpp
#define flopoco_random_func_approx_exp_stage_hpp

#include "fixed_point_exp_stage.hpp"

#include "find_close_value.hpp"

#include "FixFunctions/HOTBM.hpp"
#include "FixFunctions/HOTBM/Minimax.hh"
#include "FixFunctions/HOTBM/Exhaustive.hh"

#include "sollya.h"

namespace flopoco
{
namespace random
{

class FuncApproxExpStage
	: public FixedPointExpStage
{
protected:
	residual_type<T> m_outputResidualType;
	result_type<T> m_outputResultType;

	boost::shared_ptr<flopoco::HOTBM> m_fe;

public:
	FuncApproxExpStage(
		Target *target,
		T expMu,
		T expSigma,
		residual_type<T> inputResidualType,
		int resultFracWidth, // total result fracWidth. Output will be faithfully rounded to this
		map<string, double> inputDelays = emptyDelayMap
	)
		: FixedPointExpStage(target, inputDelays, expMu, expSigma, inputResidualType, result_type<T>(0))
		, m_outputResidualType(inputResidualType.drop_msbs(inputResidualType.Width()))
		, m_outputResultType(resultFracWidth)
	{
		if(inputResidualType.HasSign())
			throw std::invalid_argument("FuncApproxExpState - Sign not currently supported.");
		if(inputResidualType.MsbPower()>-1)
			throw std::invalid_argument("FuncApproxExpState - MSB power must be less than -1 (i.e. range of [0,0.5) or less).");
		
		m_outputResultType.Add(m_outputResultType.Floor( ReferenceExp(inputResidualType.RangeMin())));
		m_outputResultType.Add(m_outputResultType.Ceil( ReferenceExp(inputResidualType.RangeMax())));
		
		if(m_outputResultType.RangeMin()<1 || 1.5<m_outputResultType.RangeMax())
			throw std::invalid_argument("FuncApproxExpState - Result must be in range [1,1.5)");
		
		double scaleInput=pow(2.0, inputResidualType.MsbPower()+1);
		double scaleOutput=pow(2.0, m_outputResultType.FracWidth()-m_outputResultType.FracWidthNonZero()-1);
		
		std::cerr<<"in="<<m_inputResidualType<<", out="<<m_outputResultType<<"\n";
		std::cerr<<"scaleInput="<<scaleInput<<", scaleOutput="<<scaleOutput<<"\n";
		
		ostringstream function;
		function <<"(exp("<<expMu<<" + x * 2^"<<inputResidualType.MsbPower()+1<<" * "<<expSigma<<")-1)*"<<scaleOutput;
		std::cerr<<" function = "<<function.str()<<"\n";
		
		REPORT(DEBUG, "HOTBM source function = "<<function.str());
		
		
		int pp=getToolPrecision();
		setToolPrecision(96);
		int pInfPoints=Minimax::setInfNormPoints(31);	// We know it is very smooth
		Exhaustive::ScoreType pST=Exhaustive::setScoreType(Exhaustive::ScoreArea);
		int pMinAlpha=Exhaustive::setMinAlpha(std::min(inputResidualType.Width()-2, target->lutInputs()));
		int pMaxAlpha=Exhaustive::setMaxAlpha(target->lutInputs()+2);
		
		m_fe.reset(new flopoco::HOTBM(target, function.str(), "FAE", m_inputResidualType.Width(), m_outputResultType.FracWidthNonZero()+1, 2, false));
		oplist.push_back(m_fe.get());
		
		setToolPrecision(pp);
		Minimax::setInfNormPoints(pInfPoints);
		Exhaustive::setScoreType(pST);
		Exhaustive::setMinAlpha(pMinAlpha);
		Exhaustive::setMaxAlpha(pMaxAlpha);
		
		ostringstream name;
		name << "FuncApproxExpStage_"<<m_inputResidualType.DescriptionId();
		setName(shortenId(name.str()));
		setCopyrightString("Imperial College 2012");

		addInput ("iResidual" , InputResidualType().Width(), true);
		addOutput ("oResult" , OutputResultType().Width(), true);
		
		inPortMap(m_fe.get(), "X", "iResidual");
		outPortMap(m_fe.get(), "R", "result");
		
		vhdl<<instance(m_fe.get(), "fe");
		syncCycleFromSignal("result");
		
		vhdl<<tab<<"oResult <= \"";
		for(unsigned i=0;i<OutputResultType().Width()-OutputResultType().FracWidthNonZero();i++)
			vhdl<<"0";
		vhdl<<"\"  & result("<<OutputResultType().FracWidthNonZero()-1<<" downto 0);\n";
	}
	
	virtual ~FuncApproxExpStage()
	{
		// TODO : should I take m_fe out of oplist when destructing?
	}
	
	virtual residual_type<T> OutputResidualType() const
	{ return m_outputResidualType; }
	
	virtual result_type<T> OutputResultType() const
	{ return m_outputResultType; }
	
	virtual std::pair<T,T> Execute(const T &residual, const T &result) const
	{
		assert(result==1.0);
		mpz_class bits=m_inputResidualType.ToMpz(residual);
		
		boost::shared_ptr<flopoco::TestCase> tc(new flopoco::TestCase(m_fe.get()));
		
		tc->setInputValue("X", bits); // Unsafe
		m_fe->emulate(tc.get());
		const std::vector<mpz_class> &res=tc->getExpectedOutputValues("R");
		
		if((res.size()!=1) && (res.size()!=2))
			throw std::logic_error("FuncApproxExpStage::Execute - Expected one or two values from HOTBM.");
		
		T exact=ReferenceExp(residual);
		
		if(res[0]<0)
			throw std::logic_error("FuncApproxExpStage::Execute - HOTBM returned a negative result.");
		T a=m_outputResultType.FromMpz(res[0]);
		T b=a;
		if(res.size()==2){
			if(res[1]<0)
				throw std::logic_error("FuncApproxExpStage::Execute - HOTBM returned a negative result.");
			b=m_outputResultType.FromMpz(res[1]);
		}
		
		if(a>b)
			std::swap(a,b);
		
		REPORT(DEBUG, "input="<<residual<<", exact="<<exact<<", got=("<<a<<","<<b);
		
		if((exact < a) || (b<exact)){
			throw std::logic_error("FuncApproxExpStage::Execute - Output of HOTBM does not bracket exact result (not faithful).");
		}
		
		if(a!=b){
			if(m_outputResultType.Next(a)!=b)
				throw std::logic_error("FuncApproxExpStage::Execute - Output of HOTBM does not differ by one bit (not faithful).");
		}
		
		if( (exact-a) > (b-exact) ){
			return std::make_pair((T)0.0, a);
		}else{
			return std::make_pair((T)0.0, b);
		}
	}
};

}; // random
}; // flopoco

#endif
