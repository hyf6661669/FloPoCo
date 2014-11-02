#ifndef flopoco_random_poly_approx_exp_stage_hpp
#define flopoco_random_poly_approx_exp_stage_hpp

#include "fixed_point_exp_stage.hpp"

#include "find_close_value.hpp"

#include "FixedPointFunctions/HOTBM.hpp"

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
		
		double scaleInput=pow(2.0, inputResidualType.MsbPower());
		double scaleOutput=pow(2.0, m_outputResultType.FracWidth()-m_outputResultType.FracWidthNonZero());
		
		std::cerr<<"in="<<m_inputResidualType<<", out="<<m_outputResultType<<"\n";
		std::cerr<<"scaleInput="<<scaleInput<<", scaleOutput="<<scaleOutput<<"\n";
		
		// I don't really like this, as I'm not holding the reference via shared_ptr... never mind
		ostringstream function;
		function <<"(exp("<<expMu<<" + x * 2^"<<inputResidualType.MsbPower()<<" * "<<expSigma<<")-1)*"<<scaleOutput;
		std::cerr<<" function = "<<function.str()<<"\n";
		flopoco::HOTBM *fe;
		
		fe = new flopoco::HOTBM(target, function.str(), "FAE", m_inputResidualType.Width(), m_outputResultType.FracWidthNonZero()+1, 1, false);
		oplist.push_back(fe);
		
		
		ostringstream name;
		name << "FuncApproxExpStage_"<<m_inputResidualType.DescriptionId();
		setName(shortenId(name.str()));
		setCopyrightString("Imperial College 2012");

		addInput ("iResidual" , InputResidualType().Width(), true);
		addOutput ("oResult" , OutputResultType().Width(), true);
		
		inPortMap(fe, "X", "iResidual");
		outPortMap(fe, "R", "result");
		
		vhdl<<instance(fe, "fe");
		syncCycleFromSignal("result");
		
		vhdl<<tab<<"oResult <= \"";
		for(int i=0;i<OutputResultType().Width()-OutputResultType().FracWidthNonZero();i++)
			vhdl<<"0";
		vhdl<<"\"  & result("<<OutputResultType().FracWidthNonZero()-1<<" downto 0);\n";
	}
	
	virtual ~FuncApproxExpStage()
	{}
	
	virtual residual_type<T> OutputResidualType() const
	{ return m_outputResidualType; }
	
	virtual result_type<T> OutputResultType() const
	{ return m_outputResultType; }
	
	virtual std::pair<T,T> Execute(const T &residual, const T &result) const
	{
		throw std::invalid_argument("Not implemented.");
		
		/*assert(result==1.0);
		uint64_t bits=m_inputResidualType.ToBits(residual);
		
		unsigned index=bits>>(m_inputResidualType.Width()-m_tableInputBits);
		T baseResidual=m_outputResidualType.FromBits(bits & ((1ull<<(m_inputResidualType.Width()-m_tableInputBits))-1));
		
		assert(m_tableIndexType.FromBits(index)+baseResidual==residual);
		
		return std::make_pair(baseResidual+m_table[index].residualOffset, m_table[index].roundedResult);*/
	}
};

}; // random
}; // flopoco

#endif
