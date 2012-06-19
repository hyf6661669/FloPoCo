#ifndef flopoco_random_fixed_point_exp_tester_hpp
#define flopoco_random_fixed_point_exp_tester_hpp

#include "fixed_point_exp_stage.hpp"


namespace flopoco
{
namespace random
{

class FixedPointExpTester
	: public FixedPointExpStage
{
protected:
	FixedPointExpStagePtr m_impl;
public:
	FixedPointExpTester(
		Target *target,
		FixedPointExpStagePtr impl
	)
		: FixedPointExpStage(target, emptyDelayMap, impl->Mu(), impl->Sigma(), impl->InputResidualType(), impl->InputResultType())
		, m_impl(impl)
	{
		if(InputResultType().Width()>0)
			throw std::invalid_argument("FixedPointExpTester - result in must be one.");
		if(OutputResidualType().Width()>0)
			throw std::invalid_argument("FixedPointExpTester - output residual must be empty.");
		
		ostringstream name;
		name << "Test_"<<impl->getName();
		setName(name.str());
		setCopyrightString("Imperial College 2012");

		addInput ("iResidual" , InputResidualType().Width(), true);
		addOutput ("oResult" , OutputResultType().Width(), 2, true);
		inPortMap(impl.get(), "iResidual", "iResidual");
		outPortMap(impl.get(), "oResult", "result");
	
		vhdl << tab << instance(impl.get(), "impl");
		syncCycleFromSignal("result");
		vhdl << tab << "oResult <= result;\n";
	}
	
	virtual ~FixedPointExpTester()
	{}
	
	virtual residual_type<T> OutputResidualType() const
	{ return m_impl->OutputResidualType(); }
	
	virtual result_type<T> OutputResultType() const
	{ return m_impl->OutputResultType(); }
	
	virtual std::pair<T,T> Execute(const T &residual, const T &result) const
	{
		throw std::logic_error("Execute doesn't work on FixedPointExpTester.");
	}
	
	virtual void emulate(TestCase * tc) {
		T iResidual = m_inputResidualType.FromMpz(tc->getInputValue("iResidual"));
		
		T oResult=ReferenceExp(iResidual);
		
		tc->addExpectedOutput("oResult", OutputResultType().ToMpz(OutputResultType().Ceil(oResult)));
		tc->addExpectedOutput("oResult", OutputResultType().ToMpz(OutputResultType().Floor(oResult)));
	}
};

}; // random
}; // flopoco

#endif
