#ifndef flopoco_random_fixed_point_exp_stage_hpp
#define flopoco_random_fixed_point_exp_stage_hpp

#include "Operator.hpp"
#include "utils.hpp"

#include "residual_type.hpp"
#include "result_type.hpp"

#include <boost/smart_ptr.hpp>
#include <boost/lexical_cast.hpp>

namespace flopoco
{
	
extern vector<Operator *> oplist;
	
namespace random
{

class FixedPointExpStage
	: public flopoco::Operator
{
public:
	// Lazy templatisation
	typedef double T;
protected:


	T m_expMu, m_expSigma;
	residual_type<T> m_inputResidualType;
	result_type<T> m_inputResultType;

	FixedPointExpStage(Target* target,
		const std::map<std::basic_string<char>, double> &delayMap,
		T expMu, T expSigma,
		residual_type<T> inputResidualType,
		result_type<T> inputResultType
	)
		: flopoco::Operator(target, delayMap)
		, m_expMu(expMu)
		, m_expSigma(expSigma)
		, m_inputResidualType(inputResidualType)
		, m_inputResultType(inputResultType)
	{}
	
	void AddOperator(Operator *op)
	{
		oplist.push_back(op);
	}
	
	std::string shortenId(std::string x)
	{
		if(x.size()<=32)
			return x;
		
		uint32_t acc=0;
		for(unsigned i=0;i<x.size();i++){
			acc=acc*31+x[i];
		}
		return x.substr(0, 32-10)+"_"+boost::lexical_cast<std::string>(acc);
	}
public:
    virtual ~FixedPointExpStage()
	{};
	
	T ReferenceExp(const T &x) const
	{ return exp(m_expMu+m_expSigma*x); }
	
	T Mu() const
	{ return m_expMu; }
	
	T Sigma() const
	{ return m_expSigma; }
	
	const residual_type<T> &InputResidualType() const
	{ return m_inputResidualType; }
	
	result_type<T> InputResultType() const
	{ return m_inputResultType; }
	
	virtual residual_type<T> OutputResidualType() const =0;
	virtual result_type<T> OutputResultType() const=0;
	
	virtual std::pair<T,T> Execute(const T &residual, const T &result) const=0;

	virtual void emulate(TestCase * tc) {
		T iResidual = m_inputResidualType.FromMpz(tc->getInputValue("iResidual"));
		T iResult=1.0;
		if(m_inputResultType.Width()>0){
			iResult = m_inputResultType.FromMpz(tc->getInputValue("iResult"));
		}

		std::pair<T,T> output=Execute(iResidual, iResult);
		
		std::stringstream tmp;
		tmp<<"("<<iResidual<<","<<iResult<<") -> ("<<output.first<<","<<output.second<<")";
		tc->addComment(tmp.str());
		std::cerr<<tmp.str()<<"\n";
		
		if(OutputResidualType().Width()>0){
			tc->addExpectedOutput("oResidual", OutputResidualType().ToMpz(output.first));
		}
		std::cerr<<" Setting oResult, type="<<OutputResultType()<<"\n";
		tc->addExpectedOutput("oResult", OutputResultType().ToMpz(output.second));
		std::cerr<<"Done\n";
	}
	
	virtual void buildStandardTestCases(TestCaseList * tcl) {
		TestCase *tc=new TestCase(this);

		tc->addInput("iResidual", m_inputResidualType.ToMpz(m_inputResidualType.RangeMin()));
		if(m_inputResultType.Width()>0){
			tc->addInput("iResult", m_inputResultType.ToMpz(m_inputResultType.RangeMin()));
		}
		emulate(tc);
		tcl->add(tc);
		
		tc=new TestCase(this);
		tc->addInput("iResidual", m_inputResidualType.ToMpz(m_inputResidualType.RangeMax()));
		if(m_inputResultType.Width()>0){
			tc->addInput("iResult", m_inputResultType.ToMpz(m_inputResultType.RangeMin()));
		}
		emulate(tc);
		tcl->add(tc);
		
		tc=new TestCase(this);
		tc->addInput("iResidual", m_inputResidualType.ToMpz(m_inputResidualType.RangeMin()));
		if(m_inputResultType.Width()>0){
			tc->addInput("iResult", m_inputResultType.ToMpz(m_inputResultType.RangeMax()));
		}
		emulate(tc);
		tcl->add(tc);
		
		tc=new TestCase(this);
		tc->addInput("iResidual", m_inputResidualType.ToMpz(m_inputResidualType.RangeMax()));
		if(m_inputResultType.Width()>0){
			tc->addInput("iResult", m_inputResultType.ToMpz(m_inputResultType.RangeMax()));
		}
		emulate(tc);
		tcl->add(tc);
	}
	
	virtual TestCase* buildRandomTestCase(int i){
		TestCase *tc;
		tc = new TestCase(this); 
		T iResidual=m_inputResidualType.RandomElement();
		T iResult=m_inputResultType.RandomElement();
		
		
		tc->addInput("iResidual", m_inputResidualType.ToMpz(iResidual));
		if(m_inputResultType.Width()>0){
			tc->addInput("iResult", m_inputResultType.ToMpz(iResult));
		}
		emulate(tc);
		return tc;
	}
};

typedef boost::shared_ptr<FixedPointExpStage> FixedPointExpStagePtr;

}; // random
}; // flopoco

#endif
