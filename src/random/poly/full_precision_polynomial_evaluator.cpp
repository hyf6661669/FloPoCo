
#include "fixed_point_polynomial_evaluator.hpp"

namespace flopoco
{
namespace random
{

class FullPrecisionPolynomialEvaluator
	: public FixedPointPolynomialEvaluator
{
private:
	fixed_format_t m_inputFormat;
	std::vector<fixed_format_t> m_coefficientFormats;
	
	fixed_format_t m_outputFormat;
public:
	
	FullPrecisionPolynomialEvaluator(
			const fixed_format_t &inputFormat,
			const std::vector<fixed_format_t> &coefficientFormats,
			int outputLsb,
			Target *target
		)
		: FixedPointPolynomialEvaluator(target)
		, m_inputFormat(inputFormat)
		, m_coefficientFormats(coefficientFormats)
	{
		if(inputFormat.width()<=0)
			throw std::string("Input format must have positive width.");
		for(int i=0;i<(int)coefficientFormats.size();i++){
			if(coefficientFormats[i].width()<=1)
				throw std::string("Coefficient formats must have width greater than 1 (sanity check, might be situations where this isn't true).");
		}
		
		int degree=coefficientFormats.size()-1;
		
		std::stringstream name;
		name<<"FullPrecisionPolynomialEvaluator_d"<<degree<<"_uid"<<getNewUId();
		setName(name.str());
		
		useNumericStd();
		
		std::string inputName="Y";
		fixed_format_t inputType=inputFormat;
		addInput(inputName, inputType.width());
		
		fixed_format_t accType=coefficientFormats[degree];
		addInput(join("a",degree),accType.width());
		std::string accName=join("a",degree);
		
		for(int i=degree-1;i>=0;i--){
			std::string mulName = join("mul_",i);
			fixed_format_t mulType=MultiplyStatement(mulName, accName, accType, inputName, inputType);
			nextCycle();
			
			std::string stageName=join("acc_", i);
			addInput(join("a",i), coefficientFormats[i].width());
			fixed_format_t stageType=AddStatement(stageName, join("a",i), coefficientFormats[i], mulName, mulType);
			nextCycle();
			
			accName=stageName;
			accType=stageType;
		}
		
		fixed_format_t outputType=accType;
		
		if(outputType.msb <= outputLsb)
			throw std::string("FullPrecisionPolynomialEvaluator - outputType.msb <= outputLsb.");
		outputType.lsb = outputLsb;
		std::string outputName="R";
		addOutput(outputName, outputType.width());
		vhdl<<outputName <<"<= "<<RoundExpr(outputType, accName, accType)<<";\n";
		
		m_outputFormat=outputType;
	}
	
	virtual int getPolynomialDegree() const
	{ return m_coefficientFormats.size()-1; }
	
	virtual fixed_format_t getInputFormat() const
	{ return m_inputFormat; }
	
	virtual fixed_format_t getCoefficientFormat(unsigned i) const
	{ return m_coefficientFormats.at(i); }
	
	virtual fixed_format_t getOutputFormat() const
	{ return m_outputFormat; }

	    virtual void emitHLSBody(HLSContext &ctxt, HLSScope &scope) const override
  {
    throw std::runtime_error("HLS is not supported.");
  }

};

FixedPointPolynomialEvaluator *CreateFullPrecisionPolynomialEvaluator(
	const fixed_format_t &inputFormat,
	const std::vector<fixed_format_t> &coefficientFormats,
	int outputLsb,
	const mpfr::mpreal &errorBudget,
	Target *target
){
	if(errorBudget < 0){
		std::stringstream acc;
		acc<<"CreateFullPrecisionPolynomialEvaluator - Error budget less than zero makes no sense.";
		throw std::string(acc.str());
	}
	
	return new FullPrecisionPolynomialEvaluator(inputFormat, coefficientFormats, outputLsb, target);
}

}; // random
}; // flopoco
