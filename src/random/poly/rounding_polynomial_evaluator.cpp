
#include "fixed_point_polynomial_evaluator.hpp"

namespace flopoco
{
namespace random
{

class RoundingPolynomialEvaluator
	: public FixedPointPolynomialEvaluator
{
private:
	fixed_format_t m_inputFormat;
	std::vector<fixed_format_t> m_coefficientFormats;
	
	fixed_format_t m_outputFormat;
	int m_guardBits;
public:
	
	RoundingPolynomialEvaluator(
			const fixed_format_t &inputFormat,
			const std::vector<fixed_format_t> &coefficientFormats,
			int outputLsb,
			int guardBits,
			Target *target
		)
		: FixedPointPolynomialEvaluator(target)
		, m_inputFormat(inputFormat)
		, m_coefficientFormats(coefficientFormats)
		, m_guardBits(guardBits)
	{
		if(guardBits<=0)
			throw std::string("Must have non-negative guard bits.");
		
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
			nextCycle();
			
			std::string stageName=join("add_", i);
			addInput(join("a",i), coefficientFormats[i].width());
			fixed_format_t stageType=AddStatement(stageName, join("a",i), coefficientFormats[i], mulName, mulType);
			nextCycle();
			
			int roundLsb=outputLsb-guardBits-i;
			
			if(i==0 || (stageType.lsb>=roundLsb)){
				accName=stageName;
				accType=stageType;
			}else if(stageType.lsb < roundLsb){
				accName=join("round_",i);
				accType=stageType;
				accType.lsb=roundLsb;
				vhdl<<declare(accName,accType.width())<<"<="<<RoundExpr(accType,stageName,stageType)<<";\n";
			}
		      
			
		}
		
		fixed_format_t outputType=accType;
		
		if(outputType.msb <= outputLsb)
			throw std::string("RoundingPolynomialEvaluator - outputType.msb <= outputLsb.");
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

};

FixedPointPolynomialEvaluator *CreateRoundingPolynomialEvaluator(
	const fixed_format_t &inputFormat,
	const std::vector<fixed_format_t> &coefficientFormats,
	int outputLsb,
	int guardBits,
	const mpfr::mpreal &errorBudget,
	Target *target
){
	if(errorBudget < 0){
		std::stringstream acc;
		acc<<"CreateFullPrecisionPolynomialEvaluator - Error budget less than zero makes no sense.";
		throw std::string(acc.str());
	}
	
	return new RoundingPolynomialEvaluator(inputFormat, coefficientFormats, outputLsb, guardBits, target);
}

}; // random
}; // flopoco
