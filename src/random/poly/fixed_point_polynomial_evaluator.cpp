#include "fixed_point_polynomial_evaluator.hpp"

#include <assert.h>
#include <sollya.h>
#include <UtilSollya.hh>

#include "IntMultiplier.hpp"

#include <boost/lexical_cast.hpp>

#include "random/utils/operator_factory.hpp"

namespace flopoco{
namespace random{

FixedPointPolynomialEvaluator::FixedPointPolynomialEvaluator(
	Target *target
)
	: Operator(target)
{}

std::string FixedPointPolynomialEvaluator::ExtendExpr(const fixed_format_t &outType, const std::string srcExpr, const fixed_format_t &srcType) const
{
	if(outType.msb<srcType.msb)
		throw std::string("ExtendType - Attempt to do narrowing conversion.");
	if(outType.msb==srcType.msb && (outType.isSigned!=srcType.isSigned))
		throw std::string("ExtendType - Attempt to do narrowing conversion.");
	
	std::stringstream acc;
	if(outType==srcType){
		acc<<srcExpr;
	}else if(!srcType.isSigned){
		// Original is unsigned
		acc<<"(";
		if(outType.msb>srcType.msb){
			acc<<"\"";
			for(int i=outType.msb;i>srcType.msb;i--){
				acc<<"0";
			}
			acc<<"\"&";
		}
		acc<<srcExpr;
		if(outType.lsb>srcType.msb)
			throw std::string("ExtendType - Attempt to do narrowing conversion.");
		if(outType.lsb<srcType.msb){
			acc<<"&\"";
			for(int i=outType.lsb;i<srcType.lsb;i++){
				acc<<"0";
			}
			acc<<"\"";
		}
		acc<<")";
	}else{
		// Original is signed
		if(!outType.isSigned)
			throw std::string("ExtendType - Attempt to do narrowing conversion.");
		
		acc<<"(";
		if(outType.msb>srcType.msb){
			acc<<"std_logic_vector(resize(signed(";
		}
		acc<<srcExpr;
		if(outType.lsb>srcType.msb)
			throw std::string("ExtendType - Attempt to do narrowing conversion.");
		if(outType.lsb<srcType.msb){
			acc<<"&\"";
			for(int i=outType.lsb;i<srcType.lsb;i++){
				acc<<"0";
			}
			acc<<"\"";
		}
		if(outType.msb>srcType.msb){
			acc<<"),"<<outType.width()<<"))";
		}
		acc<<")";
	}
	return acc.str();
}

std::string FixedPointPolynomialEvaluator::RoundExpr(const fixed_format_t &outType, const std::string srcExpr, const fixed_format_t &srcType)
{
	if(outType.msb!=srcType.msb || outType.isSigned !=srcType.isSigned)
		throw std::string("RoundExpr - msbs do not match.");
	
	if(outType.lsb < srcType.lsb)
		return ExtendExpr(outType, srcExpr, srcType);
	
	if(outType.lsb == srcType.lsb)
		return srcExpr;
	
	int toGo=srcType.width()-outType.width();

	manageCriticalPath(getTarget()->adderDelay(outType.width()));
	
	std::stringstream acc;
	acc<<"std_logic_vector((unsigned("<<srcExpr<<range(srcType.width()-1,toGo)<<")";
		acc<<" + ";
	acc<<"unsigned("<<srcExpr<<range(toGo-1,toGo-1)<<")))";

	return acc.str();
}

fixed_format_t FixedPointPolynomialEvaluator::MultiplyType(const fixed_format_t &aType, const fixed_format_t &bType) const
{
	fixed_format_t res={aType.isSigned||bType.isSigned, aType.msb+bType.msb+1, aType.lsb+bType.lsb};
	return res;
}

/*
	UU x UU 	: {2^2-1}*{2^2-1} -> {9} : UUUU
	SUU x UU 	: {-2^2,2^2-1} * {2^2-1} -> {-12,9} -> {-16+4,9} : SUUUU
	SUU x SUU	: {-2^2,2^2-1} * {-2^2,2^2-1} -> {-12,16} : SUUUUU
*/
fixed_format_t FixedPointPolynomialEvaluator::MultiplyStatement(std::string resName, std::string aName, const fixed_format_t &aType, std::string bName, const fixed_format_t &bType)
{
	if(bType.isSigned && !aType.isSigned)
		return MultiplyStatement(resName, bName, bType, aName, aType);
	
	fixed_format_t res=MultiplyType(aType,bType);

	vhdl << "--- "<<res<<" <- "<<aType<<" * "<<bType<<"\n";
	if(0){  
	    nextCycle();

	    vhdl << declare(resName, res.width()) << " <= std_logic_vector(\n";
	    if((!aType.isSigned) && (!bType.isSigned)){
		vhdl << "  unsigned("<<aName<<") * unsigned("<<bName<<")";
	    }else if(aType.isSigned && bType.isSigned){
		vhdl << "  signed("<<aName<<") * signed("<<bName<<")";
	    }else{
		assert(!bType.isSigned);
		vhdl << "  resize(signed("<<aName<<") * signed('0' & "<<bName<<"), "<<res.width()<<")";
	    }
	    vhdl<<");\n";

	    nextCycle();
	    nextCycle();
	}else{
	    std::map<std::string,double> inputDelays;
	    inputDelays["X"]=getCriticalPath();
	    inputDelays["Y"]=getCriticalPath();

	    IntMultiplier *mult=NULL;

	    if((!aType.isSigned) && (!bType.isSigned)){
		mult=new IntMultiplier(getTarget(),aType.width(),bType.width(),0,false, 0.7, inputDelays);
		oplist.push_back(mult);

		inPortMap(mult, "X", aName);
		inPortMap(mult, "Y", bName);
	    }else if(aType.isSigned && bType.isSigned){
		mult=new IntMultiplier(getTarget(),aType.width(),bType.width(),0,true, 0.7, inputDelays);
		oplist.push_back(mult);

		inPortMap(mult, "X", aName);
		inPortMap(mult, "Y", bName);

	    }else{
		assert(!bType.isSigned);
		mult=new IntMultiplier(getTarget(),aType.width(),bType.width()+1,0,true, 0.7, inputDelays);
		oplist.push_back(mult);

		vhdl<<"  "<<declare(resName+"_conv_B",  bType.width()+1)<<" <= "<< std::string("(\"0\"&")+bName+");";

		inPortMap(mult, "X", aName);
		inPortMap(mult, "Y", resName+"_conv_B");
	    }

	    outPortMap(mult, "R", resName);
	    vhdl<<instance(mult, resName+"_theMul");
	    syncCycleFromSignal(resName, mult->getSignalDelay("R"));
	}
	return res;
}

/* UU + UU -> UUU
	SU + UU -> SUUU
	SU + UUU -> SUUUU
	SUU + UU -> SUUU
	SU + SU -> SUU
*/
fixed_format_t FixedPointPolynomialEvaluator::AddType(const fixed_format_t &a, const fixed_format_t &b) const
{
	if(b.isSigned && !a.isSigned)
		return AddType(b,a);
	
	fixed_format_t res;
	if(a.isSigned==b.isSigned){
		res.isSigned=a.isSigned;
		res.msb=std::max(a.msb,b.msb)+1;
		res.lsb=std::min(a.lsb,b.lsb);
	}else{
		res.isSigned=true;
		res.msb=std::max(a.msb,b.msb)+2;
		res.lsb=std::min(a.lsb,b.lsb);
	}
	return res;
}

fixed_format_t FixedPointPolynomialEvaluator::AddStatement(std::string resName, std::string aName, const fixed_format_t &aType, std::string bName, const fixed_format_t &bType)
{
	fixed_format_t res=AddType(aType,bType);
	
	manageCriticalPath(getTarget()->adderDelay(res.width()));

	vhdl << "--- "<<res<<" <- "<<aType<<" + "<<bType<<"\n";
	vhdl<<declare(resName, res.width())<< " <= ";
	if(res.isSigned){
		vhdl<<"std_logic_vector(signed("<<ExtendExpr(res,aName,aType)<<") + signed("<<ExtendExpr(res,bName,bType)<<"))";
	}else{
		vhdl<<"std_logic_vector(unsigned("<<ExtendExpr(res,aName,aType)<<") + unsigned("<<ExtendExpr(res,bName,bType)<<"))";
	}
	vhdl<<";\n";

	return res;
}

static mpz_class MakeStandardPattern(unsigned type, const fixed_format_t &fmt)
{
	mpz_class value;
	type=type&0x3;
	if(type==0){
		value=0;
	}else if(type==1){
		value=(mpz_class(1)<<fmt.width())-1;
	}else if(type==2){
		value=mpz_class(1)<<(fmt.width()-1);
	}else{
		value=(mpz_class(1)<<(fmt.width()-1))-1;
	}
	return value;
}


void FixedPointPolynomialEvaluator::buildStandardTestCases(TestCaseList* tcl)
{		
	int degree=getPolynomialDegree();
	
	for(unsigned mask=0;mask<(1<<(2*degree+2));mask++){
		TestCase *tc=new TestCase(this);
		
		for(int i=0;i<=degree;i++){
			tc->setInputValue( join("a",i), MakeStandardPattern(mask>>(2*i), getCoefficientFormat(i)));
		}
		tc->setInputValue("Y", MakeStandardPattern(mask>>(2*degree), getInputFormat()));			
		
		emulate(tc);
		
		tcl->add(tc);
	}
}

void FixedPointPolynomialEvaluator::emulate(TestCase *tc)
{
	int degree=getPolynomialDegree();
	int prec=std::max(getInputFormat().width()*3,(int)getToolPrecision());	// Precision for mpfr
	
	std::stringstream debug;
	
	mpfr::mpreal x=DecodeRaw(getInputFormat(), tc->getInputValue("Y"), prec);
	debug<<"x="<<x;
	mpfr::mpreal acc=DecodeRaw(getCoefficientFormat(degree), tc->getInputValue(join("a",degree)), prec);
	debug<<", a"<<degree<<"="<<acc;
	for(int i=degree-1;i>=0;i--){
		mpfr::mpreal coeff=DecodeRaw(getCoefficientFormat(i), tc->getInputValue(join("a",i)), prec);
		debug<<", a"<<i<<"="<<coeff;
		acc = coeff + x*acc;
	}
	debug<<", y="<<acc;
	
	REPORT(DEBUG, debug.str());
	
	mpfr::mpreal rUp=ldexp(ceil(ldexp(acc, -getOutputFormat().lsb)), getOutputFormat().lsb);
	mpfr::mpreal rDown=ldexp(floor(ldexp(acc, -getOutputFormat().lsb)), getOutputFormat().lsb);
	
	tc->addComment(debug.str());
	tc->addExpectedOutput("R", EncodeRaw(getOutputFormat(), rUp));
	tc->addExpectedOutput("R", EncodeRaw(getOutputFormat(), rDown));
}

FixedPointPolynomialEvaluator *CreateFullPrecisionPolynomialEvaluator(
	const fixed_format_t &inputFormat,
	const std::vector<fixed_format_t> &coefficientFormats,
	int outputLsb,
	const mpfr::mpreal &errorBudget,
	Target *target
);

FixedPointPolynomialEvaluator *CreateFixedPointPolynomialEvaluator(
	const fixed_format_t &inputFormat,
	const std::vector<fixed_format_t> &coefficientFormats,
	int outputLsb,
	const mpfr::mpreal &errorBudget,
	Target *target
){
	return CreateFixedPointPolynomialEvaluator(
		"auto",
		inputFormat, coefficientFormats, outputLsb, errorBudget,target
	);
}

FixedPointPolynomialEvaluator *CreateFixedPointPolynomialEvaluator(
	std::string selector,
	const fixed_format_t &inputFormat,
	const std::vector<fixed_format_t> &coefficientFormats,
	int outputLsb,
	const mpfr::mpreal &errorBudget,
	Target *target
){
	if(selector=="auto"){
		selector="full-precision";
		if(::flopoco::verbose>=1)
			std::cerr<<"CreateFixedPointPolynomialEvaluator - Mapping auto -> "<<selector<<"\n";
	}
	
	if(selector=="full-precision"){
		return CreateFullPrecisionPolynomialEvaluator(
			inputFormat,
			coefficientFormats,
			outputLsb,
			errorBudget,
			target
		);
	}else{
		throw std::string("CreateFixedPointPolynomialEvaluator - Didn't understand selector '")+selector+"'";
	}
}

/*
FixedPointPolynomialEvaluator *CreateFixedPointPolynomialEvaluator(
	std::string selector,
	const fixed_format_t &inputFormat,
	const std::vector<polynomial_spec_t> &polynomials,
	int outputLsb,
	const mpfr::mpreal &errorBudget,
	Target *target
){
	std::vector<fixed_format_t> coeffFormats;
	for(int i=0;i<polynomials.size();i++){
		const std::vector<mpfr::mpreal> &coeffs=polynomials[i].coefficients;
		for(int j=0;j<coeffs.size();j++){
			if(j>=coeffFormats.size()){
				coeffFormats.push_back(FixedFormatFromValue(coeffs[j]));
			}else{
				coeffFormats[j]=Union(coeffFormats[j], coeffs[j]);
			}
		}
	}
	
	return CreateFixedPointPolynomialEvaluator(
		selector,
		inputFormat,
		coeffFormats,
		outputLsb,
		errorBudget,
		target
	);
}
*/

static void FixedPointPolynomialEvaluatorFactoryUsage(std::ostream &dst)
{
	OperatorFactory::classic_OP(dst, "FixedPointPolynomialEvaluator", "selector degree targetPrec errorBudget inputFormat coeffFormat+", false);
	dst << "    Create a fixed-point polynomial evaluator with given input and coefficient formats\n";
	dst << "           selector - Somehow specifies how you want it build (see below).\n";
	dst << "	      degree - Degree of polynomial\n";
	dst << "	      targetPrec - Binary weight of the faithful rounding, i.e. target error of 2^(targetPrec-1).\n";
	dst << "           errorBudget - How much extra error this evaluator can introduce, not including final rounding (0 < errorBudget < 2^(targetPrec-1)).\n";
	dst << "           inputFormat - Format for the polynomial argument (currently must be unsigned)\n";
	dst << "           coeffFormat - Format for each coefficient from a[0] up to a[degree] (currently must be signed).\n";
	dst << "      Formats are specified as \"U|S;msb;lsb\", where u or s says whether it is two's complement,\n";
	dst << "      msb gives the weight of the msb, and lsb gives the weight of the lsb. Some examples are:,\n";
	dst << "         \"U;15;0\"  : 16 bit unsigned integer, bits represent 2^15,2^14,...,2^1,2^0\n";
	dst << "         \"S;32;0\"  : 32 bit signed integer, bits represent -2^31,2^30,...,2^1,2^0\n";
	dst << "         \"U;-1;-16\"  : 16 bit unsigned number with 16 fractional bits in range [0,1)\n";
	dst << "         \"S;0;-15\"  : 16 bit signed number with 15 fractional bits in range [-1,1)\n";
	dst << "         \"S;7;-8\"  : 16 bit signed number with 8 fraction bits in range [-128,+128)\n";
	dst << "      selector can indicate a number of different things. Known values are:\n";
	dst << "          \"full-precision\" - Evaluator which does no rounding at all until it rounds to targetPrec.\n";
	dst << "\n";
}

static Operator *FixedPointPolynomialEvaluatorFactoryParser(Target *target ,const std::vector<std::string> &args,int &consumed)
{
	if(args.size()<6)
		throw std::string("FixedPointPolynomialEvaluatorFactory - Not enough arguments, check usage.");
	
	std::string selector=args[0];
	
	int degree=boost::lexical_cast<int>(args[1]);
	if((int)args.size()<6+degree)
		throw std::string("PolynomialEvaluatorFactory - Not enough arguments, check usage.");
	
	int targetPrec=boost::lexical_cast<int>(args[2]);
	mpfr::mpreal errorBudget;
	parseSollyaConstant(get_mpfr_ptr(errorBudget), args[3]);
	fixed_format_t inputFormat=ParseFixedFormat(args[4]);
	std::vector<fixed_format_t> coeffFormats;
	for(int i=0;i<=degree;i++){
		coeffFormats.push_back(ParseFixedFormat(args.at(5+i)));
	}
	consumed=6+degree;
	
	return CreateFixedPointPolynomialEvaluator(
		selector,
		inputFormat,
		coeffFormats,
		targetPrec,
		errorBudget,
		target
	);
}

void FixedPointPolynomialEvaluator_registerFactory()
{
	DefaultOperatorFactory::Register(
		"FixedPointPolynomialEvaluator",
		"operator",
		FixedPointPolynomialEvaluatorFactoryUsage,
		FixedPointPolynomialEvaluatorFactoryParser
	);
}





}; // random
}; // flopoco
