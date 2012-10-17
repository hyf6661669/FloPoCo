#include "comparable_float_type.hpp"

namespace flopoco
{
namespace random
{
	
ComparableFloatEncoder *ComparableFloatType::MakeEncoder(Target *target, map<string, double> inputDelays)
{
	return new ComparableFloatEncoder(target, *this, inputDelays);
}	

ComparableFloatEncoder::ComparableFloatEncoder(
	Target *target,
	ComparableFloatType _type,
	map<string, double> inputDelays
)
	: Operator(target, inputDelays)
	, type(_type)
	, wE(type.wE)
	, wF(type.wF)
{
	ostringstream name;
	name << "ComparableFloatEncoded_wE"<<type.wE<<"_wF"<<type.wF<<"_uid"<<getNewUId();
	setName(name.str());
	setCopyrightString("Imperial College 2012");
	
	addFPInput ("iX" , wE, wF);
	addOutput("oY", type.Width());
	assert(type.Width()==2+wE+wF);	// Needs to be modified if non-regular numbers come in
	
	setCriticalPath( getMaxInputDelays(inputDelays) );
	vhdl << declare("prefix_in",3) << " <= iX"<<range(2+wE+wF,wE+wF)<<";\n";
	vhdl << declare("exponent_in",wE) << " <= iX"<<range(wE+wF-1,wF)<<";\n";
	vhdl << declare("fraction_in",wF) << " <= iX"<<range(wF-1,0)<<";\n";
	
	manageCriticalPath(target_->adderDelay(wF));
	
	vhdl << declare("prefix_out",2) << " <= "
		"\"11\" when prefix_in=\"010\" else\n" <<	// positive regular
		"\"10\" when prefix_in=\"000\" else\n" <<	// positive zero
		"\"01\" when prefix_in=\"001\" else\n" <<	// negative zero
		"\"00\" when prefix_in=\"011\" else\n" <<	// negative zero
		"\"XX\";\n";
	

	vhdl << declare("exponent_out",wE) << " <= "
		"(exponent_in + \"1\") when prefix_in=\"010\" else\n" <<	// positive regular
		"("<<zg(wE)<<") when (prefix_in=\"000\" or prefix_in=\"001\") else\n" <<	// both zeros
		"("<<og(wE)<<" - exponent_in - \"1\") when prefix_in=\"011\" else\n" <<	// negative zero
		xg(wE)<<";\n";
	
	vhdl << declare("fraction_out",wF) << " <= "
		"fraction_in when prefix_in=\"010\" else\n" <<	// positive regular
		"("<<zg(wF)<<") when (prefix_in=\"000\" or prefix_in=\"001\") else\n" <<	// both zeros
		"("<<og(wF)<<" - fraction_in) when prefix_in=\"011\" else\n" <<	// negative zero
		xg(wF)<<";\n";
	
	vhdl << "oY <= prefix_out & exponent_out & fraction_out;\n";
	outDelayMap["oY"] = getCriticalPath();
}

void ComparableFloatEncoder::emulate(TestCase * tc)
{
	mpz_class x=tc->getInputValue("iX");
	FPNumber fpx(wE, wF);
	fpx=x;
	
	mpfr_t tmp;
	mpfr_init2(tmp, wF+1);
	fpx.getMPFR(tmp);
	
	mpz_class y=type.ToBits(tmp);
	tc->addExpectedOutput("oY", y);
	
	mpfr_fprintf(stderr, "  emulate : iX=%Rg, oY=%Zu\n", tmp, y.get_mpz_t());
	
	mpfr_clear(tmp);
}

void ComparableFloatEncoder::buildStandardTestCases(TestCaseList* tcl)
{
	mpfr_t tmp;
	mpfr_init2(tmp, wF+1);
	
	mpfr_set_zero(tmp, +1);
	addTestCase(tcl, tmp);
	
	mpfr_set_zero(tmp, -1);
	addTestCase(tcl, tmp);
	
	mpfr_set_d(tmp, 1.0, MPFR_RNDN);
	addTestCase(tcl, tmp);
	
	mpfr_nextbelow(tmp);
	addTestCase(tcl, tmp);
	
	mpfr_nextbelow(tmp);
	addTestCase(tcl, tmp);
	
	mpfr_set_d(tmp, 1.0, MPFR_RNDN);
	mpfr_nextabove(tmp);
	addTestCase(tcl, tmp);
	
	mpfr_set_d(tmp, 1.0, MPFR_RNDN);
	mpfr_set_exp(tmp, -(1<<(wE-1)));
	addTestCase(tcl, tmp);
	
	mpfr_nextabove(tmp);
	addTestCase(tcl, tmp);
	
	mpfr_set_d(tmp, 1.0, MPFR_RNDN);
	mpfr_set_exp(tmp, (1<<(wE-1))-1);
	addTestCase(tcl, tmp);
	
	mpfr_nextabove(tmp);
	addTestCase(tcl, tmp);
	
	mpfr_set_d(tmp, 1.0, MPFR_RNDN);
	mpfr_mul_2si(tmp, tmp, 1<<(wE-1), MPFR_RNDN);
	mpfr_nextbelow(tmp);
	addTestCase(tcl, tmp);
}

void ComparableFloatEncoder::addTestCase(TestCaseList* tcl, mpfr_t val)
{
	mpfr_fprintf(stderr, " iX=%Rg\n", val);
	TestCase *tc=new TestCase(this);
	FPNumber fpx(wE, wF, val);
	mpfr_fprintf(stderr, "  excep=%Zu\n", fpx.getExceptionSignalValue().get_mpz_t());
	tc->addFPInput("iX", &fpx);
	emulate(tc);
	tcl->add(tc);
	
	mpfr_neg(val, val, MPFR_RNDN);
	
	tc=new TestCase(this);
	fpx=val;
	mpfr_fprintf(stderr, "  excep=%Zu\n", fpx.getExceptionSignalValue().get_mpz_t());
	tc->addFPInput("iX", &fpx);
	emulate(tc);
	tcl->add(tc);
	
	mpfr_neg(val, val, MPFR_RNDN);
}

TestCase* ComparableFloatEncoder::buildRandomTestCase(int i)
{
	fprintf(stderr, "Wibble\n");
	
	mpfr_t tmp, tmp2;
	mpfr_init2(tmp, wF+1);
	mpfr_init2(tmp2, wF+1);
	
	int e=(int)floor(drand48()*(1<<wE))-(1<<(wE-1));
	mpz_class fraction=getLargeRandom(wF)+(mpz_class(1)<<wF);
	
	mpfr_set_z(tmp, fraction.get_mpz_t(), MPFR_RNDN);
	mpfr_mul_2si(tmp, tmp, e-wF, MPFR_RNDN);
	
	if(drand48()>0.5)
		mpfr_setsign(tmp, tmp, -1, MPFR_RNDN);
	
	mpfr_fprintf(stderr, " random: iX=%Rg\n", tmp);
	
	mpz_class vv=type.ToBits(tmp);
	type.FromBits(tmp2, vv);
	assert(mpfr_equal_p(tmp2, tmp));
	
	TestCase *tc=new TestCase(this);
	FPNumber fpx(wE, wF, tmp);
	tc->addFPInput("iX", &fpx);
	emulate(tc);
	
	mpfr_clear(tmp);
	mpfr_clear(tmp2);
	
	return tc;
}

}; // random
}; // flopoco

