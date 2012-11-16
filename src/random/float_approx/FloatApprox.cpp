#include <iostream>
#include <sstream>

#include "range_polys.hpp"

#include "FloatApprox.hpp"

#include "mpreal.h"
#include "UtilSollya.hh"

#include "random/utils/operator_factory.hpp"
#include "random/utils/make_table.hpp"

using namespace std;

namespace flopoco{
namespace random{
  
using namespace float_approx;

class FloatApproxOperator : public Operator {
private:
  // Basic operator parameters describing the problem
  int m_wDomainE, m_wDomainF, m_wRangeE, m_wRangeF;
  mpfr::mpreal m_domainMin, m_domainMax;
  Function m_f;
  unsigned m_degree;
  mpfr::mpreal m_maxError;

  // The information we've built up about the function
  Range m_range;
  RangePolys m_polys;

  std::vector<mpz_class> m_tableContents;

  Operator *m_codec;
  Operator *m_quantiser;
  Operator *m_table;
  PolynomialEvaluator *m_poly;
public:
  FloatApproxOperator(Target* target,
    int wDomainE, int wDomainF, mpfr::mpreal domainMin, mpfr::mpreal domainMax,
    int wRangeE, int wRangeF,
    const Function &f,
    unsigned degree,
    mpfr::mpreal maxError
  )
    : Operator(target)
    // Capture properties
    , m_wDomainE(wDomainE)
    , m_wDomainF(wDomainF)
    , m_wRangeE(wRangeE)
    , m_wRangeF(wRangeF)
    , m_domainMin(domainMin)
    , m_domainMax(domainMax)
    , m_f(f)
    , m_degree(degree)
    , m_maxError(maxError)
    // Start building stuff (though not much happens yet)
    , m_range(f, wDomainF, wRangeF, domainMin.mpfr_ptr(), domainMax.mpfr_ptr())
  {       
    std::stringstream acc;
    acc<<"FloatApprox_uid"<<getNewUId();
    std::string name=acc.str();
    setName(name);
    
    REPORT(INFO, "Making range monotonic.");
    m_range.make_monotonic_or_range_flat();
    REPORT(INFO, "  -> no. of segments="<<m_range.m_segments.size());
    int nMonoSegments=m_range.m_segments.size();
    
    REPORT(INFO, "Flattening domain.");
    m_range.flatten_domain();
    REPORT(INFO, "  -> no. of segments="<<m_range.m_segments.size());
    int nDomFlatSegments=m_range.m_segments.size();
    
    REPORT(INFO, "Flattening range.");
    m_range.flatten_range();
    REPORT(INFO, "  -> no. of segments="<<m_range.m_segments.size());
    int nRanFlatSegments=m_range.m_segments.size();
    
    m_polys=RangePolys(&m_range, m_degree);
    
    REPORT(INFO, "Splitting into polynomials with error < "<<maxError);
    m_polys.split_to_error(maxError.toDouble());
    REPORT(INFO, "  -> no. of segments="<<m_range.m_segments.size());
    int nErrSplitSegments=m_range.m_segments.size();
    
    int guard=0;
    for(guard=0;guard<=9;guard++){
      if(guard==8){
        throw std::string("FloatApprox - More than 8 guard bits needed, this probably means something is going horribly wrong.");
      }
      REPORT(INFO, "Trying to find faithful polynomials with "<<guard<<" guard bits.");
      try{
        m_polys.calc_faithful_fixed_point(guard);
      }catch(std::string msg){
        if(msg!="calc_faithful_fixed_point - fixed-point poly was not faithful."){
          throw;    // All other problems should throw
        }
        continue;
      }
      break;
    }
    REPORT(INFO, "Successful, using "<<guard<<" guard bits.");
    
    REPORT(INFO, "Building fixed-point coefficient tables.");
    m_polys.build_concrete(guard);
    
    int nFinalSegments=m_range.m_segments.size();
    
    int wSegmentIndex=(int)ceil(log(m_range.m_segments.size())/log(2.0));
    
    REPORT(INFO, "Constructing static quantiser.");
    m_codec=ComparableFloatType(wDomainE,wDomainF).MakeEncoder(target);
    oplist.push_back(m_codec);
    m_quantiser=m_polys.make_static_quantiser(target, wDomainE);
    oplist.push_back(m_quantiser);
    
    REPORT(INFO, "Constructing lookup table.");
    int coeffWidths=0;
    for(unsigned i=0;i<=m_degree;i++){
      coeffWidths+=m_polys.m_concreteCoeffMsbs[i]-m_polys.m_concreteCoeffLsbs[i]+1;
    }
    int tableWidth=3+wRangeE+coeffWidths;
    m_tableContents=m_polys.build_ram_contents(guard, wRangeE);
    m_table=MakeSinglePortTable(target, name+"_table", tableWidth, m_tableContents);
    oplist.push_back(m_table);
    
    REPORT(INFO, "Constructing polynomial evaluator.");
    m_poly=m_polys.make_polynomial_evaluator(target);
    oplist.push_back(m_poly);
    REPORT(INFO, "  width of poly eval result is "<<m_poly->getRWidth()<<", weight of msb is "<<m_poly->getRWeight());
    
    REPORT(INFO, "Now constructing VHDL for FloatApprox.");
    
    addFPInput("iX", wDomainE, wDomainF);
    addFPOutput("oY", wRangeE, wRangeF);
    
    inPortMap(m_codec, "iX", "iX");
    outPortMap(m_codec, "oY", "comparable_iX");
    vhdl<<instance(m_codec, "comparable_codec");
    syncCycleFromSignal("comparable_iX");
    
    inPortMap(m_quantiser, "iX", "comparable_iX");
    outPortMap(m_quantiser, "oY", "table_index");
    vhdl<<instance(m_quantiser, "quantiser");
    syncCycleFromSignal("table_index");
    
    if(nFinalSegments>=256){
      nextCycle();  // BRAM input stage
    }
    inPortMap(m_table, "X", "table_index");
    outPortMap(m_table, "Y", "table_contents");
    vhdl<<instance(m_table, "table");
    syncCycleFromSignal("table_contents");
    if(nFinalSegments>=256){
      nextCycle();  // BRAM output register
    }
    
    ///////////////////////////////////
    // Split the coefficients into the different parts
    int offset=0;
    for(unsigned i=0;i<=m_degree;i++){
      int w=m_polys.m_concreteCoeffMsbs[i]-m_polys.m_concreteCoeffLsbs[i]+1;
      vhdl<<declare(join("coeff_",i),w)<<" <= table_contents"<<range(w+offset-1,offset)<<";\n";
      offset+=w;
    }
    vhdl<<declare("coeff_prefix",3+wRangeE)<<" <= table_contents"<<range(tableWidth-1,offset)<<";\n";
    vhdl<<declare("fraction_iX", wDomainF)<<" <= iX"<<range(wDomainF-1,0)<<";\n";
    
    ////////////////////////////////
    // Now connect to polynomial evaluator
    for(unsigned i=0;i<=m_degree;i++){
      inPortMap(m_poly, join("a",i), join("coeff_",i));
    }
    inPortMap(m_poly, "Y", "fraction_iX");
        
    outPortMap(m_poly, "R", "result_fraction");
    vhdl<<instance(m_poly, "poly");
    syncCycleFromSignal("result_fraction");
    
    // We have to convert from this format to the target format
    PolynomialEvaluator::format_t result_format=m_poly->getOutputFormat();
    REPORT(INFO, "poly-eval result_format={isSigned="<<result_format.isSigned<<",msb="<<result_format.msb<<",lsb="<<result_format.lsb);
    
    if(!result_format.isSigned)
      throw std::string("Currently FloatApprox assumes the output of PolynomialEvaluator will always be signed.");
    if(result_format.lsb > -wRangeF)
      throw std::string("The LSB of the poly-eval is higher than that requested.");
    if(result_format.msb < -1)
      throw std::string("Currently FloatApprox needs PolynomialEvaluator to return at least the interval [-0.5,0.5). This is my fault not yours, try increasing approximation range to cover entire binade.");
    
    int drop_bits=(- wRangeF) - result_format.lsb ;
    
    vhdl<<"oY <= coeff_prefix & result_fraction"<<range(wRangeF-1+drop_bits,drop_bits)<<";\n";
    
    addOutput("debug_result_fraction", wRangeF);
    vhdl<<"debug_result_fraction<=result_fraction"<<range(wRangeF-1+drop_bits,drop_bits)<<";\n";
    
    addOutput("debug_result_prefix", 3+wRangeE);
    vhdl<<"debug_result_prefix<=coeff_prefix;\n";
    
    addOutput("debug_segment", wSegmentIndex);
    vhdl<<"debug_segment<=table_index;\n";
    
    addOutput("debug_table_contents", tableWidth);
    vhdl<<"debug_table_contents<=table_contents;\n";
    
    vhdl<<"--nSeg_Monotonic = "<<nMonoSegments<<"\n";
    vhdl<<"--nSeg_FlatDomain = "<<nDomFlatSegments<<"\n";
    vhdl<<"--nSeg_FlatDomainAndRange = "<<nRanFlatSegments<<"\n";
    vhdl<<"--nSeg_ErrSplit = "<<nErrSplitSegments<<"\n";
    vhdl<<"--nSeg_Final = "<<nFinalSegments<<"\n";
  }

  void emulate(TestCase * tc)
  {
    mpz_class iX=tc->getInputValue("iX");
    FPNumber vx(m_wDomainE,m_wDomainF, iX);
    
    int prec=std::max(m_wRangeF+64, (int)getToolPrecision());
    
    mpfr_t x, exact, rounded;
    mpfr_init2(x, m_wDomainF+1);
    mpfr_init2(exact, prec);
    
    vx.getMPFR(x);
    m_f.eval(exact, x);
    
    mpfr_init2(rounded, m_wRangeF+1);
    mpfr_set(rounded, exact, MPFR_RNDD);
    FPNumber rd(m_wRangeE,m_wRangeF, rounded);
    
    mpfr_set(rounded, exact, MPFR_RNDU);
    FPNumber ru(m_wRangeE, m_wRangeF, rounded);
    
    mpz_class v=rd.getSignalValue();
    tc->addExpectedOutput("oY", v);
    tc->addExpectedOutput("debug_result_prefix", v>>m_wRangeF);
    mpz_cdiv_r_2exp(v.get_mpz_t(), v.get_mpz_t(), m_wRangeF);
    tc->addExpectedOutput("debug_result_fraction", v);
    
    v=ru.getSignalValue();
    tc->addExpectedOutput("oY", v);
    tc->addExpectedOutput("debug_result_prefix", v>>m_wRangeF);
    mpz_cdiv_r_2exp(v.get_mpz_t(), v.get_mpz_t(), m_wRangeF);
    tc->addExpectedOutput("debug_result_fraction", v);
    
    mpfr_clears(x, exact, rounded, (mpfr_ptr)0);
  }
  
  std::vector<mpz_class> emulate(const char *retSig, Operator *op, const char *n1, const std::vector<mpz_class> &v1) const
  {
    std::set<mpz_class> res;
    for(unsigned i=0;i<v1.size();i++){
      TestCase tc(op);
      tc->addInput(n1, v1[i]);
      op->emulate(tc);
      const std::vector<mpz_class> &values=tc->getExpectedOutputValues(retSig);
      res.insert(values.begin(), values.end());
    }
    return std::vector<mpz_class>(res.begin(), res.end());
  }
  
  std::vector<mpz_class> emulate(
    const char *retSig, Operator *op,
    const char *n1, const std::vector<mpz_class> &v1,
    const char *n2, const std::vector<mpz_class> &v2
  ) const
  {
    std::set<mpz_class> res;
    for(unsigned i1=0;i1<v1.size();i1++){
      for(unsigned i2=0;i2<v2.size();i2++){
        TestCase tc(op);
        tc->addInput(n1, v1[i1]);
        tc->addInput(n2, v2[i2]);
        op->emulate(tc);
        const std::vector<mpz_class> &values=tc->getExpectedOutputValues(retSig);
        res.insert(values.begin(), values.end());
      }
    }
    return std::vector<mpz_class>(res.begin(), res.end());
  }
  
  void emulate_ops(TestCase *tc)
  {
    std::vector<mpz_class> iX=std::vector<mpz_class>(1, tc->getInputValue("iX"));
    
    std::vector<mpz_class> comparable_iX=emulate("oY", m_codec, "iX", iX);
    std::vector<mpz_class> table_index=emulate("oY", m_quantiser, "iX", comparable_iX);
    std::vector<mpz_class> table_contents=emulate("Y", m_table, "X", table_index);
    
    for(unsigned i=0;i<table_contents.size();i++){
      std::vector<std::pair<std::string,std::vector<mpz_class> > > polyEvalArgs;
    
      int offset=0;
      for(unsigned i=0;i<=m_degree;i++){
        int w=m_polys.m_concreteCoeffMsbs[i]-m_polys.m_concreteCoeffLsbs[i]+1;
        vhdl<<declare(join("coeff_",i),w)<<" <= table_contents"<<range(w+offset-1,offset)<<";\n";
        offset+=w;
      }
      vhdl<<declare("coeff_prefix",3+wRangeE)<<" <= table_contents"<<range(tableWidth-1,offset)<<";\n";
      vhdl<<declare("fraction_iX", wDomainF)<<" <= iX"<<range(wDomainF-1,0)<<";\n";
  }

  void buildStandardTestCases(TestCaseList* tcl)
  {
    FPNumber dom(m_wDomainE, m_wDomainF);
    
    mpz_class index;
    Range::segment_it_t curr=m_range.m_segments.begin();
    while(curr!=m_range.m_segments.end()){
      dom=curr->domainStart;
      TestCase *tc=new TestCase(this);
      tc->addFPInput("iX", &dom);
      emulate(tc);
      tc->addExpectedOutput("debug_segment", index);
      tc->addExpectedOutput("debug_table_contents", m_tableContents[index.get_ui()]);
      tcl->add(tc);
      
      dom=curr->domainFinish;
      tc=new TestCase(this);
      tc->addFPInput("iX", &dom);
      emulate(tc);
      tc->addExpectedOutput("debug_segment", index);
      tc->addExpectedOutput("debug_table_contents", m_tableContents[index.get_ui()]);
      tcl->add(tc);
      
      ++curr;
      ++index;
    }
  }

  //TestCase* buildRandomTestCase(int i)
  //{
  //  throw std::string("Not Implemented.");
  //}
}; // FloatApproxOperator



static void FloatApproxFactoryUsage(std::ostream &dst)
{
	OperatorFactory::classic_OP(dst, "FloatApprox", "wDomE wDomF domMin domMax wRanE wRanF f degree", false);
	dst << "    Generates a float->float function approximation for y=f(x) where domMin<=x<=domMax.\n";
	dst << "	      (wDomE,wDomF) - Floating-point format of input.\n";
	dst << "	      [domMin,domMax] - Inclusive domain of approximation\n";
	dst << "	      (wRanE,wRanF) - Floating point format for output.\n";
    dst << "	      f - Any twice differentiable function.\n";
	dst << "      The transform will take a (wDomE,wDomF) input float, and produce a (wRanE,wRanF) output float.\n";
    dst << "         NOTE : currently, the domain of f must be strictly positive, i.e. 0<domMin, and\n";
    dst << "                   the image of f must also be strictly positive, i.e. f(x)>0 for domMin<=x<=domMax.\n";
}

static Operator *FloatApproxFactoryParser(Target *target ,const std::vector<std::string> &args,int &consumed)
{
	unsigned nargs = 8;
	if (args.size()<nargs)
		throw std::string("FloatApproxFactoryParser - Not enough arguments, check usage.");
	consumed += nargs;
	
	int wDomE = atoi(args[0].c_str());
	int wDomF = atoi(args[1].c_str());
    if((wDomE<1) || (wDomF<1))
      throw std::string("FloatApproxFactoryParser - wDomE and wDomF must be positive.");
    
    mpfr::mpreal domMin(0, wDomF+1);
    mpfr::mpreal domMax(0, wDomF+1);
    parseSollyaConstant(domMin.mpfr_ptr(), args[2], MPFR_RNDU);
    parseSollyaConstant(domMax.mpfr_ptr(), args[3], MPFR_RNDD);
    if(domMin<=0)
      throw std::string("FloatApproxFactoryParser - domMin must be positive.");
    if(domMax <= domMin)
      throw std::string("FloatApproxFactoryParser - Must have domMin < domMax.");
    
    int wRanE=atoi(args[4].c_str());
    int wRanF=atoi(args[5].c_str());
    if((wRanE<1) || (wRanF<1))
      throw std::string("FloatApproxFactoryParser - wRanE and wRanF must be positive.");
    
	Function f(args[6]);
    
    int degree=atoi(args[7].c_str());
    if(degree<0)
      throw std::string("FloatApproxFactoryParser - degree must be at least 0.");
    
    mpfr::mpreal maxError(pow(2.0,-wRanF-1), getToolPrecision());

	return new FloatApproxOperator(target,
      wDomE, wDomF, domMin.mpfr_ptr(), domMax.mpfr_ptr(),
      wRanE, wRanF,
      f,
      degree,
      maxError.mpfr_ptr()
    );
}

void FloatApprox_registerFactory()
{
	DefaultOperatorFactory::Register(
		"FloatApprox",
		"operator",
		FloatApproxFactoryUsage,
		FloatApproxFactoryParser
	);
}


};
};
