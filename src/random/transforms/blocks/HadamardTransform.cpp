// general c++ library for manipulating streams
#include <iostream>
#include <sstream>
#include <math.h>	// for NaN



/* header of libraries to manipulate multiprecision numbers
  There will be used in the emulate function to manipulate arbitraly large
  entries */
#include "gmp.h"
#include "mpfr.h"
#include "FPNumber.hpp"

// include the header of the Operator
#include "HadamardTransform.hpp"
#include "CLTTransform.hpp"

#include "random/distributions/sum_distribution.hpp"
#include "random/utils/fft/convolve_mpreal.hpp"

using namespace std;

namespace flopoco
{

extern vector<Operator *> oplist;	
	
namespace random
{
	

void HadamardTransform::Connect(std::string dstName, int dstIdx, std::string srcName, int srcL, int srcR, int dir, int srcW)
{	
	setCycleFromSignal(join(srcName+"_",srcL));
	nextCycle();
	
	
	
	vhdl<<declare(join(dstName+"_",dstIdx,"_left"), srcW+1) << "<="<< signExtend(join(srcName+"_",srcL),srcW+1)<<";\n";
	vhdl<<declare(join(dstName+"_",dstIdx,"_right"), srcW+1) << "<="<<signExtend(join(srcName+"_",srcR),srcW+1)<<";\n";
	
	vhdl << declare( join(dstName+"_",dstIdx), srcW+1 )
		<<" <= " << join(dstName+"_",dstIdx,"_left")
		<< (dir>0 ?" + ":" - ")
		<<join(dstName+"_",dstIdx,"_right") <<";\n";
	
	if(target_->getVendor() == "Xilinx"){
		// If we don't do this, then Xilinx tools notice that the LSB of table-transforms can be collapsed
		// (as msb(-x)==msb(x)), and starts to do strange things.
		addAttribute("equivalent_register_removal",  "string",  join(dstName+"_",dstIdx,"_d1:signal"), "no");
		addAttribute("keep",  "string",  join(dstName+"_",dstIdx,"_d1:signal"), "true");
	}
}

void HadamardTransform::HadamardStage(std::string dstName, std::string srcName, int totalSize, int size, int srcW)
{
	for(int offset=0;offset<totalSize;offset+=size){
		for(int i=0;i<size/2;i++){
			Connect(dstName, offset+i, srcName, offset+i, offset+i+size/2, +1, srcW);
		}
		for(int i=size/2;i<size;i++){
			Connect(dstName, offset+i, srcName, offset+i-size/2, offset+i, -1, srcW);
		}
	}
}

void HadamardTransform::Hadamard(std::string name, int log2size, int srcW)
{
	for(int i=1;i<=log2size;i++){
		HadamardStage(join(name+"_",i), join(name+"_",i-1), 1<<log2size, 1<<i, srcW+i-1);
		nextCycle();
	}
}
	
HadamardTransform::HadamardTransform(Target* target,
	int log2n,
	RngTransformOperator *base,
	bool extraPipeline
)
	: RngTransformOperator(target)
	, m_log2n(log2n)
	, m_n(1<<log2n)
	, m_baseWidth(base->nonUniformOutputWidth(0))
	, m_base(base)
	, m_extraPipeline(extraPipeline)
{
	std::stringstream acc;
	acc<<"HadamardTransform_n"<<m_n<<"_uid"<<getNewUId();
	setName(acc.str());
	
	oplist.push_back(base);
	
	m_uniformInputBits=base->uniformInputBits() * (1<<log2n);
	m_uniformInputName="iUniformBits";
	
	assert(base->nonUniformOutputCount()==1);	
	
	m_nonUniformOutputCount=1<<log2n;
	m_nonUniformOutputWidth=m_baseWidth+log2n;
	m_nonUniformOutputNameBase="oRng";
	
	addInput(m_uniformInputName, m_uniformInputBits);
	
	syncCycleFromSignal(m_uniformInputName);
	
	for(int i=0;i<(1<<log2n);i++){
		std::string unifBitsRange=range((i+1)*base->uniformInputBits()-1,i*base->uniformInputBits());
		vhdl<<declare(join("unif_",i), base->uniformInputBits())<<"<=iUniformBits"<<unifBitsRange<<";\n";
		inPortMap(base, base->uniformInputName(), join("unif_",i));
		outPortMap(base, base->nonUniformOutputName(0), join("stage_0_",i));
		vhdl<<instance(base, join("baseGen_",i));
	}
	
	if(!m_extraPipeline){
		for(int i=0;i<(1<<log2n);i++){
			syncCycleFromSignal(join("stage_0_",i));
			
			if(target_->getVendor() == "Xilinx"){
				// If we don't do this, then Xilinx tools notice that the LSB of table-transforms can be collapsed
				// (as msb(-x)==msb(x)), and starts to do optimisations we don't want.
				addAttribute("equivalent_register_removal",  "string",  join("stage_0_",i,":signal"), "no");
				addAttribute("keep",  "string",  join("stage_0_",i,":signal"), "true");
			}
		}
		
		Hadamard("stage", log2n, base->nonUniformOutputWidth(0));
		
		for(int i=0;i<(1<<log2n);i++){
			addOutput(join("oRng",i), m_nonUniformOutputWidth);
			vhdl<<join(m_nonUniformOutputNameBase,i)<<" <= "<<join("stage_",log2n,"_",i)<<";\n";
		}
	}else{
		// This is a Xilinx specific optimisation which uses RPMs and the extra latches in Virtex-5 and up to
		// pipeline each 2x2 hadamard so that it has both input and output registers. It makes it very
		// routable, and easy to get high clock rates without long run times.
		// Unfortunately it's easiest to do as templates VHDL, so we need to link to an ugly blob of VHDL.
		
		int slice_height;
		
		if(target->getID()=="Virtex6"){
			slice_height=4;
		}else{
			throw std::string("HadamardTransform - Haven't specialised pipelining for target '")+target->getID()+"' yet.";
		}		
		
		syncCycleFromSignal("stage_0_0");
		declare("src", m_baseWidth*(1<<log2n));
		
		for(int i=0;i<(1<<log2n);i++){
			vhdl <<"src"<<range((i+1)*m_baseWidth-1,i*m_baseWidth)<<" <= "<<join("stage_0_",i)<<";\n";
		}
		
		// Florent, forgive me!
		vhdl << "impl : entity work.rpm_nxn_hadamard_grow\n";
		vhdl << "    generic map(LOG2N=>"<<log2n<<", W=>"<<m_baseWidth<<", SLICE_HEIGHT=>"<<slice_height<<")\n";
		vhdl << "    port map(clk=>clk, src=>src";
		
		// HACK!
		for(int i=0;i<log2n;i++){
			// Each transform stage has a latency of two
			nextCycle();
			nextCycle();
		}
		
		vhdl<<", dst=>dst);\n\n";
		
		int ow=m_baseWidth+log2n;
		declare("dst", ow*(1<<log2n));		
		
		for(int i=0;i<(1<<log2n);i++){
			addOutput(join("oRng",i), m_nonUniformOutputWidth);
			vhdl<<join("oRng",i)<<" <= dst"<<range(ow*(i+1)-1,ow*i)<<";\n";
		}
	}
}

HadamardTransform::~HadamardTransform()
{}

void HadamardTransform::emulate(TestCase * tc)
{
	mpz_class bits=tc->getInputValue(m_uniformInputName);
	
	std::vector<mpz_class> curr(m_n), next(m_n);
	
	for(unsigned i=0;i<m_n;i++){
		mpz_class inBits;
		mpz_tdiv_r_2exp(inBits.get_mpz_t(), bits.get_mpz_t(), m_base->uniformInputBits());
		
		TestCase *tmp=new TestCase(m_base);
		tmp->addInput(m_base->uniformInputName(), inBits);
		m_base->emulate(tmp);
		std::vector<mpz_class> xx=tmp->getExpectedOutputValues(m_base->nonUniformOutputName(0));
		if(xx.size()!=1)
			throw std::string("HadamardTransform::emulate - base generator returned multiple expected values.");
		curr[i]=fromTwosComplement(xx[0], m_base->nonUniformOutputWidth(0));
		//mpfr_fprintf(stderr, "%d : input=%Zx, twosComp=%Zx, value=%Zd\n", i, inBits.get_mpz_t(), xx[0].get_mpz_t(), curr[i].get_mpz_t());
		
		bits=bits>>m_base->uniformInputBits();
	}
	
	for(unsigned j=1;j<=m_log2n;j++){
		unsigned size=1<<j;
		for(unsigned offset=0;offset<m_n;offset+=size){
			for(unsigned i=0;i<size/2;i++){
				// Connect(dstName, offset+i, srcName, offset+i, offset+i+size/2, +1, srcW);
				next[offset+i] = curr[offset+i] + curr[offset+i+size/2];
			}
			for(unsigned i=size/2;i<size;i++){
				//Connect(dstName, offset+i, srcName, offset+i-size/2, offset+i, -1, srcW);
				next[offset+i] = curr[offset+i-size/2] - curr[offset+i];
			}
		}
		std::swap(curr, next);
	}
	
	for(unsigned i=0;i<m_n;i++){
		mpz_class v=toTwosComplement(curr[i], m_baseWidth+m_log2n);
		//mpfr_fprintf(stderr, "  out=%Zd, twos=%Zx, fwd=%Zd\n", curr[i].get_mpz_t(), v.get_mpz_t(), fromTwosComplement(v, m_baseWidth+m_log2n).get_mpz_t());
		tc->addExpectedOutput(join(m_nonUniformOutputNameBase,i), v);
	}
}

typename Distribution<mpfr::mpreal>::TypePtr HadamardTransform::nonUniformOutputDistribution(int i, unsigned prec) const
{
	if((i<0) || (i>=(int)m_n))
		throw std::string("Distribution index out of range.");
	
	if(!m_distribution || (prec>=m_distributionPrec)){
		IRngTransformDistributions *src=dynamic_cast<IRngTransformDistributions*>(m_base);
		if(src==NULL)
			throw std::string("HadamardTransform::nonUniformOutputDistribution - Base does not know its distribution.");
		
		typename DiscreteDistribution<mpfr::mpreal>::TypePtr part(boost::dynamic_pointer_cast<DiscreteDistribution<mpfr::mpreal> >(src->nonUniformOutputDistribution(0, prec)));
		assert(part->Cdf(0).get_prec() >= (int)prec);
		assert(part->Pmf(0).get_prec() >= (int)prec);
		m_distribution=SelfAddDistributions<mpfr::mpreal>(part, (int)m_n);
		m_distributionPrec=prec;
	}
	assert(m_distribution->Pmf(0).get_prec() >= (int)prec);
	assert(m_distribution->Cdf(0).get_prec() >= (int)prec);
	return m_distribution;
}

	
TestCase* HadamardTransform::buildRandomTestCase(int i)
{
	TestCase *tc=new TestCase(this);
	
	tc->addInput(m_uniformInputName, getLargeRandom(m_uniformInputBits));
	emulate(tc);
	
  	return tc;
}

};
};
