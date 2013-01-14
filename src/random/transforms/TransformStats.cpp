#include <iostream>
#include <sstream>
#include <math.h>

#include "RngTransformOperator.hpp"

#include "random/distributions/discrete_gaussian_distribution.hpp"
#include "random/distributions/quantised_gaussian_distribution.hpp"
#include "random/distributions/chi2_partition.hpp"

#include "random/utils/operator_factory.hpp"

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/make_shared.hpp>

#include "UtilSollya.hh"

using namespace std;

namespace flopoco
{

extern vector<Operator *> oplist;	
	
namespace random
{
	
	struct segment{
		double sigma;
		mpfr::mpreal start_x;
		mpfr::mpreal worstAbsCdfErr;
		mpfr::mpreal  worstRelCdfErr;
		mpfr::mpreal worstAbsPmfErr;
		mpfr::mpreal  worstRelPmfErr;
		mpfr::mpreal rmsAbsCdfErr;
		mpfr::mpreal rmsRelCdfErr;
		mpfr::mpreal rmsAbsPmfErr;
		mpfr::mpreal rmsRelPmfErr;
	};
	
template<class T>
void DumpCdfStats(std::string prefix, std::ostream &dest_file, typename DiscreteDistribution<T>::TypePtr got, typename DiscreteDistribution<T>::TypePtr target, T start_sigma, std::ostream *dump)
{
	T one=start_sigma/start_sigma;
	
	if(!target->IsSymmetric())
		throw std::string("DumpStats - Currently only strictly symmetric target distributions are supported.");
	if(target->StandardMoment(1)!=0)
		throw std::string("DumpStats - Currently only zero mean distributions are supported.");
	
	if(!got->IsSymmetric())
		std::cerr<<"WARNING : Asymmetric distribution. This isn't handled explicitly (or well) right now.\n";
	if(got->StandardMoment(1)!=0)
		std::cerr<<"WARNING : Non zero-mean distribution. This isn't handled explicitly (or well) right now.\n";
	
	std::vector<segment> segments;
	
	int64_t curr_index=target->IndexFromRange(0)+1;
	
	T worstAbsCdfErr=0, worstRelCdfErr=0;
	T worstAbsPmfErr=0, worstRelPmfErr=0;
	T sumSqrAbsCdfErr=0, sumSqrRelCdfErr=0;
	T sumSqrAbsPmfErr=0, sumSqrRelPmfErr=0;
	T count=0;
	
	T sumAbsCdfErr=0, sumRelCdfErr=0;
	
	if(dump){
		(*dump)<<"x, sigma, target_cdf, got_cdf, cdf_abs_err, cdf_rel_err, target_pdf, got_pdf\n";
	}
	
	T target_sigma=target->StandardMoment(2);
	
	bool currently_contiguous=true;
	T max_contiguous=0.0;
	
	for(double sigma=-0.5;sigma>=start_sigma;sigma-=0.5){
		int64_t begin_index=target->ClosestIndexFromRange(sigma * target_sigma);
		int64_t end_index=curr_index;
		int64_t n=end_index-begin_index;
	
		std::vector<T> cdf_target(n), range_target(n), pmf_target(n);
	
		target->RangeByIndex(begin_index, end_index, &range_target[0]);
		target->CdfByIndex(begin_index, end_index, &cdf_target[0]);
		target->PmfByIndex(begin_index, end_index, &pmf_target[0]);
		
		for(int i=n-1;i>=0;i--){
			T cdf=got->Cdf(range_target[i]);
			T pmf=got->Pmf(range_target[i]);
			
			if(pmf==0){
				currently_contiguous=false;
			}else if(currently_contiguous){
				max_contiguous=range_target[i];
			}
			
			assert(cdf_target[i].get_prec()>=one.get_prec());
			assert(cdf.get_prec()>=one.get_prec());
			assert(pmf.get_prec()>=one.get_prec());
			
			//std::cerr<<"prec(p)="<<p.get_prec()<<", prec(range)="<<range_target[i].get_prec()<<", prec(target)="<<cdf_target[i].get_prec()<<"\n";
			
			T absCdfErr=abs(cdf-cdf_target[i]);
			T relCdfErr=absCdfErr/cdf_target[i];
			
			T absPmfErr=abs(pmf-pmf_target[i]);
			T relPmfErr=absPmfErr/pmf_target[i];
			
			if(dump){
				(*dump)<<range_target[i]<<", "<<range_target[i]*target_sigma<<", "<<cdf_target[i]<<", "<<cdf<<", "<<absCdfErr<<", "<<relCdfErr<<",   "<<pmf_target[i]<<", "<<pmf<<"\n";
			}
			
			worstAbsCdfErr=std::max(worstAbsCdfErr, absCdfErr);
			worstRelCdfErr=std::max(worstRelCdfErr, relCdfErr);
			
			worstAbsPmfErr=std::max(worstAbsPmfErr, absPmfErr);
			worstRelPmfErr=std::max(worstRelPmfErr, relPmfErr);
			
			sumSqrAbsCdfErr+= absCdfErr*absCdfErr;
			sumSqrRelCdfErr+= relCdfErr*relCdfErr;
			
			sumSqrAbsPmfErr+= absPmfErr*absPmfErr;
			sumSqrRelPmfErr+= relPmfErr*relPmfErr;
		}
		count+=n;
		
		segment seg;
		seg.sigma=sigma;
		seg.start_x=target->RangeFromIndex(begin_index);
		seg.worstAbsCdfErr=worstAbsCdfErr;
		seg.worstRelCdfErr=worstRelCdfErr;
		seg.worstAbsPmfErr=worstAbsPmfErr;
		seg.worstRelPmfErr=worstRelPmfErr;
		
		seg.rmsAbsCdfErr=sqrt(sumSqrAbsCdfErr/count);
		seg.rmsRelCdfErr=sqrt(sumSqrRelCdfErr/count);
		
		seg.rmsAbsPmfErr=sqrt(sumSqrAbsPmfErr/count);
		seg.rmsRelPmfErr=sqrt(sumSqrRelPmfErr/count);
		
		segments.push_back(seg);
		
		curr_index=begin_index;
	}
	
	for(int i=0;i<(int)segments.size();i++){
		dest_file<<prefix<<"MaxAbsCdfErr, "<<segments[i].sigma<<", "<<segments[i].worstAbsCdfErr<<"\n";
		dest_file<<prefix<<"MaxRelCdfErr, "<<segments[i].sigma<<", "<<segments[i].worstRelCdfErr<<"\n";
		dest_file<<prefix<<"MaxAbsPmfErr, "<<segments[i].sigma<<", "<<segments[i].worstAbsPmfErr<<"\n";
		dest_file<<prefix<<"MaxRelPmfErr, "<<segments[i].sigma<<", "<<segments[i].worstRelPmfErr<<"\n";
		
		dest_file<<prefix<<"RmsAbsCdfErr, "<<segments[i].sigma<<", "<<segments[i].rmsAbsCdfErr<<"\n";
		dest_file<<prefix<<"RmsRelCdfErr, "<<segments[i].sigma<<", "<<segments[i].rmsRelCdfErr<<"\n";
		dest_file<<prefix<<"RmsAbsPmfErr, "<<segments[i].sigma<<", "<<segments[i].rmsAbsPmfErr<<"\n";
		dest_file<<prefix<<"RmsRelPmfErr, "<<segments[i].sigma<<", "<<segments[i].rmsRelPmfErr<<"\n";
	}
	
	dest_file<<prefix<<"MaxContiguous, - ,"<<max_contiguous<<"\n";
	dest_file<<prefix<<"MaxContiguousSigma, - ,"<<max_contiguous*target_sigma<<"\n";
	dest_file<<prefix<<"MaxSupport, - ,"<<got->Support().first<<"\n";
	dest_file<<prefix<<"MaxSupportSigma, - ,"<<got->Support().first*target_sigma<<"\n";
}

template<class T>
void DumpMomentStats(std::string prefix, std::ostream &dest_file, typename DiscreteDistribution<T>::TypePtr got, typename DiscreteDistribution<T>::TypePtr target)
{
	for(int i=2;i<=8;i+=2){
		T e=target->StandardMoment(i), o=got->StandardMoment(i);
		dest_file<<prefix<<"MomentTarget, "<<i<<", "<<e<<"\n";
		dest_file<<prefix<<"MomentGot, "<<i<<", "<<o<<"\n";
		dest_file<<prefix<<"MomentRelErr, "<<i<<", "<<(o-e)/e<<"\n";
	}
}

struct rng_mpfr
{
	int m_prec;
	 gmp_randstate_t m_rng;
	
	rng_mpfr(int prec)
		: m_prec(prec)
	{
		gmp_randinit_default(m_rng);
	}
		
	~rng_mpfr()
	{
		gmp_randclear(m_rng);
	}
	
	mpfr::mpreal operator()()
	{ 
		mpfr::mpreal res(0, m_prec);
		mpfr_urandomb(res.mpfr_ptr(), m_rng);
		return res;
	}
};

template<class T>
void DumpChi2Stats(std::string prefix, std::ostream &dest_file, typename DiscreteDistribution<T>::TypePtr got, typename DiscreteDistribution<T>::TypePtr target)
{	
	T zero=target->Pmf(0);
	zero=zero-zero;
	T one=zero+1;
	
	rng_mpfr rng(one.get_prec());
	
	for(int k=4;k<=8192;k*=2){
		boost::shared_ptr<const Chi2Partition<T> > part(Chi2Partition<T>::CreateEqualProbability(target, k));
		bool fail=false;
		double lastGood=0;
		for(double log2n=8;log2n<=64;log2n+=0.5){
			T n=ceil(pow(2*one, log2n));
			boost::shared_ptr<const Chi2Partition<T> > part_adapt(part->AdaptForSampleSize(n));
			boost::shared_ptr<const Chi2Partition<T> > part_got(part_adapt->Resample(got));
			
			std::vector<T> sample=part_got->GenerateSample(n, rng);
			Chi2Res res=part_adapt->Chi2Test(sample);
			
			if(::flopoco::verbose >= DEBUG)
				std::cerr<<"k="<<k<<", n=2^"<<log2n<<", stat="<<res.statistic<<", pvalue="<<res.pvalue<<"\n";
			if(res.pvalue >= 0.01){
				lastGood=log2n;
			}
			if(res.pvalue < 1e-6){
				fail=true;
				break;
			}
		}
		dest_file<<prefix<<"Chi2EqualProbFail, "<<k<<", "<<(fail?"True":"False")<<"\n";
		dest_file<<prefix<<"Chi2EqualProbLastGoodLog2, "<<k<<", "<<lastGood<<"\n";
	}
	
	for(int k=4;k<=8192;k*=2){
		boost::shared_ptr<const Chi2Partition<T> > part(Chi2Partition<T>::CreateEqualRange(target, -16, +16, k));
		bool fail=false;
		double lastGood=0;
		for(double log2n=8;log2n<=64;log2n+=0.5){
			T n=ceil(pow(2*one, log2n));
			boost::shared_ptr<const Chi2Partition<T> > part_adapt(part->AdaptForSampleSize(n));
			boost::shared_ptr<const Chi2Partition<T> > part_got(part_adapt->Resample(got));
			
			std::vector<T> sample=part_got->GenerateSample(n, rng);
			Chi2Res res=part_adapt->Chi2Test(sample);
			
			if(::flopoco::verbose >= DEBUG)
				std::cerr<<"k="<<k<<", n=2^"<<log2n<<", stat="<<res.statistic<<", pvalue="<<res.pvalue<<"\n";
			if(res.pvalue >= 0.01){
				lastGood=log2n;
			}
			if(res.pvalue < 1e-6){
				fail=true;
				break;
			}
		}
		dest_file<<prefix<<"Chi2EqualRangeFail, "<<k<<", "<<(fail?"True":"False")<<"\n";
		dest_file<<prefix<<"Chi2EqualRangeLastGoodLog2, "<<k<<", "<<lastGood<<"\n";
	}
}

static void TransformStatsUsage(std::ostream &dst)
{
	OperatorFactory::classic_OP(dst, "TransformStats", "[-dump file] [-prec bits] [-startsigma sigma] distribution", false);
	dst << "       Calculates the accuracy of the previous statistical transform operator.\n";
	dst <<"          distribution - Specifies the distribution to check against.\n";
	dst<<"    distributions are:\n";
	dst<<"      QuantisedGaussian[stddev,fracbits]\n";
	dst<<"      DiscreteGaussian[stddev,fracbits]\n";
	dst<<"    optional '-dump' argument will dump distribution to the given file.\n";
	dst<<"    optional '-prec' argument will set number of bits for calculations.\n";
	dst<<"    optional '-startsigma' argument determines where distribution is examined from (should be negative).\n";
	dst<<"    optional '-prefix' argument specifies a string to come before output lines.\n";
}

DiscreteDistribution<mpfr::mpreal>::TypePtr ParseDistribution(std::string spec, int prec)
{
	std::vector<std::string> parts;
	boost::split(parts, spec, boost::is_any_of("[,]"));
	
	if(parts.size()==0)
		throw std::string("ParseDistribution('")+spec+"') - Doesn't contain any parts.";
	
	//for(int i=0;i<(int)parts.size();i++){
	//	std::cerr<<i<<" : '"<<parts[i]<<"', len="<<parts[i].size()<<"\n";
	//}
	
	if(parts[0]=="QuantisedGaussian" || parts[0]=="DiscreteGaussian"){
		if(parts.size()<3)
			throw std::string("ParseDistribution('")+spec+"') - Need exactly three parts for this distribution.";
		mpfr::mpreal stddev(0, prec);
		parseSollyaConstant(stddev.mpfr_ptr(), parts[1], MPFR_RNDN);
		int fracbits=boost::lexical_cast<int>(parts[2]);
		
		assert(stddev.get_prec()==prec);
		
		if(parts[0]=="QuantisedGaussian"){
			return boost::make_shared<QuantisedGaussianDistribution<mpfr::mpreal> >(stddev, fracbits);
		}else{
			return boost::make_shared<DiscreteGaussianDistribution<mpfr::mpreal> >(stddev, fracbits);
		}
	}else{
		throw std::string("ParseDistribution('")+spec+"') - Didn't understand first part '"+parts[0]+"'.";
	}
}

void DumpRngInfo(std::string prefix, std::ostream &dest, RngTransformOperator *op,DiscreteDistribution<mpfr::mpreal>::TypePtr gotDist, DiscreteDistribution<mpfr::mpreal>::TypePtr targetDist)
{
	dest<<prefix<<"NumOutputs, - , "<<op->nonUniformOutputCount()<<"\n";
	dest<<prefix<<"HomogenousOutputs, - , "<<op->nonUniformOutputCount()<<"\n";
	dest<<prefix<<"OutputWidth, 0 , "<<op->nonUniformOutputWidth(0)<<"\n";
	dest<<prefix<<"OutputFracBits, 0 , "<<log2(targetDist->RangeGranularity())<<"\n";
	mpfr::mpreal stddevOverRes=targetDist->StandardMoment(2) / targetDist->RangeGranularity();
	dest<<prefix<<"StdDev/Granularity, - , "<<stddevOverRes<<"\n";
	dest<<prefix<<"log2(StdDev/Granularity), - , "<<log2(stddevOverRes)<<"\n";
}

static Operator *TransformStatsParser(Target *target ,const std::vector<std::string> &args,int &consumed)
{
	consumed=0;
	
	std::ostream *dest_stream=&std::cout;
	std::ofstream dest_file;
	std::ofstream dump_stream;
	int prec=128;
	std::string startsigma_str="-12";
	std::string prefix;
	
	while(args.size()-consumed > 0){
		if(args[consumed]=="-dest"){
			consumed++;
			if(args.size()<2){
				throw std::string("TransformStatsPaser - Need at least three arguments if -dest option is used.");
			}
			std::string filename=args[consumed++];
			
			if(::flopoco::verbose > 0)
				std::cerr<<"TransformStats - Writing stats to file '"<<filename<<"'.\n";
			
			dest_file.open(filename.c_str());
			
			if(!dest_file.is_open())
				throw std::string("TranformStatsParser - Couldn't open dump target file.");
			
			dest_stream=&dest_file;
		}else if(args[consumed]=="-dump"){
			consumed++;
			if(args.size()<2){
				throw std::string("TransformStatsPaser - Need at least three arguments if -dump option is used.");
			}
			std::string filename=args[consumed++];
			
			if(::flopoco::verbose > 0)
				std::cerr<<"TransformStats - Dumping distribution to file '"<<filename<<"'.\n";
			
			dump_stream.open(filename.c_str());
			
			if(!dump_stream.is_open())
				throw std::string("TranformStatsParser - Couldn't open dump target file.");
			
			dump_stream.precision(16);
		}else if(args[consumed]=="-prec"){
			consumed++;
			
			prec=atoi(args[consumed++].c_str());
			if(::flopoco::verbose > 0)
				std::cerr<<"TransformStats - Working precision set to "<<prec<<" bits.\n";
		}else if(args[consumed]=="-startsigma"){
			consumed++;
			
			startsigma_str=args[consumed++];
		}else if(args[consumed]=="-prefix"){
			consumed++;
			
			prefix=args[consumed++];
		}else{
			break;
		}
	}
	
	if(args.size()-consumed <1)
		throw std::string("TransformStatsParser - No distribution specified.");
	
	if(::flopoco::verbose>=INFO)
		std::cerr<<prefix<<"Getting target distribution\n";
	DiscreteDistribution<mpfr::mpreal>::TypePtr targetDist=ParseDistribution(args[consumed++], prec);
	
	mpfr::mpreal startsigma(0, prec);
	parseSollyaConstant(startsigma.mpfr_ptr(), startsigma_str, MPFR_RNDN);
	
	vector<Operator*> *ops=target->getGlobalOpListRef();
	if(ops->size()==0)
		throw std::string("TransformStatsParser - No previous operator to work with.");
	
	Operator *op=ops->back();
	
	RngTransformOperator *transform=dynamic_cast<RngTransformOperator*>(op);
	if(transform==NULL)
		throw std::string("TransformStatsParser - Previous operator is not an RngTransformOperator.");
	
	IRngTransformDistributions *opDists=dynamic_cast<IRngTransformDistributions*>(op);
	if(opDists==NULL)
		throw std::string("TransformStatsParser - Previous operator does not support IRngTransformDistributions.");
	
	if(!transform->nonUniformOutputsAreHomogenous() && transform->nonUniformOutputCount() > 1)
		throw std::string("TransformStatsParser - Can't deal with multiple non-homogeneous output distributions.");
	
	if(::flopoco::verbose>=INFO)
		std::cerr<<prefix<<"Getting RNG's CDF.\n";
	Distribution<mpfr::mpreal>::TypePtr got=opDists->nonUniformOutputDistribution(0, prec);
	DiscreteDistribution<mpfr::mpreal>::TypePtr gotDist=boost::dynamic_pointer_cast<DiscreteDistribution<mpfr::mpreal> >(got);
	assert(gotDist->Cdf(0).get_prec()>=prec);
	
	int pp=dest_stream->precision();
	dest_stream->precision(16);
	if(::flopoco::verbose>=INFO)
		std::cerr<<prefix<<"Calculating Cdf Stats\n";
	DumpCdfStats(prefix, *dest_stream, gotDist, targetDist, startsigma, dump_stream.is_open() ? &dump_stream : NULL);
	if(::flopoco::verbose>=INFO)
		std::cerr<<prefix<<"Calculating moment Stats\n";
	DumpMomentStats<mpfr::mpreal>(prefix, *dest_stream, gotDist, targetDist);
	if(::flopoco::verbose>=INFO)
		std::cerr<<prefix<<"Calculating Chi2 Stats\n";
	DumpChi2Stats<mpfr::mpreal>(prefix, *dest_stream, gotDist, targetDist);
	DumpRngInfo(prefix, *dest_stream, transform, gotDist, targetDist);
	dest_stream->precision(pp);
	if(::flopoco::verbose>=INFO)
		std::cerr<<prefix<<"Finished stats.\n";
	
	return NULL;
}

void TransformStats_registerFactory()
{
	DefaultOperatorFactory::Register(
		"transform_stats",
		"utility",
		TransformStatsUsage,
		TransformStatsParser
	);
}

};
};
