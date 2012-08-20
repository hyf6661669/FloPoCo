#include "optimise_rounding.hpp"

#include "fixed_point_exp_tester.hpp"

#include "func_approx_exp_stage.hpp"
#include "close_table_exp_stage.hpp"
#include "find_close_value.hpp"
#include "chained_exp_stage.hpp"
#include "multiplier_exp_stage.hpp"
#include "table_exp_stage.hpp"

#include "residual_type.hpp"
#include "result_type.hpp"
#include "fixed_point_exp_stage.hpp"

#include "Targets/Virtex4.hpp"
#include <fstream>

#define BOOST_TEST_MODULE FixedPointExp
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(TestOptimiseRounding)
{
	::flopoco::verbose=DEBUG;
	
	typedef flopoco::random::result_type <double> result_t;
	
	typedef flopoco::random::residual_type<double> residual_t;
	typedef flopoco::random::result_type <double> result_t;
	
	boost::shared_ptr<flopoco::Target> target(new flopoco::Virtex4());
	
	residual_t inputResidual(true, 24, 6);
	
	boost::shared_ptr<flopoco::random::CloseTableExpStage> hi(new flopoco::random::CloseTableExpStage(target.get(),
		0.0, 1.0, // mu, sigma
		inputResidual, 
		6, // table address bits
		-pow(2.0, -24), pow(2.0, -24)
	));
	boost::shared_ptr<flopoco::random::CloseTableExpStage> lo(new flopoco::random::CloseTableExpStage(target.get(),
		0.0, 1.0, // mu, sigma
		hi->OutputResidualType(),
		6,
		-pow(2.0, -24), pow(2.0, -24)
	));
	
	flopoco::random::multiplier_rounding_spec<double> res=flopoco::random::FindResultBitsForMultiply(
		hi->OutputResultType().FracWidth(),
		hi->GetAllResultValues(),
		lo->OutputResultType().FracWidth(),
		lo->GetAllResultValues(),
		-pow(2.0,-23),
		pow(2.0,-23)
	);
	
	std::cerr<<"   Output residual = "<<lo->OutputResidualType()<<"\n";
}

BOOST_AUTO_TEST_CASE(TestFuncApproxExpStage)
{
	typedef flopoco::random::residual_type<double> residual_t;
	typedef flopoco::random::result_type <double> result_t;
	
	boost::shared_ptr<flopoco::Target> target(new flopoco::Virtex4());
	
	residual_t inputResidual(false, 4, -6);
	
	boost::shared_ptr<flopoco::random::FuncApproxExpStage> t(new flopoco::random::FuncApproxExpStage(target.get(),
		0.0, 1.0, // mu, sigma
		inputResidual,
		12
	));
	
	boost::shared_ptr<flopoco::random::FixedPointExpTester> tt(new flopoco::random::FixedPointExpTester(target.get(),
		t
	));
}


BOOST_AUTO_TEST_CASE(FuncApproxExpStage)
{
	typedef flopoco::random::residual_type<double> residual_t;
	typedef flopoco::random::result_type <double> result_t;
	
	boost::shared_ptr<flopoco::Target> target(new flopoco::Virtex4());
	
	residual_t inputResidual(false, 4, -6);
	
	boost::shared_ptr<flopoco::random::FuncApproxExpStage> t(new flopoco::random::FuncApproxExpStage(target.get(),
		0.0, 1.0, // mu, sigma
		inputResidual,
		16
	));
	
	std::ofstream dst("test_FuncApproxExpStage.vhdl");
	t->outputVHDLToFile(dst);
}

BOOST_AUTO_TEST_CASE(CloseTableExpStage)
{
	typedef flopoco::random::residual_type<double> residual_t;
	typedef flopoco::random::result_type <double> result_t;
	
	boost::shared_ptr<flopoco::Target> target(new flopoco::Virtex4());
	
	residual_t inputResidual(true, 8, 4);
	
	boost::shared_ptr<flopoco::random::CloseTableExpStage> t(new flopoco::random::CloseTableExpStage(target.get(),
		0.0, 1.0, // mu, sigma
		inputResidual, 
		4, // table address bits
		0.0, pow(2.0,-25)
	));
	
	std::ofstream dst("test_CloseTableExpStage.vhdl");
	t->outputVHDLToFile(dst);
}

/*
BOOST_AUTO_TEST_CASE(FindCloseValue)
{
	typedef flopoco::random::residual_type<double> residual_t;
	
	residual_t orig(true, 24, 4);
	std::cerr<<"Type="<<orig<<", Width="<<orig.Width()<<"\n";
	residual_t msbsA=orig.take_msbs(6);  
	residual_t offsetA=orig.drop_msbs(6);
	double minErr=0, maxErr=pow(2.0,-26);
	
	typedef flopoco::random::close_value_table_t<double> res_t;
	
	res_t resA=FindCloseValues(0.0, 1.0, msbsA, offsetA, minErr, maxErr);
	
	residual_t msbsB=orig.drop_msbs(5).take_msbs(6);  
	residual_t offsetB=orig.drop_msbs(5).drop_msbs(6);
	
	std::cerr<<"fracBitsA="<<resA.resultSpace.fracWidth<<"\n";
	
	res_t resB=FindCloseValues(0.0, 1.0, msbsB, offsetB, minErr, maxErr);
	
	std::cerr<<"fracBitsB="<<resB.resultSpace.fracWidth<<"\n";
	
	std::cerr<<"residual="<<orig.drop_msbs(10)<<"\n";
}
*/

BOOST_AUTO_TEST_CASE(ChainedExpStage)
{
	typedef flopoco::random::residual_type<double> residual_t;
	typedef flopoco::random::result_type <double> result_t;
	
	boost::shared_ptr<flopoco::Target> target(new flopoco::Virtex4());
	
	residual_t inputResidual(true, 8, 4);
	
	boost::shared_ptr<flopoco::random::TableExpStage> s1(new flopoco::random::TableExpStage(target.get(),
		0.0, 1.0, // mu, sigma
		inputResidual, 
		6, // table address bits
		16 	// output frac bits
	));
	boost::shared_ptr<flopoco::random::TableExpStage> s2(new flopoco::random::TableExpStage(target.get(),
		0.0, 1.0, // mu, sigma
		s1->OutputResidualType(),
		s1->OutputResidualType().Width(),
		16	// output frac bits
	));
	
	boost::shared_ptr<flopoco::random::MultiplierExpStage> m2(
		new flopoco::random::MultiplierExpStage(target.get(),
		s2,
		s1->OutputResultType(),
		12
	));
	
	std::vector<flopoco::random::FixedPointExpStagePtr> stages;
	stages.push_back(s1);
	stages.push_back(m2);
	
	boost::shared_ptr<flopoco::random::ChainedExpStage> c(
		new flopoco::random::ChainedExpStage(target.get(),
		stages
	));
	
	std::ofstream dst("test_ChainedExpStage.vhdl");
	s1->outputVHDLToFile(dst);
	s2->outputVHDLToFile(dst);
	m2->outputVHDLToFile(dst);
	c->outputVHDLToFile(dst);
}

BOOST_AUTO_TEST_CASE(MultiplierExpStage)
{
	typedef flopoco::random::residual_type<double> residual_t;
	typedef flopoco::random::result_type <double> result_t;
	
	boost::shared_ptr<flopoco::Target> target(new flopoco::Virtex4());
	
	residual_t inputResidual(true, 8, 4);
	
	boost::shared_ptr<flopoco::random::TableExpStage> s1(new flopoco::random::TableExpStage(target.get(),
		0.0, 1.0, // mu, sigma
		inputResidual, 
		6, // table address bits
		16 	// output frac bits
	));
	boost::shared_ptr<flopoco::random::TableExpStage> s2(new flopoco::random::TableExpStage(target.get(),
		0.0, 1.0, // mu, sigma
		s1->OutputResidualType(),
		s1->OutputResidualType().Width(),
		16	// output frac bits
	));
	
	boost::shared_ptr<flopoco::random::MultiplierExpStage> m(
		new flopoco::random::MultiplierExpStage(target.get(),
		s2,
		s1->OutputResultType(),
		12
	));
	
	std::ofstream dst("test_MultiplierExpStage.vhdl");
	m->outputVHDLToFile(dst);
}

BOOST_AUTO_TEST_CASE(TableExpStage)
{
	typedef flopoco::random::residual_type<double> residual_t;
	typedef flopoco::random::result_type <double> result_t;
	
	boost::shared_ptr<flopoco::Target> target(new flopoco::Virtex4());
	
	residual_t inputResidual(true, 8, 4);
	
	boost::shared_ptr<flopoco::random::TableExpStage> t(new flopoco::random::TableExpStage(target.get(),
		0.0, 1.0, // mu, sigma
		inputResidual, 
		6, // table address bits
		8 	// output frac bits
	));
	
	std::ofstream dst("test_TableExpStage.vhdl");
	t->outputVHDLToFile(dst);
}

BOOST_AUTO_TEST_CASE(ResidualTypeTests)
{
	typedef flopoco::random::residual_type<double> residual_t;
	
	flopoco::random::residual_type<double> r(false, 8);
	
	BOOST_CHECK_EQUAL(r.Width(), 8);
	BOOST_CHECK_EQUAL(r.MsbPower(), -1);
	BOOST_CHECK_EQUAL(r.LsbPower(), -1-7);
	BOOST_CHECK_EQUAL(r.TypeMin(),0.0);
	BOOST_CHECK_EQUAL(r.TypeMax(),255.0/256);
	BOOST_CHECK_EQUAL(r.RangeMin(),r.TypeMin());
	BOOST_CHECK_EQUAL(r.RangeMax(),r.TypeMax());
	BOOST_CHECK_EQUAL(r.Round(1.0/1024),0);
	BOOST_CHECK_EQUAL(r.Round(3.0/1024),1.0/256);
	
	r=flopoco::random::residual_type<double>(false, 8, -2);
	
	BOOST_CHECK_EQUAL(r.Width(), 8);
	BOOST_CHECK_EQUAL(r.MsbPower(), -3);
	BOOST_CHECK_EQUAL(r.LsbPower(), -3-7);
	BOOST_CHECK_EQUAL(r.TypeMin(),0);
	BOOST_CHECK_EQUAL(r.TypeMax(),255.0/1024);
	BOOST_CHECK_EQUAL(r.RangeMin(),r.TypeMin());
	BOOST_CHECK_EQUAL(r.RangeMax(),r.TypeMax());
	BOOST_CHECK_EQUAL(r.Round(1.0/4096),0);
	BOOST_CHECK_EQUAL(r.Round(3.0/4096),1.0/1024);
	
	r=flopoco::random::residual_type<double>(true, 6, 0);
	
	BOOST_CHECK_EQUAL(r.Width(), 7);
	BOOST_CHECK_EQUAL(r.MsbPower(), -1);
	BOOST_CHECK_EQUAL(r.LsbPower(), -1-5);
	BOOST_CHECK_EQUAL(r.TypeMin(),-1);
	BOOST_CHECK_EQUAL(r.TypeMax(),63.0/64);
	BOOST_CHECK_EQUAL(r.RangeMin(),r.TypeMin());
	BOOST_CHECK_EQUAL(r.RangeMax(),r.TypeMax());
	BOOST_CHECK_EQUAL(r.Round(1.0/256),0);
	BOOST_CHECK_EQUAL(r.Round(3.0/256),1.0/64);
	
	BOOST_CHECK_EQUAL(r.take_msbs(0), residual_t(false, 0, 0));
	BOOST_CHECK_EQUAL(r.drop_msbs(0), residual_t(true, 6, 0));
	
	BOOST_CHECK_EQUAL(r.take_msbs(1), residual_t(true, 0, 0));
	BOOST_CHECK_EQUAL(r.drop_msbs(1), residual_t(false, 6, 0));
	
	BOOST_CHECK_EQUAL(r.take_msbs(2), residual_t(true, 1, 0));
	BOOST_CHECK_EQUAL(r.drop_msbs(2), residual_t(false, 5, -1));
	
	BOOST_CHECK_EQUAL(r.take_msbs(5), residual_t(true, 4, 0));
	BOOST_CHECK_EQUAL(r.drop_msbs(5), residual_t(false, 2, -4));
	
	r=flopoco::random::residual_type<double>(false, 8, 8);
	
	BOOST_CHECK_EQUAL(r.take_msbs(0),  residual_t(false, 0, 8) );
	BOOST_CHECK_EQUAL(r.drop_msbs(0),  residual_t(false, 8, 8) );
	
	BOOST_CHECK_EQUAL(r.take_msbs(1),  residual_t(false, 1, 8) );
	BOOST_CHECK_EQUAL(r.drop_msbs(1),  residual_t(false, 7, 7) );
	
	BOOST_CHECK_EQUAL(r.take_msbs(7),  residual_t(false, 7, 8) );
	BOOST_CHECK_EQUAL(r.drop_msbs(7),  residual_t(false, 1, 1) );
	
	BOOST_CHECK_EQUAL(r.take_msbs(8),  residual_t(false, 8, 8) );
	BOOST_CHECK_EQUAL(r.drop_msbs(8),  residual_t(false, 0, 0) );
}


BOOST_AUTO_TEST_CASE(ResultTypeTests)
{
	flopoco::random::result_type<double> r(7);
	
	BOOST_CHECK_EQUAL(r.Round(1.5), 1.5);
	BOOST_CHECK_EQUAL(r.Round(1.0), 1.0);
	BOOST_CHECK_EQUAL(r.Round(255.0), 255.0);
	BOOST_CHECK_EQUAL(r.Round(255.1), 255.0);
	BOOST_CHECK_EQUAL(r.Round(249.9), 250.0);
	
	BOOST_CHECK_EQUAL(r.Width() , 7);
	
	
	r.Add(1.5);
	BOOST_CHECK_EQUAL(r.Width() , 7);
	BOOST_CHECK_EQUAL(r.expMin , 1);
	BOOST_CHECK_EQUAL(r.expMax , 1);
	BOOST_CHECK_EQUAL(r.fracMin , 0.75);
	BOOST_CHECK_EQUAL(r.fracMax , 0.75);
	BOOST_CHECK_EQUAL(r.ToBits(1.5),0x40);
	BOOST_CHECK_EQUAL(r.FromBits(r.ToBits(1.5)), 1.5);
	
	r.Add(1.75);
	BOOST_CHECK_EQUAL(r.Width() , 7);
	BOOST_CHECK_EQUAL(r.expMin , 1);
	BOOST_CHECK_EQUAL(r.expMax , 1);
	BOOST_CHECK_EQUAL(r.fracMin , 0.75);
	BOOST_CHECK_EQUAL(r.fracMax , 0.875);
	BOOST_CHECK_EQUAL(r.ToBits(1.75),0x60);
	
	BOOST_CHECK_EQUAL(r.FromBits(r.ToBits(1.5)), 1.5);
	BOOST_CHECK_EQUAL(r.FromBits(r.ToBits(1.75)), 1.75);
	
	r.Add(0.5);
	BOOST_CHECK_EQUAL(r.Width() , 8);
	BOOST_CHECK_EQUAL(r.expMin , 0);
	BOOST_CHECK_EQUAL(r.expMax , 1);
	BOOST_CHECK_EQUAL(r.fracMin , 0.5);
	BOOST_CHECK_EQUAL(r.fracMax , 0.875);
	BOOST_CHECK_EQUAL(r.ToBits(0.5),0);
	BOOST_CHECK_EQUAL(r.ToBits(255.0/256.0),127);
	BOOST_CHECK_EQUAL(r.ToBits(510.0/256.0),(1<<7) | 127);
	
	BOOST_CHECK_EQUAL(r.FromBits(r.ToBits(0.5)), 0.5);
	BOOST_CHECK_EQUAL(r.FromBits(r.ToBits(1.5)), 1.5);
	BOOST_CHECK_EQUAL(r.FromBits(r.ToBits(1.75)), 1.75);
	
	r=flopoco::random::result_type<double>(7);
	
	r.Add(0.5);
	BOOST_CHECK_EQUAL(r.FracWidthNonZero(), 0);
	
	r.Add(0.5+1.0/256);
	BOOST_CHECK_EQUAL(r.FracWidthNonZero(), 1);
	
	r.Add(0.5+8.0/256+1.0/256);
	BOOST_CHECK_EQUAL(r.FracWidthNonZero(), 4);
	
	r.Add(255/256.0);
	BOOST_CHECK_EQUAL(r.FracWidthNonZero(), r.FracWidth());
}

