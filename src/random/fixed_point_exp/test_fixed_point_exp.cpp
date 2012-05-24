#include "result_type.hpp"

#define BOOST_TEST_MODULE FixedPointExp
#include <boost/test/unit_test.hpp>

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
}
