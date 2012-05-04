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
#include "bitwise.hpp"

using namespace std;

namespace flopoco
{
namespace random
{
	
	extern vector<Operator *> oplist;

	// personalized parameter
	string bitwise::operatorInfo = "bitwise info MSB LSB m lambda";

 	bitwise * bitwise::instance=0;	

//=======================================================================================================================
	mpz_class BitwiseExpTransform::conv_f2z(mpfr_t a, int fw)
	{
				mpz_class b;
			
				mpfr_mul_2si( a, a, fw, GMP_RNDN);    //left shift fw bits
				mpfr_get_z( b.get_mpz_t(), a, GMP_RNDD);

				return (b);
	}

//=======================================================================================================================
	void BitwiseExpTransform::conv_z2f(mpfr_t b, mpz_class a, int fw)
	{
				mpfr_init2(b,fw+32);

				mpfr_set_z(b, a.get_mpz_t(), GMP_RNDN);	
				mpfr_div_2si(b, b, fw, GMP_RNDN);    //rignt shift fw bits			
	}


	BitwiseExpTransform::BitwiseExpTransform(Target* target, int MSB, int LSB, std::vector<int> m, double lambda)
	: Operator(target), MSB(MSB), LSB(LSB), lambda(lambda)
	, mVector(m)
	{

  		// definition of the name of the operator
  		ostringstream name;
  		name << "bitwise_" << MSB << "_" << LSB << "_" << lambda;
  		setName(name.str());
  		// Copyright 
  		setCopyrightString("Junfei Yan 2011");

		//printf("\nm=%d \n", m);
		n= MSB+LSB+1;
		p.resize(n);
		
		if(mVector.size()==1){
			while(mVector.size()<n)
				mVector.push_back(mVector.front());
		}else if mVector.size()!=n){
			REPORT(ERROR, "m input is neither a scalar, nor the same size as the output.");
			exit(1);
		}
		
		uni_no=GetUniformInputs();

  		// declaring inputs
  		addInput ("U" , uni_no, true);		//n = MSB+LSB+1  	width of result

  		// declaring output
  		addOutput("R" , n, 1, true);

		p=bitwise::CalculatePiForMarray(MSB, LSB, mVector, lambda);
		

	        int rangel=0, rangeh=mVector[n-1]-1;

		for (int i=-LSB, j=0; i <= MSB; i++, j++)
		{
		//start from LSB to MSB	
			if (mVector[n-j-1]==0){
				vhdl << tab <<  "R" << of(j) << "<= '0';" <<endl;
			}else
				vhdl << tab <<  "R" << of(j) << "<= '1' WHEN U"<< range(rangeh, rangel) << "< (" << "\"" << unsignedBinary(p[j], mVector[n-j-1])<<"\") ELSE '0';" << endl;
			}

			if(j<n-1)
			{
				rangel+=mVector[n-j-1];
				rangeh=rangel+mVector[n-j-2]-1;
			}

		}

	}


BitwiseExpTransform::~BitwiseExpTransform(){
}


void BitwiseExpTransform::emulate(TestCase * tc) {

  std::cerr<<">emulate\n";

  mpz_class su = tc->getInputValue("U");
  mpz_class sui;
  /* then manipulate our bit vectors in order to get the correct output*/
  mpz_class sr;	


  int max=bitwise::MaxElement( mVector, n);

/* Compute correct value */
		mpfr_t sum, temp, two;
		mpfr_init2(sum, max+32);
		mpfr_init2(temp, max+32);
		mpfr_init2(two, max+32);
		mpfr_set_d(sum, 0.0, GMP_RNDD);
		mpfr_set_d(two, 2.0, GMP_RNDD);

//keep the result random bits
		r.resize(n);
		std::cerr<<" 	b array -> ";
		mpz_class mask;

		int len=0;

//comparison starts form MSB to LSB
		for (int i=-LSB, j=0; i <=MSB; i++, j++)
		{
                        mpz_ui_pow_ui(mask.get_mpz_t(), 2, mVector[n-j-1]);		//different mask values for each bits
                        mask = mask-1;
			mpz_fdiv_q_2exp(sui.get_mpz_t(), su.get_mpz_t(), len);		//su right shift len bits
			mpz_and(sui.get_mpz_t(), mask.get_mpz_t(), sui.get_mpz_t());	//mask the high order bits

			if (p[j] > sui) {   r[j] = 1;	
					    mpfr_pow_si(temp, two, i, GMP_RNDN);
					    mpfr_add(sum, sum, temp, GMP_RNDN);}

				else 	{   r[j] = 0;}

			len+=mVector[n-j-1];			//update No. of shift bits
		}

		sr = bitwise::conv_f2z(sum, LSB);

//output sum


  /* at the end, we indicate to the TestCase object what is the expected
    output corresponding to the inputs */
  tc->addExpectedOutput("R",sr);
  mpfr_clears(sum, temp, two, NULL);
}



TestCase* BitwiseExpTransform::buildRandomTestCase(int i) {

		TestCase* tc = new TestCase(this);
		mpz_class a = getLargeRandom(uni_no);
		tc->addInput("U", a);
			
  	/* Get correct outputs */
  	emulate(tc);

  	return tc;
}

};
};
