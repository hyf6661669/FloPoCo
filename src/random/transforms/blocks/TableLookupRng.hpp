/* Each Operator declared within the flopoco framework has 
to inherit the class Operator and overload some functions listed below*/
#ifndef random_table_lookup_transform_hpp
#define random_table_lookup_transform_hpp
#include <vector>
#include <sstream>
#include <string>
#include <fstream>
#include "gmp.h"
#include "mpfr.h"
#include <math.h>
#include <assert.h>
#include <iomanip>
#include <time.h>

#include "RngTransformOperator.hpp"

/* This file contains a lot of useful functions to manipulate vhdl */
#include "utils.hpp"



/*  All flopoco operators and utility functions are declared within
  the flopoco namespace.
*/
namespace flopoco{


// new operator class declaration
class bitwise : public RngTransformOperator {
  public:

    static string operatorInfo;
    int MSB; 
    int LSB;
    std::string m;

    int n;
    double lambda;
    std::vector<mpz_class> p;		//global variable to keep p & q arrays
    std::vector<int> r;
    unsigned int j;

//m vector used in CPP
    std::vector<int> mVector;
    
//No. of uniform bits
    int uni_no;

  public:
    // definition of some function for the operator    

static mpz_class conv_f2z(mpfr_t a, int fw);
static void conv_z2f(mpfr_t b, mpz_class a, int fw);
static int MaxElement( vector<int> mi, int size)
{
			int max=mi.at(0);
			for(int i=1; i< size; i++)
			{
				if(mi.at(i)>max)	max=mi.at(i);
			}

	return (max);

}

//=======================================================================================================================

static vector<mpz_class> CalculatePiForMarray(int MSB, int LSB, vector<int> mi, double lambda)
{
			int width=MSB+LSB+1;
			vector<mpz_class> p(width);
			mpfr_t pi, u, b;
			int max=MaxElement( mi, width);

				mpfr_init2 (pi, max+32);
				mpfr_init2 (u, max+32);
				mpfr_init2 (b, max+32);

	  			mpfr_set_d (u, 1.0, GMP_RNDD);
	  			mpfr_set_d (b, 2.0, GMP_RNDD);

//pi 0-LSB		max-MSB
			
			//work out values for p
			for (int i=-LSB, j=0; i <= MSB; i++, j++)
			{
				mpfr_pow_si(pi, b, i, GMP_RNDN);	//2^i
				mpfr_mul_d(pi, pi, lambda, GMP_RNDN);	//lambda*2^i
				mpfr_exp(pi, pi, GMP_RNDN);
				mpfr_add(pi, pi, u, GMP_RNDN);
				mpfr_div(pi, u, pi, GMP_RNDN);
				p.at(j) = bitwise::conv_f2z(pi, mi.at(width-j-1));		//p[j] mpz_class value
			}
		  mpfr_clears(pi, u, b, NULL);

		return(p);
}

//=======================================================================================================================
static void ReverseSearch(double tau, int MSB, int LSB, int m_max, double lambda, vector<unsigned long> &order, vector<int> &m)
{
	assert(MSB>=(-LSB));
	int i=0, set=0;
	int width=MSB+LSB+1;

//hold the current i for testing of error<tau
	int j=0;
	unsigned long length=(1<<width);
		std::vector<double> errorarray;
		errorarray.resize(length);

		std::vector<mpz_class>  p,q;			
		p.resize(width);
		q.resize(width);

                ofstream fm;
                fm.open("reverse.txt");

		stringstream fnamee;
		fnamee << "error_2_"<< MSB << "_" << LSB <<"_"<<lambda << "_" << tau*100000000 << ".txt";
		ofstream ferror;
		ferror.open(fnamee.str().c_str());

		fm << "i-> " << i;
		for(int q=0; q<width; q++)
		{fm << "\nm[" << q << "]-> " << m[q];}
		fm << "\n=======================================";

//================================================================================================================
	while(i<width)
	{	
		if(i>j) j=i;
		printf("i=%d	j=%d\n", i,j);
		fm << "\nj-> " << j << "	i-> "<< i;
// always test for the current i position
		if (bitwise::DoesMeetTau(tau, i, j, MSB, LSB, lambda, order, m, p, q, errorarray)){
			if(set==0){
				m[i]=0; set=1;		//(1)
			}
			else{
				i++; set=0;		//(2)
				//print m values

				for(int q=0; q<width; q++)
				{fm << "\nm[" << q << "]-> " << m[q];}
				fm << "\n=======================================";
			}

		}else{
			if (m[i]==m_max && i==0){	//(5)

			std::cerr << " CalculateMarrayFor Tau: fail with MSB reaches m_max, need more integer bits" << endl;
			exit (EXIT_FAILURE);//error message		back to MSB, the end. or add more integer bits
			
			}else if (m[i]==m_max){		//(3)
			i--;
			set=1;
				
			}else{				//(4)
				m[i]++;
			}
		}

	}
//================================================================================================================
//second optimisation: output uniform random number directly when pi=0.5
		//int max=MaxElement( mi, width);
		for(int i=0; i<width; i++)
		{
			mpfr_t p_t;
			//mpfr_init2(p_t, max+32);
			bitwise::conv_z2f(p_t, p.at(i), m.at(i));
			if(mpfr_cmp_d(p_t, 0.5)==0)	m.at(i)=0;
		}


		fm << "\n==========================================\n";
		fm << setw(8)<<"tau->" << setw(8)<< tau << "\n";
		fm <<setw(8)<<"lambda->" << setw(8)<< lambda << "\n";
		fm << setw(8)<< "MSB->" << setw(8)<<MSB <<setw(8)<<"LSB-> " <<setw(8)<< -LSB << "\n";
		fm <<setw(8)<<"width->" <<setw(8)<< width << "\n";

		fm << setw(8)<<"m_max->" <<setw(8)<< m_max << "\n";
		int twidth=0;
		for(int q=0; q<width; q++)
		{twidth+=m.at(q);}

		fm <<setw(8)<<"totalw->" <<setw(8)<< twidth <<"\n";

		fm<<flush;
		fm.close();

			for(unsigned long i=0; i<length; i++)
			{ferror<< errorarray.at(i)<< endl;}
			ferror<< flush;
			ferror.close();

}


//=======================================================================================================================
static bool DoesMeetTau(double tau, int m_MSB, int m_LSB, int MSB, int LSB, double lambda, vector<unsigned long> &order, vector<int> &m, vector<mpz_class> &p, vector<mpz_class> &q, vector<double> &errorarray)
{
			//DoesMeetTau(tau, i, j, MSB, LSB, lambda, order, m, p, q))
assert(m_MSB<=m_LSB);

bool result=true;
			
int width=MSB+LSB+1;

unsigned long lmax=(1<<width);

vector<mpz_class> Ecdf(lmax+1, 0);


		mpz_class xx, xt, pdf_t, pdf, cdf;
		mpfr_t u, b, error, temp, temp2, xf, cdf_1,cdf_2;		//cdf_2:ideal	cdf_1:empirical


		double errord=0.0;

//total Uniform random number needed
		int totalw;

		totalw=0;
		for(int i=0; i<width; i++ )			//m array, 0->MSB	max->LSB
		{totalw= totalw+m.at(i);}


			mpfr_init2 (xf, LSB);
			mpfr_init2 (error, totalw);
			mpfr_init2 (temp, totalw);
		        mpfr_init2 (temp2, totalw);
			mpfr_init2 (cdf_1, totalw);
			mpfr_init2 (cdf_2, totalw);


			//calculate empirical discrete cdf of each iteration
			//pr[order.at(y)] include the current order.at(y)

			mpfr_set_d (error, 0.0, GMP_RNDD);

			//function call,calculate p values
			p=bitwise::CalculatePiForMarray(MSB, LSB, m, lambda);


			//calculate q values
			for (int j=0; j<width; j++)
			{
			    mpz_ui_pow_ui(q.at(j).get_mpz_t(), 2 , m.at(width-j-1));
			    q.at(j)=q.at(j)-p.at(j);			//q[j]=2^m-p[j]
			}


			for (unsigned long i=0; i<lmax; i++)
			{	pdf=1;
				mpz_set_ui(xt.get_mpz_t(),i);

				for(int j=0; j<width; j++)
				{
					if(mpz_tstbit(xt.get_mpz_t(),j)==1) {mpz_mul(pdf_t.get_mpz_t(), pdf.get_mpz_t(), p.at(j).get_mpz_t());}
					else	{mpz_mul(pdf_t.get_mpz_t(),pdf.get_mpz_t(),q.at(j).get_mpz_t());}
					pdf=pdf_t;
				}
				Ecdf.at(i+1)=Ecdf.at(i)+pdf;

			}


//================================================================================================================================================
//reverse detection algorithm
			int y=(1<<m_MSB)-1;	//index traces order array elements
			for(int g=m_MSB; (g<=m_LSB) && result; g++)		//g traces mi array elements
			{
							//initialisation with floating point precision based on m
							mpfr_init2 (u, m.at(g)+32);
							mpfr_init2 (b, m.at(g)+32);

				  			mpfr_set_d (u, 1.0, GMP_RNDD);
				  			mpfr_set_d (b, 2.0, GMP_RNDD);



			//============================================================================				

				for(unsigned long n=0; (n<(unsigned long)(1<<g)) && result; n++, y++)
				{	
				// calculate empirical CDF


							bitwise::conv_z2f(cdf_1, Ecdf.at(order.at(y)+1), totalw);		// practical value

				//ideal cdf value
							mpz_set_ui(xx.get_mpz_t(), order.at(y));			// ideal value
							bitwise::conv_z2f(xf, xx, LSB);
							mpfr_mul_d(xf, xf, lambda, GMP_RNDN);		//lambda*x
							mpfr_neg(xf, xf, GMP_RNDN);			//-xf

							mpfr_exp(temp, xf, GMP_RNDN);			// exp(-xf)
							mpfr_sub(cdf_2, u, temp, GMP_RNDN);		//1-exp(-xf)

							//mpfr_printf("cdf_2-> %Rf\n", cdf_2);

							mpfr_sub(temp, cdf_1, cdf_2, GMP_RNDN);
							mpfr_abs(error, temp, GMP_RNDN);			//|temp|


					if(mpfr_cmp_d(error,tau)>=0)	result=false;

					else {
						//print out error values
						errord=mpfr_get_d(error, GMP_RNDN);
						errorarray.at(order.at(y))=errord;
						}
			
			}

			//============================================================================	

				printf("\ng-> %d\n", g);
		}

		//call free the space
		  mpfr_clears(u, b, error, xf, temp, temp2, cdf_1,cdf_2, NULL);

	return(result);

}



//=======================================================================================================================
	static vector<int> CalculateMarrayForTau(int MSB, int LSB, double tau, double lambda)	//array of different values of m in regards of each random bit
	{

 //operation time measurement
		clock_t a=clock();
// define m according to the Random number quality tau
			int width=MSB+LSB+1;


//calculate m_max & intialise mi array
			int m_max;
			m_max=floor(-log2(tau))+10;
			vector<int> mi(width,m_max);

			unsigned long xmax;	
			xmax=(1<<width)-1;	


		stringstream fname;
		fname << "summary_2_"<< MSB << "_" << LSB <<"_"<<lambda << "_" << tau*100000000<< ".txt";

//output report file
		ofstream fout;
		fout.open(fname.str().c_str());

//generate the correct order of x to be processed
		unsigned long v=(1<<width)-2;
		unsigned long j=1, s=0;		//j index of order. s index of previous array


		vector<unsigned long> order(xmax);
		order.at(0)=(v>>1);	//first element

		for(int i=1; i<width; i++)
		{
			for(unsigned long n=0; n<((unsigned long)(1<<(i-1)));n++,j++,s++)
			{
				order.at(j)=(order.at(s)-1)>>1;
			}

			for(unsigned long n=0; n<((unsigned long)(1<<(i-1)));n++,j++,s++)
			{
				order.at(j)=v-order.at(s);
			}

			s=s-(1<<(i-1));
		}


//call reverse search algorithm
	bitwise::ReverseSearch( tau, MSB, LSB, m_max, lambda, order, mi);

//execution time measurement
	clock_t bb=clock();

	ofstream fc;
        fc.open("time_fr_elpase.txt",ios::app);
        fc<< (bb-a) <<endl;
        fc<< flush;
        fc.close();

//=======================================================================================================================
// data printing
		fout << "============================================\n";
		fout << setw(8)<< "m_max->" << setw(8) << m_max<< endl;

		fout << setw(8)<< "MSB->" << setw(8)<<MSB << endl;
		fout <<setw(8)<<"LSB-> " <<setw(8)<< -LSB << endl;


		int twidth=0;
				for(int q=0; q<width; q++)
				{twidth+=mi.at(q);}
		
		fout <<setw(8)<<"totalw->" <<setw(8)<< twidth <<"\n";


		fout <<"\n===========================================\nM array\n";

		for(int i=0; i<width;i++)
		{
			fout<< i<< "->	"<<mi.at(i)<< "\n";
		}
		fout <<"\n===========================================\n";


		fout<<flush;
		fout.close();



		return (mi);

	}

//=======================================================================================================================

	
    // constructor, defined there with two parameters (default value 0 for each)
    	bitwise(Target* target,int MSB, int LSB, string m, double lambda=1.0);		//private?


    // destructor
    	~bitwise();

	static bitwise *instance;		

	static bitwise *CreateOperatorForTolerance(Target *target, int MSB, int LSB, double tau, double lambda)
	{
		vector <int> mj(MSB+LSB+1);
		mj=CalculateMarrayForTau(MSB, LSB, tau, lambda);
		

		std::string m; 
		stringstream ss;

		for(int i=0; i<MSB+LSB+1; i++)
		{ ss << mj[i] <<" ";}

		m=ss.str();

		cout<< m << "\n";

		instance = new bitwise(target, MSB, LSB, m, lambda);

		return (instance);
	}

	int GetUniformInputs(void)		//total uniform bits
	{
		int totalw=0;
		//mj=bitwise::CalculateMarrayForTau(MSB, LSB, tau, lambda);

		for(int i=0; i<n; i++ )			//mi array, 0->MSB	max->LSB
		{totalw+=mVector[i];}

		return (totalw);
	}



    void emulate(TestCase * tc);

    void buildRandomTestCases(TestCaseList* tcl, int n);

    TestCase* buildRandomTestCase(int i);
};

};
#endif
