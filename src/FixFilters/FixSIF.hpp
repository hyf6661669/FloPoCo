#ifndef FIXSIF_HPP
#define FIXSIF_HPP

#include "Operator.hpp"
#include "utils.hpp"

namespace flopoco{ 

	
	class FixSIF : public Operator {
	  
	public:
		/** normal constructor, building the SIF out of the specificatin file */
		FixSIF(Target* target, int lsb, string file,  map<string, double> inputDelays = emptyDelayMap); 

		/* Destructor */
		~FixSIF();

		/** The method that does the bulk of operator construction, isolated to enable sub-classes such as FixHalfSine etc */
		void buildVHDL();

		// Below all the functions needed to test the operator
		/* the emulate function is used to simulate in software the operator
			 in order to compare this result with those outputed by the vhdl opertator */
		void emulate(TestCase * tc);

		int p;							/**< The precision (opposite of LSB weight) of inputs and outputs */ 
		//TODO: update definitions to msb+lsb
		int n;							/**< number of taps */

		int wO;

	private:
		string file;
		mpfr_t mpcoeff[10000];			/**< the absolute values of the coefficients as MPFR numbers */
		bool coeffsign[10000];			/**< the signs of the coefficients */

		vector<mpz_class*> xHistories; //history of xs used by emulate
		vector<mpfr_t*> yHistories; //history of ys used by emulate
		//mpz_class xHistory[10000]; // history of x used by emulate
		//int currentIndexA;
		//int currentIndexB;
		//mpfr_t yHistory[10000]; // history of y (result) used by emulate

		vector < vector<string> > coeffs;			/**< the coefficients as strings */
		int nt; /**< number of intermediate computations */
		int nx; /**< size of the state-space */
		int ny; /**< number of outputs */
		int nu; /**< numnber of inputs */
		int ne=20; //number of taps TODO: calculate this

		int parseFile(); //fills nt, nx, nu, ny and coeffs out of the input file. Returns 0 if succeed, error else.
		int readMatrices( vector< vector <vector<string> >**> &Z, ifstream &openedFile, int &lc);
		int readMatrix(string &header, string JKLMNPQRS, vector < vector <string> > * &toFill, ifstream &openedFile, int lc = 0 );
		int readPrecision( vector <int> &msbs, vector<int> &lsbs, bool inFile=1 );
	};

}

#endif
