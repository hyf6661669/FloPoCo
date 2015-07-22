#ifndef FIXSIF_HPP
#define FIXSIF_HPP

#include "Operator.hpp"
#include "utils.hpp"
#include "FixSOPC.hpp"

#include <boost/numeric/ublas/matrix.hpp>

namespace flopoco{ 

	
	/** SIF building 
	  */
	class FixSIF : public Operator {
	  
	public:
		/** normal constructor, building the SIF out of the specificatin file */
		FixSIF(Target* target, int lsb, std::string file,  std::map<std::string, double> inputDelays = emptyDelayMap); 

		/* Destructor */
		~FixSIF();

		/** The method that does the bulk of operator construction, isolated to enable sub-classes such as FixHalfSine etc */
		void buildVHDL();

		// Below all the functions needed to test the operator
		/* the emulate function is used to simulate in software the operator
			 in order to compare this result with those outputed by the vhdl opertator */
		void emulate(TestCase * tc);

		//int p;							/**< The precision (opposite of LSB weight) of inputs and outputs */ 
		//TODO: update definitions to msb+lsb
		int n;							/**< number of taps */

		int wO;

	private:
		std::string file;
		mpfr_t mpcoeff[10000];			/**< the absolute values of the coefficients as MPFR numbers */
		bool coeffsign[10000];			/**< the signs of the coefficients */

		mpfr_t *xPrec; //history of xs used by emulate (only 2 )
		mpfr_t *xCurrent; //history of xs used by emulate (only 2 )
		//mpz_class xHistory[10000]; // history of x used by emulate
		//int currentIndexB;
		//mpfr_t yHistory[10000]; // history of y (result) used by emulate

		std::vector < std::vector<std::string> > coeffs;			/**< the coefficients as strings */
		std::vector <FixSOPC*> sopcs; //list of sopcs (corresponding to matrix lines)

		std::vector <int> msbIn, lsbIn, msbOut, lsbOut, msbInter, lsbInter;

		uint32_t nt; /**< number of intermediate computations */
		uint32_t nx; /**< size of the state-space */
		uint32_t ny; /**< number of outputs */
		uint32_t nu; /**< numnber of inputs */
		uint32_t ne=20; //number of taps TODO: calculate this

		/** parseFile parses files corresponding to the following format:
			" header
			  matrix
			"
				with header="matrixname lines columns"
					example: "J 0 0"
					with matrixname="{J,K,L,M,N,P,Q,R,S}"
						 lines=* (integer)
						 columns=* (integer)
					 matrix=lines and columns corresponding to the coefficients of the matrix (integers or floating point format)
					 see examples in /SIFexamples
			Precisions for inputs and outputs have to be specified in an other file.
		 */
		int parseFile(); //fills nt, nx, nu, ny and coeffs out of the input file. Returns 0 if succeed, error else.

		/**
			computeMSBsLSBs: old function to read precisions in a file
			output:
				-0 if success
		  */
		int computeMSBsLSBs();

	int computeABCD(boost::numeric::ublas::matrix<double> const &bJ, boost::numeric::ublas::matrix<double> const &bK, boost::numeric::ublas::matrix<double> const &bL, boost::numeric::ublas::matrix<double> const &bM, boost::numeric::ublas::matrix<double> const &bN, boost::numeric::ublas::matrix<double> const &bP, boost::numeric::ublas::matrix<double> const &bQ, boost::numeric::ublas::matrix<double> const &bR, boost::numeric::ublas::matrix<double> const &bS, boost::numeric::ublas::matrix<double> &bA, boost::numeric::ublas::matrix<double> &bB, boost::numeric::ublas::matrix<double> &bC, boost::numeric::ublas::matrix<double> &bD);
		/**
			readPrecision: reads precisions in a file
			input:
				-msbsIn: the msbs of inputs for the operator
				-lsbsIn: the lsbs of inputs for the operator
				-msbsOut: the desired msbs of outputs for the operator
				-lsbsOut: the desired lsbs of outputs for the operator
				-inFile: 1 to read in the file "precisions.txt", 0 for taking default precision
				TODO: give the ability to user- define the precision file
			output:
				-0 if success
				-1 if error
				TODO: handle error propagation in the parsing framework
		  */
		int readPrecision( std::vector <int> &msbsIn, std::vector<int> &lsbsIn, std::vector<int> &msbsOut, std::vector<int> &lsbsOut, bool inFile=1 );

		/**
			bMToDoubleM: fills the matrix doubleM from the coefficients stored in bM.
			Note: this assumes that doubleM is instanciated
			input: 
				-bM: boost matrix to convert in an array-matrix j
				-doubleM: the array to fill
			Ouptut: 0 if success.
		*/
		int bMToDoubleM( boost::numeric::ublas::matrix<double> const &bM, double * &doubleM );

		/**
			vvToBoostMatrix: fills the matrix bM from the string coefficients stored in sM.
			Note: this assumes that bM is instanciated
			Input:
				-sM: std::vector of std::vector of strings (it won't be modified)
				-bM: boost matrix to be filled
			Output:
				-0 if success
		  */
		int vvToBoostMatrix( std::vector< std::vector <std::string> > const &sM, boost::numeric::ublas::matrix<double> &bM );


		int readMatrices( std::vector< std::vector < std::vector<std::string> >**> &Z, std::ifstream &openedFile, int &lc);

		int readMatrix(std::string &header, std::string JKLMNPQRS, std::vector < std::vector <std::string> > * &toFill, std::ifstream &openedFile, int lc = 0 );

		/**
			displayBMatrix: displays the matrix M using the macro REPORT
		  */
		void displayBMatrix( boost::numeric::ublas::matrix<double> const &M , std::string const name ) const;


	};

}

#endif
