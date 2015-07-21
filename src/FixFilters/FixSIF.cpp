#include <iostream>
#include <sstream>
#include <queue>

#include "gmp.h"
#include "mpfr.h"
#include "sollya.h"
#include "soplex.h"

#include "FixSIF.hpp"

#include "ShiftReg.hpp"

#include <boost/numeric/ublas/lu.hpp> 
#include <boost/numeric/ublas/io.hpp>
extern "C" {
#include "/usr/local/include/wcpg_API.h"
}
//#include "/usr/local/include/clapack_headers_config.h"	

#define BUF_SIZE 256 //default buffer size
#define DEFAULT_LSB -64 //default precision for sopcs
#define DEFAULT_MSB 8 //default precision for sopcs

//using namespace std;
namespace ublas = boost::numeric::ublas; 
template<class T> 
bool invertMatrix (const ublas::matrix<T>& input, ublas::matrix<T>& inverse) { 
	std::cout<<"entering invertMatrix"<<std::endl;
	typedef ublas::permutation_matrix<std::size_t> pmatrix; 
	// create a working copy of the input 
	ublas::matrix<T> A(input); 
	// create a permutation matrix for the LU-factorization 
	pmatrix pm(A.size1()); 
	// perform LU-factorization 
	int res = lu_factorize(A,pm); 
	if( res != 0 ){
		std::cout<<"exiting invertMatrix by short path"<<std::endl;
		return false; 
	}
	// create identity matrix of "inverse" 
	inverse.assign(ublas::identity_matrix<T>(A.size1())); 
	// backsubstitute to get the inverse 
	lu_substitute(A, pm, inverse); 
	std::cout<<"exiting invertMatrix by long path"<<std::endl;
	return true; 
}
bool matMul (ublas::matrix<double> &res, ublas::matrix<double> const &X, ublas::matrix<double> const &Y) { 
	if(X.size2()!=Y.size1())
		return false;

	for (size_t i=1; i<X.size1(); i++) {
		for (size_t j=1; j<X.size1(); j++) {
			res(i,j)=0.0;
		}
	}
	
	for (size_t i=1; i<X.size1(); i++) {
		for (size_t j=1; j<Y.size2(); j++) {
			for (size_t k=1; k<X.size2(); k++) {
				res(i,j) += X(i,k) * Y(k,j);
			}
		}
	}
	return true;
}

bool matAdd (ublas::matrix<double> &res, ublas::matrix<double> const &X, ublas::matrix<double> const &Y) { 
	for (size_t i=1; i<X.size1(); i++) {
		for (size_t j=1; j<X.size1(); j++) {
			res(i,j) = X(i,j) + Y(i,j);
		}
	}
	return true;
}

//double maxi( double *wcpg, int l, int c, int cs) {
//	
//}

	int getFullLine( std::ifstream &f, std::string &s, int &lc ){
		std::cout<<"entered getFullLine"<<std::endl;

			getline(f, s);
			lc++;
			//check for empty lines
			if ( s.empty() ) {
				while ( s.empty() && (!f.eof()))
				getline(f, s);
				lc++;
			}

			//check for empty lines filled with spaces or other non senses
			for (uint32_t i = 0; i<s.size(); i++) {
				if ( isdigit(s[i]) ){
		std::cout<<"exiting getFullLine by short path"<<std::endl;
					return 0;
				}

				if ((!isdigit(s[i])) && (i==s.size()-1) && (!f.eof()) ) {
					getline(f, s);
					lc++;
					i=-1;
				}
			}
		std::cout<<"exiting getFullLine by long path"<<std::endl;
			return 1;
	}

//TODO: build a policy to detect suitable use of shiftRegisters.
namespace flopoco {
	const int veryLargePrec = 6400;  /*6400 bits should be enough for anybody */

	int displayMatrix( std::vector< std::vector<std::string> > &m, std::string name = "matrix"){
		std::string srcFileName="FixSIF";
		//srcFileName<<"TOTO.cpp";
		REPORT(0, "Displaying Matrix "<<name<<" "<<m.size()<<" "<<m[0].size());
		for (uint32_t i=0; i<m.size(); i++){
			std::stringstream line;
			for (uint32_t  j=0; j<m[i].size(); j++){
				line<<m[i][j];
				if (j!=m[i].size()-1){
					line<<"	";
				}
			}
			REPORT(0,line.str());
		}
		return 0;
	}


	FixSIF::FixSIF(Target* target, int lsb_, std::string file_, std::map<std::string, double> inputDelays) : 
		Operator(target), file(file_)
	{
		target->setNotPipelined();

		srcFileName="FixSIF";
		setCopyrightString ( "Antoine Martinet, Florent de Dinechin (2014)" );

		std::ostringstream name;
		//name << "FixSIF_" << p << "_" << coeffs.size() << "_uid" << getNewUId();
		name << "FixSIF_uid" << getNewUId();
		setNameWithFreq( name.str() );

		parseFile();

//		int msbOut=1;//TODO: calculate this
//		int g=-128;//TODO: calculate this
//
//		int hugePrec = 10*(1+msbOut+p+g);
//
//		for (int i = 0; i<nx; i++){
//			xHistories.push_back( (mpz_class*)malloc(sizeof(mpz_class)*10000));
//		}
//
//		for (int i = 0; i<ny; i++){
//			yHistories.push_back( (mpfr_t*)malloc(sizeof(mpz_class)*10000));
//		}
//
//		for (int i = 0; i<ny; i++){
//			for (int j = 0; j<ne; j++)
//			{
//				mpfr_init2 (yHistories[i][j], hugePrec);
//				mpfr_set_d (yHistories[i][j], 0.0, GMP_RNDN);
//			}
//		}
//
//		for(int i=0; i<nx; i++) {
//			for(int j=0; j<ne; j++) {
//				xHistories[i][j]=0;
//			}
//		}

		buildVHDL();
	};

	

	// The method that does the work once coeffs[][] is known
	void FixSIF::buildVHDL(){
//	//TODO: beginning of new code
		displayMatrix(coeffs, "coeffs");
		n = nt + nx + nu; //TODO: update this definition with a clever number
		readPrecision(msbIn, lsbIn, msbOut, lsbOut, 0);
		int precCompBool=computeMSBsLSBs();
		if (precCompBool) {
			THROWERROR("Unknown error occured during precision computation");
		}

		useNumericStd_Unsigned();
		//if(p<1) {
		//	THROWERROR("Can't build an architecture for this value of LSB")
		//}

		// initialize stuff for emulate
		//xPrec= (mpfr_t *) malloc ( nx * sizeof(mpfr_t*));
		xPrec= new mpfr_t[nx];
		//xCurrent= (mpfr_t *) malloc ( nx * sizeof(mpfr_t*));
		xCurrent= new mpfr_t[nx];
		for (uint32_t i=0; i<nx; i++){
			mpfr_init2( xPrec[i], veryLargePrec );
			mpfr_init2( xCurrent[i], veryLargePrec );
			mpfr_set_ui( xPrec[i], 0, GMP_RNDN );
			mpfr_set_ui( xCurrent[i], 0, GMP_RNDN );
		}

		REPORT(0,"Emulate init done");

		//declare intermediate T(k+1)
		for (uint32_t i=0; i<nt; i++) {
			declare(join("T",i), 1-lsbInter[i]);
		}

		//nu inputs numbering from 0 to nu-1 FIXME: check the relevance of this convention dealing with SIF
		for (uint32_t i=0; i<nu; i++) {
			std::ostringstream input;
			input << "U" <<i;
			addFixInput(std::string(input.str()), true, 1+msbIn[i], lsbIn[i]);
			REPORT(0,"Popopo"<<i);
		}

		//ny outputs numbering from 0 to ny-1 FIXME: check the relevance of this convention dealing with SIFs
		for (uint32_t i=0; i<ny; i++) {
			std::ostringstream output;
			//output << "Y";// << i;
			std::string name=join("Y",i);
			addFixOutput(name, 1, 1+msbOut[i], lsbOut[i]);
		}
	
		//declare intermediate X and X(k+1)
		for (uint32_t i=0; i<nx; i++) {
			//declare(join("Xplus",i));
			declare(join("X",i), 1-lsbInter[nt+i]);
		}

		std::vector < std::vector <std::string> > coefDel(nt+nx+ny); //copy of coeffs which will be deleted.
		//FIXME: refactor algorithm with iterators not to use copies.
		for ( uint32_t i=0; i<coeffs.size(); i++ ) {
			for ( uint32_t j=0; j<coeffs.size(); j++ ){
				coefDel[i].push_back(coeffs[i][j]);
			}
		}

		while ( !coefDel.empty() ) {

			std::vector <std::string> nonZeros(n);
			std::queue <int> coefIndst;
			std::queue <int> coefIndsx;
			std::queue <int> coefIndsu;

			//copy first line of coeffs
				REPORT(0,"non zeros coefficients number="<<nonZeros.size());
			for (uint32_t i = 0; i<coefDel[0].size(); i++){
				nonZeros[i] = coefDel[0][i];
				REPORT(0,"coeffs["<<i<<"]="<<coefDel[0][i]);
			}
			//output matrix std::vector
			//three steps (t, x and u) because all are not treated the same way
			//drop zeros and (implicit) ones in the diagonal but keep indices for consistency
			int ib = 0; //bias in std::vector indices induced by erase operations
			uint32_t i = 0;
			std::vector <int> lsbInC, msbInC;
			for ( ; i<nt; i++ ){
				if ( ( std::stof(nonZeros[i+ib].c_str()) == 0.0 ) || ( (std::stof(nonZeros[i+ib].c_str()) == 1.0)&&(n-coefDel.size()==i)) ){
					REPORT(0,"no coef found at "<<i+ib);
					REPORT(0,"numerical value of non zero coeff at "<<i+ib<<":"<<std::stof(nonZeros[i+ib].c_str()));
					nonZeros.erase(nonZeros.begin()+i+ib);
					ib--;
				}
				else {
					coefIndst.push(i);
					lsbInC.push_back(lsbIn[i]);
					msbInC.push_back(msbIn[i]);
					REPORT(0,"coef found at "<<i+ib);
					REPORT(0,"numerical value of non zero coeff at "<<i+ib<<":"<<std::stof(nonZeros[i+ib].c_str()));
					REPORT(0,"pushing "<<i<<" in list of indices for T");
				}
			}

					REPORT(0,"number of non zeros coefficients:"<<nonZeros.size()<<", bias="<<ib<<", iteration="<<i<<", n=nu+nx+nt="<<nu+nx+nt);
			
			//only drop zeros for x
			for ( ; i<nt+nx; i++ ){
				if ( std::stof(nonZeros[i+ib].c_str()) == 0.0 ){
					REPORT(0,"no coef found at "<<i+ib);
					REPORT(0,"numerical value of non zero coeff at "<<i+ib<<":"<<std::stof(nonZeros[i+ib].c_str()));
					nonZeros.erase(nonZeros.begin()+i+ib);
					ib--;
				}
				else{
					coefIndsx.push(i-nt);
					lsbInC.push_back(lsbIn[i]);
					msbInC.push_back(msbIn[i]);
					REPORT(0,"coef found at "<<i+ib);
					REPORT(0,"numerical value of non zero coeff at "<<i+ib<<":"<<std::stof(nonZeros[i+ib].c_str()));
					REPORT(0, "pushing "<<i-nt<<" in list of indices for x");
				}
			}

			REPORT(0, "cleaning u from implicit zeros, bias="<<ib<<", iteration="<<i<<", n=nu+nx+nt="<<nu+nx+nt);

			//and for u
			for ( ; i<nu+nx+nt; i++ ){
				if ( std::stof(nonZeros[i+ib].c_str()) == 0.0 ){
					REPORT(0,"no coef found at "<<i+ib);
					REPORT(0,"numerical value of non zero coeff at "<<i+ib<<":"<<std::stof(nonZeros[i+ib].c_str()));
					nonZeros.erase(nonZeros.begin()+i+ib);
					ib--;
				}
				else{
					coefIndsu.push(i-nt-nx);
					lsbInC.push_back(lsbIn[i]);
					msbInC.push_back(msbIn[i]);
					REPORT(0,"coef found at "<<i+ib);
					REPORT(0,"numerical value of non zero coeff at "<<i+ib<<":"<<std::stof(nonZeros[i+ib].c_str()));
					REPORT(0, "pushing "<<i-nt<<" in list of indices for x");
				}
			}

			//if only one coeff, just write a multiplier, if only ones, just write an adder
			//if ( nonZeros.size() > 1 ) 
					
					REPORT(0, "Calling FixSOPC with arguments:");
					REPORT(0, "	target="<<target_);
					for ( unsigned int deb=0; deb<msbIn.size(); deb++ ) {
						REPORT(0, "	msbIn["<<deb<<"]="<<msbIn[deb]);
						REPORT(0, "	lsbIn["<<deb<<"]="<<lsbIn[deb]);
					}
					REPORT(0, "	sizeof lsbIn="<<lsbIn.size());
					REPORT(0, "	msbOut="<<*msbOut.begin());
					REPORT(0, "	lsbOut="<<*lsbOut.begin());
					for ( unsigned int deb=0; deb<nonZeros.size(); deb++ ) {
						REPORT(DEBUG, "	nonZeros["<<deb<<"]="<<nonZeros[deb]);
					}
					sopcs.push_back(new FixSOPC( target_, msbInC, lsbInC, *msbInter.begin(), *lsbInter.begin(), nonZeros));

					addSubComponent(*sopcs.rbegin());
					msbOut.erase(msbInter.begin());
					lsbOut.erase(lsbInter.begin());

					//wiring
					i=0;
					for ( ; i<coefIndst.size(); i++) {
					REPORT(0, "inporting T, with parameters X"<<i<<" and T"<<coefIndst.front());
						inPortMap(*sopcs.rbegin(), join("X", i), join("T", coefIndst.front()));
						coefIndst.pop();
					}

					for ( ; i-coefIndst.size()<coefIndsx.size(); i++) {
					REPORT(0, "inporting X, with parameters X"<<i<<" and X"<<coefIndsx.front());
						inPortMap(*sopcs.rbegin(), join("X", i), join("X", coefIndsx.front()));
						coefIndsx.pop();
					}
					for ( ; i-( coefIndst.size() + coefIndsx.size() )<coefIndsu.size(); i++) {
					REPORT(0, "inporting U, with parameters X"<<i<<" and U"<<coefIndsu.front());
						inPortMap(*sopcs.rbegin(), join("X", i), join("U", coefIndsu.front()));
						coefIndsu.pop();
					}
					REPORT(0, "inport std::map done");

					if (sopcs.size() <= nt) {
					REPORT(0, "outporting T, with parameters R"<<" and T"<<sopcs.size()-1);
						outPortMap(*sopcs.rbegin(), "R", join("T", sopcs.size()-1));
						vhdl << instance(*sopcs.rbegin(), join("fixSOPC_t",sopcs.size()-1));
					}
					else if (sopcs.size()<=nt+nx){
					REPORT(0, "outporting Xplus, with parameters R"<<" and Xplus"<<sopcs.size()-nt);
						outPortMap(*sopcs.rbegin(), "R", join("Xplus", sopcs.size()-nt));
						vhdl << instance(*sopcs.rbegin(), join("fixSOPC_y",sopcs.size()));
					}
					else if (sopcs.size()<=nt+nx+ny){
					REPORT(0, "outporting Y, with parameters R"<<" and Y"<<sopcs.size()-nt-nx);
						outPortMap(*sopcs.rbegin(), "R", join("Y", sopcs.size()-(nt+nx)),0);
						vhdl << instance(*sopcs.rbegin(), join("fixSOPC_y",sopcs.size()));
					}
					else {
						THROWERROR("Unexpected error happened during internal building");
						}
					REPORT(0, "outport std::map done");
					


					coefDel.erase(coefDel.begin());
					
					REPORT(0, "finished line. Beginning step "<<coefDel.size()<<", coeffs.empty()="<<coefDel.empty());
		}
		for (uint32_t i = 0; i < nx; i++){
			std::stringstream sig;
			sig<<"Xplus"<<i;
//			REPORT(0,"export 	X"<<i<<" <= "<<delay(sig.str(),1)<<";");
//			vhdl << tab << join("X",i)<<" <= "<<delay(sig.str(),1)<<";"<<endl;
		}
		REPORT(0, "time shifting done");
	//TODO: end of new code
		
	};

	FixSIF::~FixSIF(){

	};

	void FixSIF::emulate(TestCase * tc){

		std::cout<<"	toto is ok -6"<<std::endl;
		std::vector<mpz_class> uValues(nu);
		std::cout<<"	toto is ok -5"<<std::endl;
		std::vector<mpfr_t> uVals(nu);
		std::cout<<"	toto is ok -4"<<std::endl;
		std::vector<mpfr_t> t(nt);
		std::cout<<"	toto is ok -3"<<std::endl;
		mpfr_t coeff, resMul, result;
		std::cout<<"	toto is ok -2"<<std::endl;
		mpfr_init2(resMul, veryLargePrec);
		std::cout<<"	toto is ok -1"<<std::endl;
		mpfr_init2(result, veryLargePrec);
		std::cout<<"	toto is ok 0"<<std::endl;
		for (uint32_t i=0; i<nt; i++) {
			mpfr_init2(t[i], veryLargePrec);
			mpfr_init2(uVals[i], veryLargePrec);
		}

		std::cout<<"	toto is ok 1"<<std::endl;
		for (uint32_t i = 0; i<nu; i++){
		std::cout<<"	toto is ok 1.0"<<std::endl;
				uValues[i]=tc->getInputValue(join("U",i));
		std::cout<<"	toto is ok 1.1"<<std::endl;
		std::cout<<"uValues["<<i<<"]="<<uValues[i]<<std::endl;
		std::cout<<"lsbIn["<<i<<"]="<<lsbIn[i]<<std::endl;
			mpfr_set_z( uVals[i],uValues[i].get_mpz_t(),GMP_RNDD );
			mpfr_mul_2si( uVals[i], uVals[i], lsbIn[i], GMP_RNDD );
		std::cout<<"uVals["<<i<<"]=";
			mpfr_printf( "%.65RNf\n", uVals[i]);
		std::cout<<"	toto is ok 1.2"<<std::endl;
				
		}
		std::cout<<"	toto is ok 2"<<std::endl;

		for (uint32_t l=0; l<nt+nx+ny; l++){
			for (uint32_t c=0; c<nt+nx+nu; c++){

		std::cout<<"	toto is ok 3"<<std::endl;
				sollya_obj_t stringCoeff;
				//get coeff as sollya object
				stringCoeff = sollya_lib_parse_string(coeffs[l][c].c_str());
				// If conversion did not succeed (i.e. parse error)
				if(stringCoeff == 0)	{
						std::ostringstream error;
						error << srcFileName << ":"<<"emulate(): Unable to parse std::string " << coeffs[l][c] << " as a numeric constant" << std::endl;
						throw error.str();
				}
		std::cout<<"	toto is ok 4"<<std::endl;
				mpfr_init2(coeff, veryLargePrec);
				//get coeff as mpfr
				sollya_lib_get_constant(coeff, stringCoeff);
				sollya_lib_clear_obj(stringCoeff);
		std::cout<<"	toto is ok 5"<<std::endl;
				if (coeff != 0){
					if ( !(l==c && l<nt) ) { //ones on the J diagonal are implicit
						if (c<nt){
		std::cout<<"	toto is ok 5.1"<<std::endl;

							mpfr_mul(resMul, t[c], coeff, GMP_RNDN);
		std::cout<<"	toto is ok 5.2"<<std::endl;
							mpfr_add(result, result, resMul, GMP_RNDN);
						}
						else if (c<nt+nx){
		std::cout<<"	toto is ok 6"<<std::endl;
							mpfr_mul(resMul, xPrec[c], coeff, GMP_RNDN);
		std::cout<<"	toto is ok 6.1"<<std::endl;
							mpfr_add(result, result, resMul, GMP_RNDN);
		std::cout<<"	toto is ok 7"<<std::endl;

						}
						else {
		std::cout<<"	toto is ok 7.1"<<std::endl;

							mpfr_mul(resMul, uVals[c-nt-nx], coeff, GMP_RNDN);
		std::cout<<"	toto is ok 7.2"<<std::endl;
							mpfr_add(result, result, resMul, GMP_RNDN);
		std::cout<<"	toto is ok 8"<<std::endl;
						}
					}
				}
			}
			if (l<nt){
				mpfr_set(t[l],result, GMP_RNDN);
		std::cout<<"	toto is ok 9"<<std::endl;
			}
			else if (l<nt+nx) {
				mpfr_set(xCurrent[l], result, GMP_RNDN);
		std::cout<<"	toto is ok 10"<<std::endl;
			}
			if (l>=nx+nt){
		std::cout<<"	toto is ok 11.0"<<std::endl;
				mpfr_t rd, ru;
				mpfr_init2 (rd, 1+msbIn[l]-lsbOut[l]);
				mpfr_init2 (ru, 1+msbIn[l]-lsbOut[l]);		
				mpz_class rdz, ruz;
		std::cout<<"	toto is ok 11"<<std::endl;

				mpfr_get_z (rdz.get_mpz_t(), result, GMP_RNDD); 					// there can be a real rounding here
				rdz=signedToBitVector(rdz, 1-lsbIn[l]-lsbOut[l]);
//TODO: here we should have computed msbIns
//TODO: replace lsbIn by MsbIn
		std::cout<<"	toto is ok 12"<<std::endl;

				mpfr_get_z (ruz.get_mpz_t(), result, GMP_RNDU); 					// there can be a real rounding here	
				ruz=signedToBitVector(ruz, 1-lsbIn[l]-lsbOut[l]);
				
				mpfr_clears ( result, rd, ru, NULL);

		std::cout<<"	toto is ok 13"<<std::endl;
				tc->addExpectedOutput(join("Y",l-nt-nx), rdz);
				tc->addExpectedOutput(join("Y",l-nt-nx), ruz);
		std::cout<<"	toto is ok 14"<<std::endl;
			}

		}
		for ( uint32_t i=0; i<nx; i++ ) {
		std::cout<<"	toto is ok 15"<<std::endl;
			mpfr_set( xPrec[i], xCurrent[i], GMP_RNDN );
		}



//		mpfr_set_d(s, 0.0, GMP_RNDN); // initialize s to 0
//
//		for (int i=0; i< n; i++)
//		{
//			u_ = xHistory[(currentIndexB+n-i)%n];		// get the input bit std::vector as an integer		
//			u_ = bitVectorToSigned(sx, 1+p); 						// convert it to a signed mpz_class		
//			mpfr_set_z (x, sx.get_mpz_t(), GMP_RNDD); 				// convert this integer to an MPFR; this rounding is exact
//			mpfr_div_2si (x, x, p, GMP_RNDD); 						// multiply this integer by 2^-p to obtain a fixed-point value; this rounding is again exact
//
//			mpfr_mul(t, x, mpcoeffb[i], GMP_RNDN); 					// Here rounding possible, but precision used is ridiculously high so it won't matter
//
//			if(coeffsignb[i]==1)
//				mpfr_neg(t, t, GMP_RNDN); 
//
//			mpfr_add(s, s, t, GMP_RNDN); 							// same comment as above
//			
//		}
//
//		for (int i=0; i<m; i++)
//		{
//
//			mpfr_mul(u, yHistory[(currentIndexA+m-i-1)%m], mpcoeffa[i], GMP_RNDN); 					// Here rounding possible, but precision used is ridiculously high so it won't matter
//
//			if(coeffsigna[i]==1)
//				mpfr_neg(u, u, GMP_RNDN); 
//
//			mpfr_add(s, s, u, GMP_RNDN); 							// same comment as above
//
//		}
//
//		mpfr_set(yHistory[currentIndexA], s, GMP_RNDN);
//
//
//		// now we should have in s the (exact in most cases) sum
//		// round it up and down
//
//		// make s an integer -- no rounding here 
//		mpfr_mul_2si (s, s, p, GMP_RNDN);
//
//		
//		// We are waiting until the first meaningful value comes out of the IIR
//
//		mpz_class rdz, ruz;
//
//		mpfr_get_z (rdz.get_mpz_t(), s, GMP_RNDD); 					// there can be a real rounding here
//		rdz=signedToBitVector(rdz, wO);
//		tc->addExpectedOutput ("R", rdz);
//
//		mpfr_get_z (ruz.get_mpz_t(), s, GMP_RNDU); 					// there can be a real rounding here	
//		ruz=signedToBitVector(ruz, wO);
//		tc->addExpectedOutput ("R", ruz);
//
//		
//		mpfr_clears (x, t, u, s, NULL);
//
//		currentIndexB = (currentIndexB +1)%n; // We use a circular buffer to store the inputs
//		currentIndexA = (currentIndexA +1)%m;
//
//
//#else
//		static int idx = 0;
//		static bool full = false; 							// set to true when the fir start to output valid data (after n input) 
//		static TestCase * listTC [10000]; // should be enough for everybody
//
//
//		listTC[idx] = tc;
//
//		if(n == 1)					// if the fir part has only one tap we don't wait to get the output
//			full = true; 
//
//		// We are waiting until the first meaningful value comes out of the FIR
//		if (full) {
//			mpfr_t x, t, s, rd, ru;
//			mpfr_init2 (x, 1+p);
//			mpfr_init2 (t, 10*(1+p));
//			mpfr_init2 (s, 10*(1+p));
//			mpfr_init2 (rd, 1+p);
//			mpfr_init2 (ru, 1+p);		
//
//			mpfr_set_d(s, 0.0, GMP_RNDN); // initialize s to 0
//
//
//			int k = idx; // We start to sum from the last input
//
//			for (int i=0; i< n; i++)
//			{
//
//				mpz_class sx = listTC[k]->getInputValue("X"); 		// get the input bit std::vector as an integer
//				sx = bitVectorToSigned(sx, 1+p); 						// convert it to a signed mpz_class
//				mpfr_set_z (x, sx.get_mpz_t(), GMP_RNDD); 				// convert this integer to an MPFR; this rounding is exact
//				mpfr_div_2si (x, x, p, GMP_RNDD); 						// multiply this integer by 2^-p to obtain a fixed-point value; this rounding is again exact
//
//				mpfr_mul(t, x, mpcoeff[i], GMP_RNDN); 					// Here rounding possible, but precision used is ridiculously high so it won't matter
//
//				if(coeffsign[i]==1)
//					mpfr_neg(t, t, GMP_RNDN); 
//
//				mpfr_add(s, s, t, GMP_RNDN); 							// same comment as above
//			
//				k = (k+1)%n;	
//			}
//
//			k = (k-1+n)%n; //to get the corresponding testCase to the outputed value
//
//			// now we should have in s the (exact in most cases) sum
//			// round it up and down
//
//			// make s an integer -- no rounding here
//			mpfr_mul_2si (s, s, p, GMP_RNDN);
//
//			mpz_class rdz, ruz;
//
//			mpfr_get_z (rdz.get_mpz_t(), s, GMP_RNDD); 					// there can be a real rounding here
//			rdz=signedToBitVector(rdz, wO);
//			listTC[k]->addExpectedOutput ("R", rdz);
//			// tc->addExpectedOutput ("R", rdz);
//
//			mpfr_get_z (ruz.get_mpz_t(), s, GMP_RNDU); 					// there can be a real rounding here	
//			ruz=signedToBitVector(ruz, wO);
//			listTC[k]->addExpectedOutput ("R", ruz);
//
//			mpfr_clears (x, t, s, rd, ru, NULL);
//		}
//		
//		idx = (idx-1+n)%n; // We use a circular buffer to store the inputs
//
//		if (idx ==  1) {
//			full = true;
//		}
//
//#endif
	};

	int FixSIF::readMatrices( std::vector< std::vector <std::vector<std::string> >**> &Z, std::ifstream &openedFile, int &lc){

		std::string line;
		for ( int i=0; i<9; i++ ){

			int gotLine=getFullLine(openedFile, line, lc);
			if(gotLine)
				return 1;

			REPORT(0,"Got header at line "<<lc<<": \""<<line<<"\"")
			if (Z.size()!=9)
				return 2;

			switch (line[0]){
				case 'J':
					lc = readMatrix(line, "J", *Z[0], openedFile, lc);
					break;
				case 'K':
					lc = readMatrix(line, "K", *Z[1], openedFile, lc);
					break;
				case 'L':
					lc = readMatrix(line, "L", *Z[2], openedFile, lc);
					break;
				case 'M':
					lc = readMatrix(line, "M", *Z[3], openedFile, lc);
					break;
				case 'N':
					lc = readMatrix(line, "N", *Z[4], openedFile, lc);
					break;
				case 'P':
					lc = readMatrix(line, "P", *Z[5], openedFile, lc);
					break;
				case 'Q':
					lc = readMatrix(line, "Q", *Z[6], openedFile, lc);
					break;
				case 'R':
					lc = readMatrix(line, "R", *Z[7], openedFile, lc);
					break;
				case 'S':
					lc = readMatrix(line, "S", *Z[8], openedFile, lc);
					break;
				default:
					return 1;
					break;
			}
					REPORT(0,"Popped back to readMatrices");
		}
		return 0;



	}

	int FixSIF::readMatrix(std::string &header, std::string JKLMNPQRS, std::vector < std::vector <std::string> > * &toFill, std::ifstream &openedFile, int lc){

			std::string line;

			int nc=0;
			int nl=0;
			REPORT(0, "Reading matrix given the header: \""<<header<<"\". Matrix descriptor: "<<JKLMNPQRS);
			if (header.find_first_of(JKLMNPQRS) != std::string::npos){
				size_t bNum = header.find_first_of(" ");
				size_t aNum = header.find_first_of(" ",bNum+1);
				if (bNum != std::string::npos){
					nl = stoi(header.substr(bNum+1,aNum-1));
					size_t bNum = header.find_first_of(" ",aNum);
					if (bNum != std::string::npos){
						nc = stoi(header.substr(aNum+1,bNum-1));
						toFill = new std::vector< std::vector<std::string> >(nl);
					}
					else {
							std::stringstream error;
							error<<"line "<<lc<<": space missing";
							THROWERROR(error.str());
					}
					
				}
				else {
						std::stringstream error;
						error<<"line "<<lc<<": no space found";
						THROWERROR(error.str());
				}
				REPORT(0, nc<<" columns and "<<nl<<" lines found in matrix from header: \""<<header<<"\"");
			}
			else{
				std::stringstream error;
				error<<"line "<<lc<<": Matrix descriptor "<< JKLMNPQRS <<" not found at its expected place";
				THROWERROR(error.str());
			}


			if (JKLMNPQRS == "J"){
				nt=nl;
			}
			else if (JKLMNPQRS == "M"){
				nx=nc;
			}
			else if (JKLMNPQRS == "N"){
				nu=nc;
			}
			else if (JKLMNPQRS == "L"){
				ny=nl;
			}

			//if empty matrix (column or line number=0)
			if (!nc ||!nl)
				return lc;

			for (int i=0; i<nl; i++){
				getFullLine(openedFile, line, lc);
				REPORT(0,"Got line "<<lc<<": \""<<line<<"\"")

				size_t bNum = 0;
				size_t aNum = 0;
				for (int j = 0; j<nc; j++) {
					aNum = line.find_first_of("	", bNum);
					if (aNum != std::string::npos){
						 (*toFill)[i].push_back(line.substr(bNum,aNum-bNum));
						 bNum=aNum+1;
					}

					else if (j==nc-1) {
						 (*toFill)[i].push_back(line.substr(bNum,line.size()-bNum));
					}
					else {
						std::stringstream error;
						error<<"line "<<lc<<": no tab found";
						THROWERROR(error.str());
					}
				}

			}
			
			return lc;

	}



	int FixSIF::parseFile(){

		int lineCounter=0;
		std::ifstream input;
		input.open( file, std::ifstream::in );
		std::vector < std::vector<std::string> > *J;
		std::vector < std::vector<std::string> > *K;
		std::vector < std::vector<std::string> > *L;
		std::vector < std::vector<std::string> > *M;
		std::vector < std::vector<std::string> > *N;
		std::vector < std::vector<std::string> > *P;
		std::vector < std::vector<std::string> > *Q;
		std::vector < std::vector<std::string> > *R;
		std::vector < std::vector<std::string> > *S;

		std::vector < std::vector< std::vector<std::string> >** > Z;
		Z.push_back( &J );
		Z.push_back( &K );
		Z.push_back( &L );
		Z.push_back( &M );
		Z.push_back( &N );
		Z.push_back( &P );
		Z.push_back( &Q );
		Z.push_back( &R );
		Z.push_back( &S );

		if ( input.is_open() ){

	//		lineCounter = readMatrix( "J", J, input );
	//		lineCounter = readMatrix( "K", K, input, lineCounter );
	//		lineCounter = readMatrix( "L", L, input, lineCounter );
	//		lineCounter = readMatrix( "M", M, input, lineCounter );
	//		lineCounter = readMatrix( "P", P, input, lineCounter );
	//		lineCounter = readMatrix( "R", R, input, lineCounter );
	//		lineCounter = readMatrix( "N", N, input, lineCounter );
	//		lineCounter = readMatrix( "Q", Q, input, lineCounter );
	//		lineCounter = readMatrix( "S", S, input, lineCounter );
			REPORT(0, "Parsing matrices");
			readMatrices(Z, input, lineCounter);
			input.close();
		}
		else {
			std::stringstream error;
			error << "Cannot open file "<<file<<std::endl;
			THROWERROR(error.str());
		}

		coeffs = std::vector< std::vector<std::string> > (nt+nx+ny); //FIXME: check if nt+nx+ny is effectiveley relevant, regarding nt+nx+nu.

		uint32_t i=0, j;
		for (; i<nt; i++){
			for (j=0; j<nt; j++){
				coeffs[i].push_back((*J)[i][j]);
			}
			for (; j<nt+nx; j++){
				coeffs[i].push_back((*M)[i][j-nt]);
			}
			for (; j<nt+nx+nu; j++){
				coeffs[i].push_back((*N)[i][j-nt-nx]);
			}
		}
		for (; i<nt+nx; i++){
			for (j=0; j<nt; j++){
				coeffs[i].push_back((*K)[i-nt][j]);
			}
			for (; j<nt+nx; j++){
				coeffs[i].push_back((*P)[i-nt][j-nt]);
			}
			for (; j<nt+nx+nu; j++){
				coeffs[i].push_back((*Q)[i-nt][j-nt-nx]);
			}
		}
		for (; i<nt+nx+nu; i++){
			for (j=0; j<nt; j++){
				coeffs[i].push_back((*L)[i-nt-nx][j]);
			}
			for (; j<nt+nx; j++){
				coeffs[i].push_back((*R)[i-nt-nx][j-nt]);
			}
			for (; j<nt+nx+nu; j++){
				coeffs[i].push_back((*S)[i-nt-nx][j-nt-nx]);
			}
		}
		displayMatrix(*J,"J");
		displayMatrix(*K,"K");
		displayMatrix(*L,"L");
		displayMatrix(*M,"M");
		displayMatrix(*P,"P");
		displayMatrix(*R,"R");
		displayMatrix(*N,"N");
		displayMatrix(*Q,"Q");
		displayMatrix(*S,"S");
		displayMatrix(coeffs,"coeffs");

		return 0;
	}

	/**
		bMToDoubleM: fills the matrix doubleM from the coefficients stored in bM.
		Note: this assumes that doubleM is instanciated
		input: 
			-bM: boost matrix to convert in an array-matrix j
			-doubleM: the array to fill
		Ouptut: 0 if success.
	*/
	int FixSIF::bMToDoubleM( ublas::matrix<double> const &bM, double * &doubleM ) {
		std::cout<<"entering bMToDoubleM"<<std::endl;	
		for (size_t i=0; i<bM.size1(); i++) {
			for (size_t j=0; j<bM.size2(); j++) {
				doubleM [i*bM.size1()+j] = bM(i,j);
			}
		}
		std::cout<<"existing bMToDoubleM"<<std::endl;	
		return 0;
	}


	/**
		vvToBoostMatrix: fills the matrix bM from the std::string coefficients stored in sM.
		Note: this assumes that bM is instanciated
		Input:
			-sM: std::vector of std::vector of std::strings (it won't be modified)
			-bM: boost matrix to be filled
		Output:
			-0 if success
	  */
	int FixSIF::vvToBoostMatrix( std::vector< std::vector <std::string> > const &sM, ublas::matrix<double> &bM ) {

		std::cout<<"entering vvToBoostMatrix"<<std::endl;

		for ( unsigned int i=0; i<sM.size(); i++) {
			for ( unsigned int j=0; j<sM[i].size(); j++ ) {
				std::cout<<"iteration: i="<<i<<", j="<<j<<std::endl;
				bM(i,j)=atof( sM[i][j].c_str() );
			}
		}
		std::cout<<"exiting vvToBoostMatrix"<<std::endl;
		return 0;
	}

	int FixSIF::computeABCD(ublas::matrix<double> const &bJ, ublas::matrix<double> const &bK, ublas::matrix<double> const &bL, ublas::matrix<double> const &bM, ublas::matrix<double> const &bN, ublas::matrix<double> const &bP, ublas::matrix<double> const &bQ, ublas::matrix<double> const &bR, ublas::matrix<double> const &bS, ublas::matrix<double> &bA, ublas::matrix<double> &bB, ublas::matrix<double> &bC, ublas::matrix<double> &bD){

		std::cout<<"entering computeABCD"<<std::endl;
		
		ublas::matrix<double>invbJ(bJ.size1(),bJ.size2());
		invertMatrix(bJ,invbJ);

		ublas::matrix<double>A1(bK.size1(), bJ.size2());
		ublas::matrix<double>A2(bP.size1(), bP.size2());

		ublas::matrix<double>B1(bK.size1(), bJ.size2());
		ublas::matrix<double>B2(bQ.size1(), bQ.size2());

		ublas::matrix<double>C1(bK.size1(), bJ.size1());
		ublas::matrix<double>C2(bR.size1(), bR.size2());

		ublas::matrix<double>D1(bK.size1(), bJ.size1());
		ublas::matrix<double>D2(bS.size1(), bS.size2());

		std::cout<<"computeABCD false for now"<< std::endl;

		if (!matMul(A1, bK, invbJ))
			return false;
		if (!matMul(A2, A1, bM   ))
			return false;
		if (!matAdd(bA, A2, bP   ))
			return false;
                                  
		if (!matMul(B1, bK, invbJ))
			return false;
		if (!matMul(B2, B1, bN   ))
			return false;
		if (!matAdd(bB, B2, bQ   ))
			return false;
                                  
		if (!matMul(C1, bL, invbJ))
			return false;
		if (!matMul(C2, C1, bM   ))
			return false;
		if (!matAdd(bC, C2, bR   ))
			return false;
                                  
		if (!matMul(D1, bL, invbJ))
			return false;
		if (!matMul(D2, D1, bN   ))
			return false;
		if (!matAdd(bD, D2, bS   ))
			return false;
		std::cout<<"exiting computeABCD ok"<<std::endl;
		return true;
	}
	
	int FixSIF::computeMSBsLSBs(){
		//A, B, C, D are double pointers
		//J ,K, L, M, N, P, Q, R, and S are boost double matrices
		//prefix b in bX denotes a boost matrix
		//First: computing MSBs

		std::cout<<"entering computeMSBsLSBs"<<std::endl;

		//Init ABCD
		double *A = new double[nx*nt] ;
		double *B = new double[nx*nu] ;
		double *C = new double[ny*nx] ;
		double *D = new double[ny*nu] ;

		//Declare Z
		ublas::matrix<double> Z(nt+nx+ny,nt+nx+nu);
		//Init Z from coeffs
		vvToBoostMatrix(coeffs, Z);
		

		//Declare JKLMNPQRS
		ublas::matrix<double> J(nt,nt);
		ublas::matrix<double> K(nx,nt);
		ublas::matrix<double> L(ny,nt);
		ublas::matrix<double> M(nt,nx);
		ublas::matrix<double> N(nt,nu);
		ublas::matrix<double> P(nx,nx);
		ublas::matrix<double> Q(nx,nu);
		ublas::matrix<double> R(ny,nx);
		ublas::matrix<double> S(ny,nu);

		//Init JKLMNPQRS from Z
		unsigned int i;
		unsigned int j;
		for ( i = 0; i<nt; i++ ){
			for ( j = 0; j<nt; j++ ){
				J(i,j) = Z(i,j);
			}

			for ( j = 0; j<nx; j++ ){
				M(i,j) = Z(i,j+nt);
			}

			for ( j = 0; j<nu; j++ ){
				N(i,j) = Z(i,j+nt+nx);
			}
		}
		
		for ( i = 0 ; i<nx; i++ ){
			for ( j = 0; j<nt; j++ ){
				K(i,j) = Z(i+nt,j);
			}

			for ( j = 0; j<nx; j++ ){
				P(i,j) = Z(i+nt,j+nt);
			}

			for ( j = 0; j<nu; j++ ){
				Q(i,j) = Z(i+nt,j+nt+nx);
			}
		}
		
		
		for ( i = 0 ; i<ny; i++ ){
			for ( j = 0; j<nt; j++ ){
				L(i,j) = Z(i+nt+nx,j);
			}

			for ( j = 0; j<nx; j++ ){
				R(i,j) = Z(i+nt+nx,j+nt);
			}

			for ( j = 0; j<nu; j++ ){
				S(i,j) = Z(i+nt+nx,j+nt+nx);
			}
		}

		//Declare boost ABCD representation
		ublas::matrix<double> bA(nx,nx);
		ublas::matrix<double> bB(nx,nu);
		ublas::matrix<double> bC(ny,nx);
		ublas::matrix<double> bD(ny,nu);

		//Compute ABCD from JKLMNPQRS into boost ABCD
		computeABCD(J, K, L, M, N, P, Q, R, S, bA, bB, bC, bD);

		//Transfert into double* format
		bMToDoubleM( bA, A);
		bMToDoubleM( bB, B);
		bMToDoubleM( bC, C);
		bMToDoubleM( bD, D);

		//Declare wcpg of the filter
		double *wcpgF=new double [(nx+ny)*(nx+nu)];

		//Compute WCPG
		//WCPG_ABCD( wcpgF, A, B, C, D, uint64_t(nx), uint64_t(nu), uint64_t(ny) );
		WCPG_ABCD( wcpgF, A, B, C, D, nx, nu, ny );

		int maxMSB=0;
		for ( unsigned int i=0; i<ny; i++) {

			if (ceil(log((pow(2,msbIn[j]+1)-1)*wcpgF[i*nu+j])/log(2)) > maxMSB)
				maxMSB=ceil(log((pow(2,msbIn[j]+1)-1)*wcpgF[i*nu+j])/log(2));
		}

		//Init all msbsOut to 0
		//FIXME: find a better upper bound for MSBs
		for ( unsigned int i=0; i<nt+nx; i++) {
			msbInter.push_back(maxMSB);
		}

		//msbsOut=wcpgF*max(U)
		//max(U)=2^(msbIn+1)-1
		for ( unsigned int i=nt+nx; i<ny; i++ ) {
			for ( unsigned int j=0; j<nu; j++) {
				msbInter.push_back(ceil(log((pow(2,msbIn[j]+1)-1)*wcpgF[i*nu+j])/log(2)));
			}
		}



		//Declare pseudo identities
		ublas::matrix<double> Mt(nt,nt+nx+ny);
		ublas::matrix<double> Mx(nx,nt+nx+ny);
		ublas::matrix<double> My(ny,nt+nx+ny);

		//Init Mt
		for ( unsigned int i=0; i<nt; i++ ) {
			for ( unsigned int j=0; j<nt+nx+ny; j++ ) {
				if ( (j<nt) && (i=j) ) {
					Mt(i,j)=1;
				}
				else {
					Mt(i,j)=0;
				}
			}
		}

		//Init Mx
		for ( unsigned int i=0; i<nx; i++ ) {
			for ( unsigned int j=0; j<nt+nx+ny; j++ ) {
				if ( (j>=nt) && (j<nt+nx) && (i=j) ) {
					Mx(i,j)=1;
				}
				else {
					Mx(i,j)=0;
				}
			}
		}

		//Init My
		for ( unsigned int i=0; i<ny; i++ ) {
			for ( unsigned int j=0; j<nt+nx+ny; j++ ) {
				if ( (j>=nt+nx) && (i=j) ) {
					My(i,j)=1;
				}
				else {
					My(i,j)=0;
				}
			}
		}

		//Declare boost ABCD representation for the error filter
		ublas::matrix<double> bAe(nx,nx);
		ublas::matrix<double> bBe(nx,nu);
		ublas::matrix<double> bCe(ny,nx);
		ublas::matrix<double> bDe(ny,nu);

		//Compute ABCD error filter
		computeABCD(J, K, L, M, Mt, P, Mx, R, My, bAe, bBe, bCe, bDe);

		//Init ABCD for error filter
		double *Ae = new double[nx*nt] ;
		double *Be = new double[nx*nu] ;
		double *Ce = new double[ny*nx] ;
		double *De = new double[ny*nu] ;

		//Transfert into double* format
		bMToDoubleM( bAe, Ae );
		bMToDoubleM( bBe, Be );
		bMToDoubleM( bCe, Ce );
		bMToDoubleM( bDe, De );

		//Declare wcpgE of the error filter
		double *wcpgE=new double [(nt+nx+ny)*(ny)];

		//Compute WCPG error
		WCPG_ABCD( wcpgE, A, B, C, D, nx, nt+nx+ny, ny );



		soplex::SoPlex mysoplex;

	   /* set the objective sense */
	   mysoplex.setIntParam(soplex::SoPlex::OBJSENSE, soplex::SoPlex::OBJSENSE_MAXIMIZE);

	   /* we first add variables */
	   soplex::DSVector dummycol(0);
	   //set weight to 1 in the objective function and upper bound to 1
	   //FIXME: compute a better lower bound for better performance
	   for ( unsigned int i =0; i<nt+nx+ny; i++){
	   mysoplex.addColReal(soplex::LPCol(1.0, dummycol, 1.0, -(soplex::infinity)));
	   }

	   /* then constraints one by one */
	   for (unsigned int i=0; i<ny; i++){
		   soplex::DSVector row(nt+nx+ny);
		   for (unsigned int j=0; j<nt+nx+ny; j++) {
			   row.add(j, wcpgE[j+i*(nt+nx+ny)]);
		   }
		   mysoplex.addRowReal( soplex::LPRow( -(soplex::infinity), row, pow(2,lsbOut[i]-1) ) );
	   }

	   /* NOTE: alternatively, we could have added the matrix nonzeros in dummycol already; nonexisting rows are then
		* automatically created. */

	   /* write LP in .lp format */
	   mysoplex.writeFileReal("dump.lp", NULL, NULL, NULL);

	   /* solve LP */
	   soplex::SPxSolver::Status stat;
	   soplex::DVector prim(nt+nx+ny);
	   soplex::DVector dual(nt+nx+ny);
	   stat = mysoplex.solve();

	   /* get solution */
	   if( stat == soplex::SPxSolver::OPTIMAL )
	   {
		  mysoplex.getPrimalReal(prim);
		  mysoplex.getDualReal(dual);
		  std::cout << "LP solved to optimality.\n";
		  std::cout << "Objective value is " << mysoplex.objValueReal() << ".\n";
		  std::cout << "Primal solution is [" << prim[0] << ", " << prim[1] << "].\n";
		  std::cout << "Dual solution is [" << dual[0] << "].\n";
		  //std::cout << "LP solved to optimality.\n";
		  //std::cout << "Objective value is " << mysoplex.objValueReal() << ".\n";
		  //std::cout << "Primal solution is [" << prim[0] << ", " << prim[1] << "].\n";
		  //std::cout << "Dual solution is [" << dual[0] << "].\n";
	   }
	   for ( unsigned int i=0; i<nt+nx+ny; i++ ) {
		   lsbInter[i]=ceil(log(prim[i])/log(2));
	   }
	   //real toto=0.7;
	   //toto+=41.3;



		



		std::cout<<"exiting computeMSBsLSBs"<<std::endl;
		return 0;
	}

	int FixSIF::readPrecision( std::vector<int> &msbsIn, std::vector<int> &lsbsIn, std::vector<int> &msbsOut, std::vector<int> &lsbsOut, bool inFile ){
		//TODO:return wcpg*error min
		REPORT(0,"Entering readPrecision");

		if ( inFile ){
			REPORT(0,"Filling with file information");
			std::vector < std::vector<std::string> > *msbslsbsIn, *msbslsbsOut;
			std::ifstream precisions;
			precisions.open( "precisions.txt", std::ifstream::in );
			std::stringstream head;
			head<<"P "<<nt+nx+ny<<" "<<2;
			std::string header=head.str();
			readMatrix( header, "IN", msbslsbsIn, precisions );
			for ( uint32_t i=0; i<msbslsbsIn->size(); i++ ){
				msbsIn.push_back(stoi((*msbslsbsIn)[i][0]));
				lsbsIn.push_back(stoi((*msbslsbsIn)[i][1]));
			}
			readMatrix( header, "OUT", msbslsbsOut, precisions );
			for ( uint32_t i=0; i<msbslsbsOut->size(); i++ ){
				msbsOut.push_back(stoi((*msbslsbsOut)[i][0]));
				lsbsOut.push_back(stoi((*msbslsbsOut)[i][1]));
			}
			REPORT(0,"exiting by file reading path");
			return 0;
			precisions.close();
		}
		else {
			REPORT(0,"Filling with default precision");
			for ( uint32_t i=0; i<nu+nt+nx+ny; i++ ) {
				msbsIn.push_back(DEFAULT_MSB);
				lsbsIn.push_back(DEFAULT_LSB);
				msbsOut.push_back(DEFAULT_MSB);
				lsbsOut.push_back(DEFAULT_LSB);
			}
			REPORT(0,"exiting by default precision path");
			return 0;
		}
		REPORT(0,"exiting readPrecision with error");
		return 1;
	}
	
}
	
