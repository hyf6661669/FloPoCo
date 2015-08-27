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
#define WCPG_F //disable WCPG until it works

#ifndef WCPG_F
extern "C" {
#include "/usr/local/include/wcpg_API.h"
}
#endif
//#include "/usr/local/include/clapack_headers_config.h"	

#define BUF_SIZE 256 //default buffer size for parsing the specification file
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

		std::cout<<"entering matMul: X.size1()="
			<<X.size1()
			<<", X.size2()="
			<<X.size2()
			<<", Y.size1()="
			<<Y.size1()
			<<", Y.size2()="
			<<Y.size2()
		<<std::endl;

	if(X.size2()!=Y.size1()) {
		std::cout<<"matMul: failed to apply matrix product"<<std::endl;
		return false;
	}

	for (size_t i=0; i<X.size1(); i++) {
		for (size_t j=0; j<Y.size2(); j++) {
			res(i,j)=0.0;
		}
	}
	
	for (size_t i=0; i<X.size1(); i++) {
		for (size_t j=0; j<Y.size2(); j++) {
			for (size_t k=0; k<X.size2(); k++) {
				res(i,j) += X(i,k) * Y(k,j);
			}
		}
	}
	std::cout<<"exiting matMul successfully"<<std::endl;
	return true;
}

bool matAdd (ublas::matrix<double> &res, ublas::matrix<double> const &X, ublas::matrix<double> const &Y) { 

		std::cout<<"entering matAdd: X.size1()="
			<<X.size1()
			<<", X.size2()="
			<<X.size2()
			<<", Y.size1()="
			<<Y.size1()
			<<", Y.size2()="
			<<Y.size2()
		<<std::endl;


	if ( ( X.size1() != Y.size1() ) || ( X.size2() != X.size2() ) ) {
		std::cout<<"matAdd: failed to apply matrix addition"<<std::endl;
		return false;
	}

	for (size_t i=0; i<X.size1(); i++) {
		for (size_t j=0; j<X.size2(); j++) {
			res(i,j) = X(i,j) + Y(i,j);
			std::cout<<"matAdd: adding "<<X(i,j)<<" with " <<Y(i,j)<<". Result gives: "<<X(i,j)+Y(i,j)<<std::endl;
		}
	}
	std::cout<<"exiting matAdd successfully"<<std::endl;
	return true;
}

//double maxi( double *wcpg, int l, int c, int cs) {
//	
//}

	int getFullLine( std::ifstream &f, std::string &s, int &lc ){

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
					return 0;
				}

				if ((!isdigit(s[i])) && (i==s.size()-1) && (!f.eof()) ) {
					getline(f, s);
					lc++;
					i=-1;
				}
			}
			return 1;
	}

//TODO: build a policy to detect suitable use of shiftRegisters.
namespace flopoco {
	const int veryLargePrec = 6400;  /*6400 bits should be enough for anybody */

	int FixSIF::displayMatrix( std::vector< std::vector<std::string> > &m, std::string name ){
		REPORT(3, "Displaying Matrix "<<name<<" "<<m.size()<<" "<<m[0].size());
		for (uint32_t i=0; i<m.size(); i++){
			std::stringstream line;
			for (uint32_t  j=0; j<m[i].size(); j++){
				line<<m[i][j];
				if (j!=m[i].size()-1){
					line<<"	";
				}
			}
			REPORT(3,line.str());
		}
		return 0;
	}


	FixSIF::FixSIF(Target* target, int lsb_, std::string file_, std::map<std::string, double> inputDelays) : 
		Operator(target), file(file_)
	{
		//target->setNotPipelined();

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

		std::ostringstream precDisp;
		precDisp<<"msbInter = [";
		for ( unsigned int i=0; i<msbInter.size(); i++) {
			precDisp<<" "<<msbInter[i];
		}
		precDisp<<" ]";
		REPORT(3,precDisp.str());
		precDisp.str(std::string());

		precDisp<<"lsbInter = [";
		for ( unsigned int i=0; i<lsbInter.size(); i++) {
			precDisp<<" "<<lsbInter[i];
		}
		precDisp<<" ]";
		REPORT(3,precDisp.str());


		useNumericStd();
		//useNumericStd_Unsigned();
		//useNumericStd_Signed();
		//if(p<1) {
		//	THROWERROR("Can't build an architecture for this value of LSB")
		//}

		// initialize stuff for emulate
		//xPrec= (mpfr_t *) malloc ( nx * sizeof(mpfr_t*));
		xPrec= new mpfr_t[nx];
		//xCurrent= (mpfr_t *) malloc ( nx * sizeof(mpfr_t*));
				std::cout<<"xCurrent="<<xCurrent<<std::endl;
		xCurrent= new mpfr_t[nx];
				std::cout<<"xCurrent="<<xCurrent<<std::endl;
		for (uint32_t i=0; i<nx; i++){
			mpfr_init2( xPrec[i], veryLargePrec );
			mpfr_init2( xCurrent[i], veryLargePrec );
			mpfr_set_ui( xPrec[i], 0, GMP_RNDN );
			mpfr_set_ui( xCurrent[i], 0, GMP_RNDN );
		}

		REPORT(3,"Emulate init done");

		//declare intermediate T(k+1)
		for (uint32_t i=0; i<nt; i++) {
			declare(join("T",i), msbInter[i]-lsbInter[i]+1);
		}

		//nu inputs numbering from 0 to nu-1 FIXME: check the relevance of this convention dealing with SIF
		for (uint32_t i=0; i<nu; i++) {
			std::ostringstream input;
			input << "U" <<i;
			//addFixInput(std::string(input.str()), true, 1+msbIn[i], lsbIn[i]);
			addInput(std::string(input.str()), 1+msbIn[i]-lsbIn[i]);
			declare(join("Ulv",i), msbInter[nt+i]-lsbInter[nt+i]+1);
			vhdl<< tab << join("Ulv",i) << " <= " << "std_logic_vector(U"<<i<<");"<<std::endl;
		}

		//ny outputs numbering from 0 to ny-1 FIXME: check the relevance of this convention dealing with SIFs
		for (uint32_t i=0; i<ny; i++) {
			std::ostringstream output;
			//output << "Y";// << i;
			std::string name=join("Y",i);
			addFixOutput(name, 1, msbOut[i], lsbOut[i]);
			declare(join("Ylv",i), 1+msbOut[i]-lsbOut[i]);
		}
	
		//declare intermediate X and X(k+1)
		for (uint32_t i=0; i<nx; i++) {
			//declare(join("Xplus",i));
			declare(join("X",i), msbInter[nt+i]-lsbInter[nt+i]+1);
			declare(join("Xplus",i), msbInter[nt+i]-lsbInter[nt+i]+1);
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
				REPORT(3,"non zeros coefficients number="<<nonZeros.size());
			for (uint32_t i = 0; i<coefDel[0].size(); i++){
				nonZeros[i] = coefDel[0][i];
				REPORT(3,"coeffs["<<i<<"]="<<coefDel[0][i]);
			}
			//output matrix std::vector
			//three steps (t, x and u) because all are not treated the same way
			//drop zeros and (implicit) ones in the diagonal but keep indices for consistency
			int ib = 0; //bias in std::vector indices induced by erase operations
			uint32_t i = 0;
			std::vector <int> lsbInC, msbInC;
			for ( ; i<nt; i++ ){
				if ( ( std::stof(nonZeros[i+ib].c_str()) == 0.0 ) || ( (std::stof(nonZeros[i+ib].c_str()) == 1.0)&&(n-coefDel.size()==i)) ){
					REPORT(3,"no coef found at "<<i+ib);
					REPORT(3,"numerical value of non zero coeff at "<<i+ib<<":"<<-std::stof(nonZeros[i+ib].c_str()));
					nonZeros.erase(nonZeros.begin()+i+ib);
					ib--;
				}
				else {
					coefIndst.push(i);
					lsbInC.push_back(lsbIn[i]);
					msbInC.push_back(msbIn[i]);
					REPORT(3,"coef found at "<<i+ib);
					REPORT(3,"numerical value of non zero coeff at "<<i+ib<<":"<<-std::stof(nonZeros[i+ib].c_str()));
					//hadling "-" for the use of -J in the calculations
					//this operation is done here because we need the positive version J for all calculations concerning precisions.
					if ( nonZeros[i+ib][0]!='-' ) {
						nonZeros[i+ib]="-"+nonZeros[i+ib];
					}
					else {
						nonZeros[i+ib]=nonZeros[i+ib].substr(1,std::string::npos);
					}
					REPORT(3,"pushing "<<i<<" in list of indices for T"); } } REPORT(3,"number of non zeros coefficients:"<<nonZeros.size()<<", bias="<<ib<<", iteration="<<i<<", n=nu+nx+nt="<<nu+nx+nt); //only drop zeros for x
			for ( ; i<nt+nx; i++ ){
				if ( std::stof(nonZeros[i+ib].c_str()) == 0.0 ){
					REPORT(3,"no coef found at "<<i+ib);
					REPORT(3,"numerical value of non zero coeff at "<<i+ib<<":"<<std::stof(nonZeros[i+ib].c_str()));
					nonZeros.erase(nonZeros.begin()+i+ib);
					ib--;
				}
				else{
					coefIndsx.push(i-nt);
					lsbInC.push_back(lsbIn[i]);
					msbInC.push_back(msbIn[i]);
					REPORT(3,"coef found at "<<i+ib);
					REPORT(3,"numerical value of non zero coeff at "<<i+ib<<":"<<std::stof(nonZeros[i+ib].c_str()));
					REPORT(3, "pushing "<<i-nt<<" in list of indices for x");
				}
			}

			REPORT(3, "cleaning u from implicit zeros, bias="<<ib<<", iteration="<<i<<", n=nu+nx+nt="<<nu+nx+nt);

			//and for u
			for ( ; i<nu+nx+nt; i++ ){
				if ( std::stof(nonZeros[i+ib].c_str()) == 0.0 ){
					REPORT(3,"no coef found at "<<i+ib);
					REPORT(3,"numerical value of non zero coeff at "<<i+ib<<":"<<std::stof(nonZeros[i+ib].c_str()));
					nonZeros.erase(nonZeros.begin()+i+ib);
					ib--;
				}
				else{
					coefIndsu.push(i-nt-nx);
					lsbInC.push_back(lsbIn[i]);
					msbInC.push_back(msbIn[i]);
					REPORT(3,"coef found at "<<i+ib);
					REPORT(3,"numerical value of non zero coeff at "<<i+ib<<":"<<std::stof(nonZeros[i+ib].c_str()));
					REPORT(3, "pushing "<<i-nt-nx<<" in list of indices for x");
				}
			}

			unsigned int indstS = coefIndst.size();
			unsigned int indsxS = coefIndsx.size();
			unsigned int indsuS = coefIndsu.size();

			std::stringstream indOut;
			indOut<< "coefIndst ("<< coefIndst.size() << ")";
			REPORT(3,indOut.str());
			indOut.str( std::string() );

			indOut<< "coefIndsx ("<< coefIndsx.size() << ")";
			REPORT(3,indOut.str());
			indOut.str( std::string() );


			indOut<< "coefIndsu ("<< coefIndsu.size() << ")";
			REPORT(3,indOut.str());


			//if only one coeff, just write a multiplier, if only ones, just write an adder
			//if ( nonZeros.size() > 1 ) 
					
					REPORT(3, "Calling FixSOPC with arguments:");
					REPORT(3, "	target="<<target_);
					for ( unsigned int deb=0; deb<msbIn.size(); deb++ ) {
						REPORT(3, "	msbIn["<<deb<<"]="<<msbIn[deb]);
						REPORT(3, "	lsbIn["<<deb<<"]="<<lsbIn[deb]);
					}
					REPORT(3, "	sizeof lsbIn="<<lsbIn.size());
					REPORT(3, "	msbOut="<<*msbInter.begin());
					REPORT(3, "	lsbOut="<<*lsbInter.begin());
					for ( unsigned int deb=0; deb<nonZeros.size(); deb++ ) {
						REPORT(DEBUG, "	nonZeros["<<deb<<"]="<<nonZeros[deb]);
					}
					sopcs.push_back(new FixSOPC( target_, msbInC, lsbInC, *msbInter.begin(), *lsbInter.begin(), nonZeros));

					addSubComponent(*sopcs.rbegin());
					msbInter.erase(msbInter.begin());
					lsbInter.erase(lsbInter.begin());

					//wiring
					i=0;
					for ( ; i<indstS; i++) {
					REPORT(3, "inporting T, with parameters X"<<i<<" and T"<<coefIndst.front());
						inPortMap(*sopcs.rbegin(), join("X", i), join("T", coefIndst.front()));
						coefIndst.pop();
					}

					for ( ; i-indstS<indsxS; i++) {
					REPORT(3, "inporting X, with parameters X"<<i<<" and X"<<coefIndsx.front());
						inPortMap(*sopcs.rbegin(), join("X", i), join("X", coefIndsx.front()));
						coefIndsx.pop();
					}
					for ( ; i-( indstS + indsxS )<indsuS; i++) {
					REPORT(3, "inporting Ulv, with parameters X"<<i<<" and Ulv"<<coefIndsu.front());
						inPortMap(*sopcs.rbegin(), join("X", i), join("Ulv", coefIndsu.front()));
						coefIndsu.pop();
					}
					REPORT(3, "inport map done"<<i);

					if (sopcs.size() <= nt) {
					REPORT(3, "outporting T, with parameters R"<<" and T"<<sopcs.size()-1);
						outPortMap(*sopcs.rbegin(), "R", join("T", sopcs.size()-1),0);
						vhdl << instance(*sopcs.rbegin(), join("fixSOPC_t",sopcs.size()-1));
					}
					else if (sopcs.size()<=nt+nx){
					REPORT(3, "outporting Xplus, with parameters R"<<" and Xplus"<<sopcs.size()-nt-1); //sopcs.size()-nt-1 because X starts from 0
						outPortMap(*sopcs.rbegin(), "R", join("Xplus", sopcs.size()-nt-1),0);
						vhdl << instance(*sopcs.rbegin(), join("fixSOPC_y",sopcs.size()));
					}
					else if (sopcs.size()<=nt+nx+ny){
					REPORT(3, "outporting Y, with parameters R"<<" and Y"<<sopcs.size()-nt-nx-1);
						outPortMap(*sopcs.rbegin(), "R", join("Ylv", sopcs.size()-(nt+nx)-1),0); //spocs.size()-(nt+nx)-1 because Y starts to 0
						vhdl << instance(*sopcs.rbegin(), join("fixSOPC_y",sopcs.size()));
						vhdl << tab << join("Y",sopcs.size()-nt-nx-1) << " <= " << join("signed(Ylv",sopcs.size()-nt-nx-1) << ");" << std::endl;
					}
					else {
						THROWERROR("Unexpected error happened during internal building");
						}
					REPORT(3, "outport std::map done");
					


					coefDel.erase(coefDel.begin());
					
					REPORT(0, "finished line. Beginning step "<<coefDel.size()<<", coeffs.empty()="<<coefDel.empty());
		}
		for (uint32_t i = 0; i < nx; i++){
			std::stringstream sig;
			sig<<"Xplus"<<i;
			REPORT(0,"export 	X"<<i<<" <= "<<delay(sig.str(),1)<<";");
			vhdl << tab << join("X",i)<<" <= "<<delay(sig.str(),1)<<";"<<std::endl;
		}
		REPORT(0, "time shifting done");
	//TODO: end of new code
		
	};

	FixSIF::~FixSIF(){

	};

	void FixSIF::emulate(TestCase * tc){

		std::vector<mpz_class> uValues(nu);
		std::vector<mpfr_t> uVals(nu);
		std::vector<mpfr_t> t(nt);
		mpfr_t coeff, resMul, result;
		mpfr_init2(resMul, veryLargePrec);
		mpfr_init2(result, veryLargePrec);
		mpfr_init2(coeff, veryLargePrec);
		for (uint32_t i=0; i<nt; i++) {
			mpfr_init2(t[i], veryLargePrec);
		}
		for ( uint32_t i=0;i<nu; i++)
			mpfr_init2(uVals[i], veryLargePrec);

		std::stringstream uValGot, lsbInGot;

		uValGot<<"uValues ("<<uValues.size()<< ") = [";
		lsbInGot<<"lsbIn ("<<uValues.size()<< ") = [";
		for (uint32_t i = 0; i<nu; i++){
				uValues[i]=tc->getInputValue(join("U",i));
				uValGot << " " << uValues[i];
				lsbInGot << " " << lsbIn[i];
			mpfr_set_z( uVals[i],uValues[i].get_mpz_t(),GMP_RNDD );
			mpfr_mul_2si( uVals[i], uVals[i], lsbIn[i], GMP_RNDD );
		//std::cout<<" FixSIF: uVals["<<i<<"]=";
			//mpfr_printf( "%.65RNf\n", uVals[i]);
				
		}
		uValGot << " ]";
		lsbInGot << " ]";
		REPORT(0, uValGot.str());
		REPORT(3, lsbInGot.str());

		for (uint32_t l=0; l<nt+nx+ny; l++){
			mpfr_set_ui(result,0, GMP_RNDN);
			for (uint32_t c=0; c<nt+nx+nu; c++){
		//std::cout<<"GETTING A NEW COEFF"<<std::endl;

				sollya_obj_t stringCoeff;
				//get coeff as sollya object
				stringCoeff = sollya_lib_parse_string(coeffs[l][c].c_str());
				// If conversion did not succeed (i.e. parse error)
				if(stringCoeff == 0)	{ std::ostringstream error;
						error << srcFileName << ":"<<"emulate(): Unable to parse std::string " << coeffs[l][c] << " as a numeric constant" << std::endl;
						throw error.str();
				}
				//get coeff as mpfr
				sollya_lib_get_constant(coeff, stringCoeff);
				sollya_lib_clear_obj(stringCoeff);
				if (coeff != 0){ //if coeff is not null, compute the product and add it to the result for the line

						if (c>0 && c<nt) {
						}
						if (c<nt){ //if multiplicand is a t
							if ( ( (l>=nt) || (l>c) ) ) {
							 //ones on the J diagonal are implicit
								REPORT(3, "computing t, c="<<c<<",l="<<l);
								mpfr_mul(resMul, t[c], coeff, GMP_RNDN);
								mpfr_sub(result, result, resMul, GMP_RNDN);
							}
						}
						else if (c<nt+nx){ //if multiplicand is a x
							mpfr_mul(resMul, xPrec[c-nt], coeff, GMP_RNDN);
							mpfr_add(result, result, resMul, GMP_RNDN);

						}
						else { //if not t and not x, multiplicand is a u

							mpfr_mul(resMul, uVals[c-nt-nx], coeff, GMP_RNDN);
							mpfr_add(result, result, resMul, GMP_RNDN);
						}
				}
			}
			//for ( unsigned int i=0; i<nx; i++) {
			//			std::cout<<"FixSIF:emulate: xPrec["<<i<<"]=";
			//		mpfr_out_str(NULL, 10, 0, xPrec[i], GMP_RNDN);
			//			std::cout<<std::endl;
			//}
			//for ( unsigned int i=0; i<nx; i++) {
			//			std::cout<<"FixSIF:emulate: xCurrent["<<i<<"]=";
			//		mpfr_out_str(NULL, 10, 0, xCurrent[i], GMP_RNDN);
			//			std::cout<<std::endl;
			//}
			//for ( unsigned int i=0; i<nt; i++) {
			//			std::cout<<"FixSIF:emulate: t["<<i<<"]=";
			//		mpfr_out_str(NULL, 10, 0, t[i], GMP_RNDN);
			//			std::cout<<std::endl;
			//}
			if (l<nt){
				mpfr_set(t[l],result, GMP_RNDN);
			}
			else if (l<nt+nx) {
				mpfr_set(xCurrent[l-nt], result, GMP_RNDN);
			}
			else if (l>=nx+nt){
				//mpfr_t rd, ru;
				//mpfr_init2 (rd, veryLargePrec);
				//mpfr_init2 (ru, veryLargePrec);		
				mpz_class rdz, ruz;

				//mpfr_get_z (rdz.get_mpz_t(), result, GMP_RNDD); 					// there can be a real rounding here

				mpfr_get_z (rdz.get_mpz_t(), result, GMP_RNDD); 					// there can be a real rounding here

				rdz=signedToBitVector(rdz, msbOut[l-(nt+nx)]-lsbOut[l-(nt+nx)]);
//TODO: here we should have computed msbIns
//TODO: replace lsbIn by MsbIn

				mpfr_get_z (ruz.get_mpz_t(), result, GMP_RNDU); 					// there can be a real rounding here	
				ruz=signedToBitVector(ruz, msbOut[l-(nt+nx)]-lsbOut[l-(nt+nx)]);
				
				//mpfr_clears ( result, rd, ru, NULL);
				mpfr_clears ( result, NULL);

				REPORT(0, "adding expected output for Y"<<l-nt-nx<<":"<<rdz);
				REPORT(0, "adding expected output for Y"<<l-nt-nx<<":"<<ruz);

				tc->addExpectedOutput(join("Y",l-nt-nx), rdz);
				tc->addExpectedOutput(join("Y",l-nt-nx), ruz);
			}

		}
		for ( uint32_t i=0; i<nx; i++ ) {
			mpfr_set( xPrec[i], xCurrent[i], GMP_RNDN );
		}

		for( uint32_t i=0;i<nt; i++ ) {
			mpfr_clears(t[i], NULL);
		}
		for( uint32_t i=0;i<nu; i++ ) {
			mpfr_clears(uVals[i], NULL);
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

			REPORT(3,"Got header at line "<<lc<<": \""<<line<<"\"")
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
					REPORT(3,"Popped back to readMatrices");
		}
		return 0;



	}

	int FixSIF::readMatrix(std::string &header, std::string JKLMNPQRS, std::vector < std::vector <std::string> > * &toFill, std::ifstream &openedFile, int lc){

			std::string line;

			int nc=0;
			int nl=0;
			REPORT(3, "Reading matrix given the header: \""<<header<<"\". Matrix descriptor: "<<JKLMNPQRS);
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
				REPORT(3, nc<<" columns and "<<nl<<" lines found in matrix from header: \""<<header<<"\"");
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
				REPORT(3,"Got line "<<lc<<": \""<<line<<"\"")

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
			REPORT(3, "Parsing matrices");
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
		std::cout<<"exiting bMToDoubleM"<<std::endl;	
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
				//std::cout<<"iteration: i="<<i<<", j="<<j<<std::endl;
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

		ublas::matrix<double>C1(bL.size1(), bJ.size1());
		ublas::matrix<double>C2(bR.size1(), bR.size2());

		ublas::matrix<double>D1(bL.size1(), bJ.size1());
		ublas::matrix<double>D2(bS.size1(), bS.size2());

		std::cout<<"computeABCD false for now"<< std::endl;

		std::cout<<"computing KJ^(-1)M+P"<< std::endl;
		if (!matMul(A1, bK, invbJ))
			return false;
		if (!matMul(A2, A1, bM   ))
			return false;
		if (!matAdd(bA, A2, bP   ))
			return false;
                                  
		std::cout<<"computing KJ^(-1)N+Q"<< std::endl;
		if (!matMul(B1, bK, invbJ))
			return false;
		if (!matMul(B2, B1, bN   ))
			return false;
		if (!matAdd(bB, B2, bQ   ))
			return false;
                                  
		std::cout<<"computing LJ^(-1)M+R"<< std::endl;
		if (!matMul(C1, bL, invbJ))
			return false;
		if (!matMul(C2, C1, bM   ))
			return false;
		if (!matAdd(bC, C2, bR   ))
			return false;
                                  
		std::cout<<"computing LJ^(-1)N+S"<< std::endl;
		if (!matMul(D1, bL, invbJ))
			return false;
		if (!matMul(D2, D1, bN   ))
			return false;
		if (!matAdd(bD, D2, bS   ))
			return false;

		displayBMatrix(bJ, "computeABCD: J");
		displayBMatrix(bK, "computeABCD: K");
		displayBMatrix(bL, "computeABCD: L");
		displayBMatrix(bM, "computeABCD: M");
		displayBMatrix(bN, "computeABCD: N");
		displayBMatrix(bP, "computeABCD: P");
		displayBMatrix(bQ, "computeABCD: Q");
		displayBMatrix(bR, "computeABCD: R");
		displayBMatrix(bS, "computeABCD: S");

		displayBMatrix(invbJ, "computeABCD: J^(-1)");
		
		displayBMatrix(A1, "computeABCD: KJ^(-1)");
		displayBMatrix(A2, "computeABCD: KJ^(-1)M");
		displayBMatrix(bA, "computeABCD: KJ^(-1)M+P");
		
		displayBMatrix(B1, "computeABCD: KJ^(-1)");
		displayBMatrix(B2, "computeABCD: KJ^(-1)N");
		displayBMatrix(bB, "computeABCD: KJ^(-1)N+Q");
		
		displayBMatrix(C1, "computeABCD: LJ^(-1)");
		displayBMatrix(C2, "computeABCD: LJ^(-1)M");
		displayBMatrix(bC, "computeABCD: LJ^(-1)M+R");
		
		displayBMatrix(D1, "computeABCD: LJ^(-1)");
		displayBMatrix(D2, "computeABCD: LJ^(-1)N");
		displayBMatrix(bD, "computeABCD: LJ^(-1)N+S");

		std::cout<<"exiting computeABCD ok"<<std::endl;
		
		return true;
	}
	
	int FixSIF::computeMSBsLSBs(){
		//A, B, C, D are double pointers
		//J ,K, L, M, N, P, Q, R, and S are boost double matrices
		//prefix b in bX denotes a boost matrix
		//First: computing MSBs

		std::cout<<"entering computeMSBsLSBs"<<std::endl;
#ifdef WCPG_F
		for ( unsigned int i=0; i<nt+nx+ny; i++) {
			lsbInter.push_back(DEFAULT_LSB);
			msbInter.push_back(DEFAULT_MSB);
		}
		return 0;
#endif

#ifndef WCPG_F
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

		displayBMatrix(Z, "Z");
		displayBMatrix(J, "J");
		displayBMatrix(K, "K");
		displayBMatrix(L, "L");
		displayBMatrix(M, "M");
		displayBMatrix(N, "N");
		displayBMatrix(P, "P");
		displayBMatrix(Q, "Q");
		displayBMatrix(R, "R");
		displayBMatrix(S, "S");

		displayBMatrix(bA, "bA");
		displayBMatrix(bB, "bB");
		displayBMatrix(bC, "bC");
		displayBMatrix(bD, "bD");

		//Transfert into double* format
		bMToDoubleM( bA, A);
		bMToDoubleM( bB, B);
		bMToDoubleM( bC, C);
		bMToDoubleM( bD, D);

		//Declare wcpg of the filter
		double *wcpgF=new double [(ny)*(nu)];
		//double *wcpgF=(double*) malloc( sizeof(double)*(nx+ny)*(nx+nu) );

		//Compute WCPG
		//WCPG_ABCD( wcpgF, A, B, C, D, uint64_t(nx), uint64_t(nu), uint64_t(ny) );
		if ( !WCPG_ABCD( wcpgF, A, B, C, D, nx, nu, ny ) ) {
			THROWERROR("Unable to compute filter WCPG");
		}

		//displaying wcpg
		int nanErr=0;
		std::cout<<"computateABCD: nu="<<nu<<std::endl;
		for ( unsigned int i=0; i<ny; i++ ) {
			std::ostringstream wcpgFDisp;
			wcpgFDisp<<"wcpgF["<<i<<"]:";
			for ( unsigned int j=0; j<nu; j++ ) {
				wcpgFDisp<<" "<<wcpgF[i*nu+j];
				if ( isnan((wcpgF[i*(nu)+j])) )
					nanErr++;
			}
			REPORT(0,wcpgFDisp.str());
		}
		if ( nanErr ){
			REPORT(0,"CRITICAL ERROR: One of the wcpg coefficients is a nan. Any further computation would be garbage. Aborting here.");
			throw "Aborted.";
		}

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


		REPORT(0,"MSBs stored.");

		//Declare pseudo identities
		ublas::matrix<double> Mt(nt,nt+nx+ny);
		ublas::matrix<double> Mx(nx,nt+nx+ny);
		ublas::matrix<double> My(ny,nt+nx+ny);

		REPORT(0,"Initializing Mt matrix");
		//Init Mt
		for ( unsigned int i=0; i<nt; i++ ) {
			for ( unsigned int j=0; j<nt+nx+ny; j++ ) {
				if ( (j<nt) && (i==j) ) {
					Mt(i,j)=1;
				}
				else {
					Mt(i,j)=0;
				}
			}
		}

		REPORT(0,"Initializing Mx matrix");
		//Init Mx
		for ( unsigned int i=0; i<nx; i++ ) {
			for ( unsigned int j=0; j<nt+nx+ny; j++ ) {
				if ( (j>=nt) && (j<nt+nx) && (i==(j-nt)) ) {
					Mx(i,j)=1;
				}
				else {
					Mx(i,j)=0;
				}
			}
		}

		REPORT(0,"Initializing My matrix");
		//Init My
		for ( unsigned int i=0; i<ny; i++ ) {
			for ( unsigned int j=0; j<nt+nx+ny; j++ ) {
				if ( (j>=nt+nx) && (i==j) ) {
					My(i,j)=1;
				}
				else {
					My(i,j)=0;
				}
			}
		}

		//Declare boost ABCD representation for the error filter
		ublas::matrix<double> bAe(nx,nx);
		ublas::matrix<double> bBe(nx,nt+nx+ny);
		ublas::matrix<double> bCe(ny,nx);
		ublas::matrix<double> bDe(ny,nt+nx+ny);

		//Compute ABCD error filter
		computeABCD(J, K, L, M, Mt, P, Mx, R, My, bAe, bBe, bCe, bDe);

		displayBMatrix(J, "J");
		displayBMatrix(K, "K");
		displayBMatrix(L, "L");
		displayBMatrix(M, "M");
		displayBMatrix(Mt, "Mt");
		displayBMatrix(P, "P");
		displayBMatrix(Mx, "Mx");
		displayBMatrix(R, "R");
		displayBMatrix(My, "My");

		displayBMatrix(bAe, "bAe");
		displayBMatrix(bBe, "bBe");
		displayBMatrix(bCe, "bCe");
		displayBMatrix(bDe, "bDe");

		//Init ABCD for error filter
		double *Ae = new double[nx*nx] ;
		double *Be = new double[nx*nt+nx+ny] ;
		double *Ce = new double[ny*nx] ;
		double *De = new double[ny*nt+nx+ny] ;

		//Transfert into double* format
		bMToDoubleM( bAe, Ae );
		bMToDoubleM( bBe, Be );
		bMToDoubleM( bCe, Ce );
		bMToDoubleM( bDe, De );

		//displaying Ae
		for ( unsigned int i=0; i<nx; i++ ) {
			std::ostringstream doubleMDisp;
			doubleMDisp<<"Ae["<<i<<"]:";
			for ( unsigned int j=0; j<nx; j++ ) {
				doubleMDisp<<" "<<Ae[i*(nx)+j];
			}
			REPORT(0,doubleMDisp.str());
		}
		
		//displaying Be
		for ( unsigned int i=0; i<nx; i++ ) {
			std::ostringstream doubleMDisp;
			doubleMDisp<<"Be["<<i<<"]:";
			for ( unsigned int j=0; j<nt+nx+ny; j++ ) {
				doubleMDisp<<" "<<Be[i*(nt+nx+ny)+j];
			}
			REPORT(0,doubleMDisp.str());
		}

		//displaying Ce
		for ( unsigned int i=0; i<ny; i++ ) {
			std::ostringstream doubleMDisp;
			doubleMDisp<<"Ce["<<i<<"]:";
			for ( unsigned int j=0; j<nx; j++ ) {
				doubleMDisp<<" "<<Ce[i*(nx)+j];
			}
			REPORT(0,doubleMDisp.str());
		}
		
		//displaying De
		for ( unsigned int i=0; i<ny; i++ ) {
			std::ostringstream doubleMDisp;
			doubleMDisp<<"De["<<i<<"]:";
			for ( unsigned int j=0; j<nt+nx+ny; j++ ) {
				doubleMDisp<<" "<<De[i*(nt+nx+ny)+j];
			}
			REPORT(0,doubleMDisp.str());
		}

		//Declare wcpgE of the error filter
		double *wcpgE=new double [(nt+nx+ny)*(ny)];

		//Compute WCPG error
		if ( !WCPG_ABCD( wcpgE, A, B, C, D, nx, nt+nx+ny, ny ) ) {
			THROWERROR("Unable to compute wcpg for error filter.");
		}

		//displaying wcpg
		nanErr=0;
		for ( unsigned int i=0; i<ny; i++ ) {
			std::ostringstream wcpgEDisp;
			wcpgEDisp<<"wcpgE["<<i<<"]:";
			for ( unsigned int j=0; j<nt+nx+ny; j++ ) {
				wcpgEDisp<<" "<<wcpgE[i*(nt+nx+ny)+j];
				if ( isnan((wcpgE[i*(nt+nx+ny)+j])) )
					nanErr++;
			}
			REPORT(0,wcpgEDisp.str());
		}
		if ( nanErr )
			THROWERROR("One of the wcpg coefficients is a nan. Any further computation would be garbage. Aborting here.");



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
		   lsbInter.push_back(ceil(log(prim[i])/log(2)));
	   }
	   //real toto=0.7;
	   //toto+=41.3;



		



		std::cout<<"exiting computeMSBsLSBs"<<std::endl;
		return 0;
#endif
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
	

	void FixSIF::displayBMatrix( ublas::matrix<double> const &M, std::string const name ) const {
		for ( unsigned int i=0; i<M.size1(); i++ ) {
			std::ostringstream display;
			display<<name<<"["<<i<<"]:";
			for ( unsigned int j=0; j<M.size2(); j++ ) {
				display<<" "<<M(i,j);
			}
			REPORT(0, display.str());
		}
	}
}
