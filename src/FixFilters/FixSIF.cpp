#include <iostream>
#include <sstream>
#include <queue>

#include "gmp.h"
#include "mpfr.h"

#include "FixSIF.hpp"

#include "ShiftReg.hpp"
#include "FixSOPC.hpp"

#define BUF_SIZE 256 //default buffer size
#define DEFAULT_PREC -8 //default precision for sopcs

using namespace std;

int displayMatrix( vector< vector<string> > &m, string name = "matrix"){
	cerr<<"Beginning matrix display"<<endl;
	cerr<<name<<" "<<m.size()<<" "<<m[0].size()<<endl;
	for (int i=0; i<m.size(); i++){
		for (int j=0; j<m[i].size(); j++){
			cerr<<m[i][j];
			if (j!=m[i].size()-1){
				cerr<<"	";
			}
			else
				cerr<<endl;

		}
		cerr<<endl;
	}
	return 0;
}

//TODO: build a policy to detect suitable use of shiftRegisters.
namespace flopoco {

	FixSIF::FixSIF(Target* target, int lsb_, string file_, map<string, double> inputDelays) : 
		Operator(target), p(-lsb_), file(file_)
	{

		srcFileName="FixSIF";
		setCopyrightString ( "Antoine Martinet, Florent de Dinechin (2014)" );

		ostringstream name;
		name << "FixSIF_" << p << "_" << coeffs.size() << "_uid" << getNewUId();
		setNameWithFreq( name.str() );

		parseFile();
		buildVHDL();
	};

	

	// The method that does the work once coeffs[][] is known
	void FixSIF::buildVHDL(){
//	//TODO: beginning of new code
		n = nt + nx + nu; //TODO: update this definition with a clever number
		vector <int> lsbIn, lsbOut;
		readPrecision(lsbIn, lsbOut, 0);

		for ( int i = 0; i<lsbOut.size(); i++)
		{
			cout<<"lsbOut["<<i<<"]="<<lsbOut[i]<<endl;
		}
		cout<<"nt="<<nt<<", nx="<<nx<<", ny="<<ny<<", nu="<<nu<<endl;
		useNumericStd_Unsigned();
		if(p<1) {
			THROWERROR("Can't build an architecture for this value of LSB")
		}

		//declare intermediate T(k+1)
		for (int i=0; i<nt; i++) {
			vhdl<< tab << declare(join("T",i));
		}

		//nu inputs numbering from 0 to nu-1 FIXME: check the relevance of this convention dealing with SIF
		for (int i=0; i<nu; i++) {
			ostringstream input;
			input << "U" <<i;
			addInput(input.str(), 1+p, true);
		}

		//ny outputs numbering from 0 to ny-1 FIXME: check the relevance of this convention dealing with SIFs
		for (int i=0; i<ny; i++) {
			ostringstream output;
			output << "Y" << i;
			addOutput(output.str(), 1+p, true);
		}
	
		//declare intermediate X and X(k+1)
		for (int i=0; i<nx; i++) {
			//vhdl<< tab << declare(join("Xplus",i));
			vhdl<< tab << declare(join("X",i));
		}

		int txy = 0; //integer to determine wether we are comcomputing an
					//intermediate variable, a state or an output
		vector <FixSOPC*> sopcs;//(n);

		while ( !coeffs.empty() ) {

			vector <string> nonZeros(n);
			queue <int> coefIndst;
			queue <int> coefIndsx;
			queue <int> coefIndsu;

			cout<<"copy first line of coeffs"<<endl;
			//copy first line of coeffs
			for (int i = 0; i<coeffs[0].size(); i++){
				nonZeros[i] = coeffs[0][i];
				cout<<"coeffs[0]["<<i<<"]="<<coeffs[0][i]<<endl;
				cout<<"	nonZeros["<<i<<"]="<<nonZeros[i]<<", nonZeros.size()="<<nonZeros.size()<<endl;
			}
			//output matrix vector
			cout<<"nonZeros feed"<<endl;
			for (int i = 0; i<nonZeros.size(); i++)
				cout<<"	nonZeros["<<i<<"]="<<nonZeros[i]<<endl;


			cout<<"cleaning t from implicit ones and zeros"<<endl;

			//three steps (t, x and u) because all are not treated the same way
			//drop zeros and (implicit) ones in the diagonal but keep indices for consistency
			int ib = 0; //bias in vector indices induced by erase operations
			int i = 0;
			for ( ; i<nt; i++ )
			{
				if ( ( stof(nonZeros[i+ib].c_str()) == 0.0 ) || ( (stof(nonZeros[i+ib].c_str()) == 1.0)&&(n-coeffs.size()==i)) ){
					nonZeros.erase(nonZeros.begin()+i+ib);
					cout<<"deleting nonZeros["<<i<<"]"<<endl;
					ib--;
				}
				else {
					cout<<"pushing "<<i<<" in coefIndst"<<endl;
					coefIndst.push(i);
				}
			}

			cout<<"cleaning x from implicit zeros, nonZeros of size:"<<nonZeros.size()<<endl;
					cout<<"ib="<<ib<<", i="<<i<<endl;
			
			//output matrix vector
			for (int i = 0; i<nonZeros.size(); i++)
				cout<<"	nonZeros["<<i<<"]="<<nonZeros[i]<<endl;

			//only drop zeros for x
			for ( ; i<nt+nx; i++ )
			{
				if ( stof(nonZeros[i+ib].c_str()) == 0.0 ){
			cout<<"no coef found at "<<i+ib<<endl;
			cout<<"numerical value of nonZeros["<<i+ib<<"]="<<stof(nonZeros[i+ib].c_str())<<endl;
					nonZeros.erase(nonZeros.begin()+i+ib);
					ib--;
				}
				else{
					coefIndsx.push(i-nt);
			cout<<"coef found at "<<i+ib<<endl;
			cout<<"numerical value of nonZeros["<<i+ib<<"]="<<stof(nonZeros[i+ib].c_str())<<endl;
					cout<<"pushing "<<i-nt<<" in coefIndsx"<<endl;
				}
			}

			cout<<"cleaning u from implicit zeros, ib="<<ib<<", i="<<i<<", nu+nx+nt="<<nu+nx+nt<<endl;

			//output matrix vector
			for (int i = 0; i<nonZeros.size(); i++)
				cout<<"	nonZeros["<<i<<"]="<<nonZeros[i]<<endl;

			//and for u
			for ( ; i<nu+nx+nt; i++ )
			{
				if ( stof(nonZeros[i+ib].c_str()) == 0.0 ){
			cout<<"no coef found at "<<i+ib<<endl;
			cout<<"numerical value of nonZeros["<<i+ib<<"]="<<stof(nonZeros[i+ib].c_str())<<endl;
					nonZeros.erase(nonZeros.begin()+i+ib);
					ib--;
				}
				else{
					coefIndsu.push(i-nt-nx);
					cout<<"pushing "<<i-nt-nx<<" in coefIndsu"<<endl;
				}
			}

			//if only one coeff, just write a multiplier, if only ones, just write an adder
			//if ( nonZeros.size() > 1 ) {
					
				cout<<"trying to push SOPC with following parameters:"<<endl;
				cout<<"	lsbOut="<<*lsbOut.begin()<<endl;
				cout<<"	nonZeros.size()="<<nonZeros.size()<<endl;
				for ( int i = 0; i<nonZeros.size(); i++)
				{
					cout<<"	nonZeros["<<i<<"]=";
					cout<<nonZeros[i]<<endl;
				}
				cout<<"starting push"<<endl;

					sopcs.push_back(new FixSOPC( target_, *lsbIn.begin(), *lsbOut.begin(), nonZeros));
				cout<<"SOPC pushed back"<<endl;

					addSubComponent(*sopcs.rbegin());
					lsbIn.erase(lsbIn.begin());
					lsbOut.erase(lsbOut.begin());

					cout<<"coefIndst.size()="<<coefIndst.size()<<endl;
					cout<<"coefIndsx.size()="<<coefIndsx.size()<<endl;
					cout<<"coefIndsu.size()="<<coefIndsu.size()<<endl;

					//wiring
					i=0;
					for ( ; i<coefIndst.size(); i++) {
					cout<<"inporting T, with parameters X"<<i<<" and T"<<coefIndst.front()<< endl;
						inPortMap(*sopcs.rbegin(), join("X", i), join("T", coefIndst.front()));
						coefIndst.pop();
					}

					for ( ; i-coefIndst.size()<coefIndsx.size(); i++) {
					cout<<"inporting X, with parameters X"<<i<<" and X"<<coefIndsx.front()<< endl;
						inPortMap(*sopcs.rbegin(), join("X", i), join("X", coefIndsx.front()));
						coefIndsx.pop();
					}
					for ( ; i-( coefIndst.size() + coefIndsx.size() )<coefIndsu.size(); i++) {
					cout<<"inporting U, with parameters X"<<i<<" and U"<<coefIndsu.front()<< endl;
						inPortMap(*sopcs.rbegin(), join("X", i), join("U", coefIndsu.front()));
						coefIndsu.pop();
					}
					cout<<"inport map done"<<endl;

					if (sopcs.size() <= nt) {
					cout<<"outporting T, with parameters R"<<" and T"<<sopcs.size()<< endl;
						outPortMap(*sopcs.rbegin(), "R", join("T", sopcs.size()));
						vhdl << instance(*sopcs.rbegin(), join("fixSOPC_t",sopcs.size()));
					}
					else if (sopcs.size()<=nt+nx){
					cout<<"outporting Xplus, with parameters R"<<" and Xplus"<<sopcs.size()<< endl;
						outPortMap(*sopcs.rbegin(), "R", join("Xplus", sopcs.size()-nt));
						vhdl << instance(*sopcs.rbegin(), join("fixSOPC_y",sopcs.size()));
					}
					else if (sopcs.size()<=nt+nx+ny){
					cout<<"outporting Y, with parameters R"<<" and Y"<<sopcs.size()<< endl;
						outPortMap(*sopcs.rbegin(), "R", join("Y", sopcs.size()-(nt+nx)));
						vhdl << instance(*sopcs.rbegin(), join("fixSOPC_y",sopcs.size()));
					}
					else {
						THROWERROR("Unexpected error happened during internal building");
						}
					cout<<"outport map done"<<endl;
					
					for (int i = 0; i < nx; i++){
						stringstream sig;
						sig<<"Xplus"<<i;
						vhdl<<"X"<<i<<" <= "<<delay(sig.str(),1);
					}
					cout<<"time shifting done"<<endl;


					coeffs.erase(coeffs.begin());
					
					cout<<"finished. beginning step "<<coeffs.size();
						cout<<", coeffs.empty()="<<coeffs.empty()<<endl;
		}
		cout<<"exited while loop"<<endl;
	//TODO: end of new code
		
};

	FixSIF::~FixSIF(){

	};

	void FixSIF::emulate(TestCase * tc){

#if 1
		mpz_class sx;

		sx = tc->getInputValue("X"); 		// get the input bit vector as an integer
		xHistory[currentIndex] = sx;

		mpfr_t x, t, s, rd, ru;
		mpfr_init2 (x, 1+p);
		mpfr_init2 (t, 10*(1+p));
		mpfr_init2 (s, 10*(1+p));
		mpfr_init2 (rd, 1+p);
		mpfr_init2 (ru, 1+p);		
		mpfr_set_d(s, 0.0, GMP_RNDN); // initialize s to 0
		
		for (int i=0; i< n; i++)	{
			sx = xHistory[(currentIndex+n-i)%n];
			sx = bitVectorToSigned(sx, 1+p); 						// convert it to a signed mpz_class
			mpfr_set_z (x, sx.get_mpz_t(), GMP_RNDD); 	 // convert this integer to an MPFR; this rounding is exact
			mpfr_div_2si (x, x, p, GMP_RNDD); 						// multiply this integer by 2^-p to obtain a fixed-point value; this rounding is again exact
			
			mpfr_mul(t, x, mpcoeff[i], GMP_RNDN); 					// Here rounding possible, but precision used is ridiculously high so it won't matter

			if(coeffsign[i]==1)
				mpfr_neg(t, t, GMP_RNDN); 
			
			mpfr_add(s, s, t, GMP_RNDN); 							// same comment as above
			}

			// now we should have in s the (exact in most cases) sum
			// round it up and down

			// make s an integer -- no rounding here
			mpfr_mul_2si (s, s, p, GMP_RNDN);

			mpz_class rdz, ruz;

			mpfr_get_z (rdz.get_mpz_t(), s, GMP_RNDD); 					// there can be a real rounding here
			rdz=signedToBitVector(rdz, wO);
			tc->addExpectedOutput ("R", rdz);
			// tc->addExpectedOutput ("R", rdz);

			mpfr_get_z (ruz.get_mpz_t(), s, GMP_RNDU); 					// there can be a real rounding here	
			ruz=signedToBitVector(ruz, wO);
			tc->addExpectedOutput ("R", ruz);

			mpfr_clears (x, t, s, rd, ru, NULL);
			
			currentIndex=(currentIndex+1)%n; //  circular buffer to store the inputs

	
#else
		static int idx = 0;
		static bool full = false; 							// set to true when the fir start to output valid data (after n input) 
		static TestCase * listTC [10000]; // should be enough for everybody


		listTC[idx] = tc;

		if(n == 1)					// if the fir part has only one tap we don't wait to get the output
			full = true; 

		// We are waiting until the first meaningful value comes out of the FIR
		if (full) {
			mpfr_t x, t, s, rd, ru;
			mpfr_init2 (x, 1+p);
			mpfr_init2 (t, 10*(1+p));
			mpfr_init2 (s, 10*(1+p));
			mpfr_init2 (rd, 1+p);
			mpfr_init2 (ru, 1+p);		

			mpfr_set_d(s, 0.0, GMP_RNDN); // initialize s to 0


			int k = idx; // We start to sum from the last input

			for (int i=0; i< n; i++)
			{

				mpz_class sx = listTC[k]->getInputValue("X"); 		// get the input bit vector as an integer
				sx = bitVectorToSigned(sx, 1+p); 						// convert it to a signed mpz_class
				mpfr_set_z (x, sx.get_mpz_t(), GMP_RNDD); 				// convert this integer to an MPFR; this rounding is exact
				mpfr_div_2si (x, x, p, GMP_RNDD); 						// multiply this integer by 2^-p to obtain a fixed-point value; this rounding is again exact

				mpfr_mul(t, x, mpcoeff[i], GMP_RNDN); 					// Here rounding possible, but precision used is ridiculously high so it won't matter

				if(coeffsign[i]==1)
					mpfr_neg(t, t, GMP_RNDN); 

				mpfr_add(s, s, t, GMP_RNDN); 							// same comment as above
			
				k = (k+1)%n;	
			}

			k = (k-1+n)%n; //to get the corresponding testCase to the outputed value

			// now we should have in s the (exact in most cases) sum
			// round it up and down

			// make s an integer -- no rounding here
			mpfr_mul_2si (s, s, p, GMP_RNDN);

			mpz_class rdz, ruz;

			mpfr_get_z (rdz.get_mpz_t(), s, GMP_RNDD); 					// there can be a real rounding here
			rdz=signedToBitVector(rdz, wO);
			listTC[k]->addExpectedOutput ("R", rdz);
			// tc->addExpectedOutput ("R", rdz);

			mpfr_get_z (ruz.get_mpz_t(), s, GMP_RNDU); 					// there can be a real rounding here	
			ruz=signedToBitVector(ruz, wO);
			listTC[k]->addExpectedOutput ("R", ruz);

			mpfr_clears (x, t, s, rd, ru, NULL);
		}
		
		idx = (idx-1+n)%n; // We use a circular buffer to store the inputs

		if (idx ==  1) {
			full = true;
		}

#endif
	};

	int FixSIF::readMatrix(string JKLMNPQRS, vector < vector <string> > * &toFill, ifstream &openedFile, int lc){

			string line;
			getline(openedFile, line);
			lc++;
			//check for empty lines
			if ( line.empty() ) {
				while ( line.empty() && (!openedFile.eof()))
				getline(openedFile, line);
				lc++;
			}

			//check for empty lines filled with spaces or other non senses
			for (int i = 0; i<line.size(); i++) {
				if ( isdigit(line[i]) )
					break;

				if ((!isdigit(line[i])) && (i==line.size()-1)) {
					getline(openedFile, line);
					lc++;
					i=-1;
				}
			}

			int nc=0;
			int nl=0;
			if (line.find_first_of(JKLMNPQRS) != string::npos){
				size_t bNum = line.find_first_of(" ");
				size_t aNum = line.find_first_of(" ",bNum+1);
				if (bNum != string::npos){
					nl = stoi(line.substr(bNum+1,aNum-1));
					size_t bNum = line.find_first_of(" ",aNum);
					if (bNum != string::npos){
						nc = stoi(line.substr(aNum+1,bNum-1));
						toFill = new vector< vector<string> >(nl);
					}
					else {
							stringstream error;
							error<<"line "<<lc<<": space missing";
							THROWERROR(error.str());
					}
					
				}
				else {
						stringstream error;
						error<<"line "<<lc<<": no space found";
						THROWERROR(error.str());
				}
			}
			else{
				stringstream error;
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
			for (int i=0; i<nl; i++){
				getline(openedFile, line);
				lc++;
				size_t bNum = 0;
				size_t aNum = 0;
				for (int j = 0; j<nc; j++) {
					aNum = line.find_first_of("	", bNum);
					if (aNum != string::npos){
						 (*toFill)[i].push_back(line.substr(bNum,aNum-bNum));
						 bNum=aNum+1;
					}

					else if (j==nc-1) {
						 (*toFill)[i].push_back(line.substr(bNum,line.size()-bNum));
					}
					else {
						stringstream error;
						error<<"line "<<lc<<": no tab found";
						THROWERROR(error.str());
					}
				}

			}
			
			return lc;

	}



	int FixSIF::parseFile(){

		int buf_size = ceil(log10((2<<p)-1));
		int lineCounter=0;
		ifstream input;
		input.open( file, ifstream::in );
		vector < vector<string> > *J;
		vector < vector<string> > *K;
		vector < vector<string> > *L;
		vector < vector<string> > *M;
		vector < vector<string> > *N;
		vector < vector<string> > *P;
		vector < vector<string> > *Q;
		vector < vector<string> > *R;
		vector < vector<string> > *S;

		if ( input.is_open() ){

			lineCounter = readMatrix( "J", J, input );
			lineCounter = readMatrix( "K", K, input, lineCounter );
			lineCounter = readMatrix( "L", L, input, lineCounter );
			lineCounter = readMatrix( "M", M, input, lineCounter );
			lineCounter = readMatrix( "P", P, input, lineCounter );
			lineCounter = readMatrix( "R", R, input, lineCounter );
			lineCounter = readMatrix( "N", N, input, lineCounter );
			lineCounter = readMatrix( "Q", Q, input, lineCounter );
			lineCounter = readMatrix( "S", S, input, lineCounter );
			input.close();
		}
		else {
			stringstream error;
			error << "Cannot open file "<<file<<endl;
			THROWERROR(error.str());
		}

		coeffs = vector< vector<string> > (nt+nx+ny);

		int i=0, j;
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
//		displayMatrix(*J,"J");
//		displayMatrix(*K,"K");
//		displayMatrix(*L,"L");
//		displayMatrix(*M,"M");
//		displayMatrix(*P,"P");
//		displayMatrix(*R,"R");
//		displayMatrix(*N,"N");
//		displayMatrix(*Q,"Q");
//		displayMatrix(*S,"S");
//		displayMatrix(coeffs,"coeffs");

	}

	int FixSIF::readPrecision( vector<int> &msbs, vector<int> &lsbs, bool inFile ){
		if ( inFile ){
			vector < vector<string> > *msbslsbs;
			ifstream precisions;
			precisions.open( "precisions.txt", ifstream::in );
			readMatrix( "P", msbslsbs, precisions );
			for (int i=0; i<msbslsbs->size(); i++){
				msbs.push_back(stoi((*msbslsbs)[i][0]));
				lsbs.push_back(stoi((*msbslsbs)[i][1]));
			}
			return 0;
			precisions.close();
		}
		else {
			for ( int i=0; i<nt+nx+ny; i++ ) {
				msbs.push_back(DEFAULT_PREC);
				lsbs.push_back(DEFAULT_PREC);
			}
			return 0;
		}
	}
	
}
	
