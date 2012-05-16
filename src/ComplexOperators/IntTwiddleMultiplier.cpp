#include <fstream>
#include <sstream>
#include "IntTwiddleMultiplier.hpp"

using namespace std;

namespace flopoco{

	extern vector<Operator *> oplist;

	//TODO: explore implementation using multiply-accumulate operators
	//TODO: explore the use of both KCM and Shift-And-Add techniques ->
	//		or just a simple signed multiplier	
	//FIXME: correct timing of the circuit
	//FIXME: verify behaviour in case of negative re/im twiddle parts
	//FIXME: correct the size of the output and intermediary computations
	//		 for now it's fixed at 2*w, achieved through padding (and with 1 all around)
	//FIXME: correct the emulate function
	IntTwiddleMultiplier::IntTwiddleMultiplier(Target* target, int wI_, int wF_, int twiddleExponent_, int n_, bool signedOperator_, bool reducedMultiplications)
		: Operator(target), wI(wI_), wF(wF_), twiddleExponent(twiddleExponent_), n(n_), signedOperator(signedOperator_)
	{
		bool completeExecutionPath = true;
		
		signedOperator ? w = 1 + wI + wF : w = wI + wF;
		
		ostringstream name;

		setCopyrightString ( "Istoan Matei, Florent de Dinechin (2008-2012)" );
		if(target->isPipelined())
			name << "IntTwiddleMultiplier_" << w << "_w_exp_" << twiddleExponent << "_f"<< target->frequencyMHz() << "_uid" << getNewUId();
		else
			name << "IntTwiddleMultiplier_" << w << "_w_exp_" << twiddleExponent << "_uid" << getNewUId();
		setName( name.str() );

		addInput("Xi", 		w, true);
		addInput("Xr", 		w, true);
		if((twiddleExponent==0) || ((double)twiddleExponent == (double)n/4.0))
			addOutput("Zi",   w, 2);
		else
			addOutput("Zi",   2*w, 2);
		if((twiddleExponent==0) || ((double)twiddleExponent == (double)n/4.0))
			addOutput("Zr",   w, 2);
		else
			addOutput("Zr",   2*w, 2);
		
		
		if(twiddleExponent == 0){
						
			vhdl << tab << "Zi <= Xi;" << endl;
			vhdl << tab << "Zr <= Xr;" << endl;
			
			completeExecutionPath = false;
		} else if((double)twiddleExponent == (double)n/4.0){
			
			vhdl << tab << declare("neg_Xi", w) << " <= Xi xor (" << w-1 << " downto 0 => \'1\');" << endl;
			
			IntAdder* addOperator =  new IntAdder(target, w, inDelayMap("X",getCriticalPath()));
			oplist.push_back(addOperator);
			
			inPortMap 	(addOperator, "X", 	 "neg_Xi");
			inPortMapCst(addOperator, "Y", 	 zg(w, 0));
			inPortMapCst(addOperator, "Cin", "\'1\'");
			outPortMap	(addOperator, "R", 	 "Zr", false);
			vhdl << instance(addOperator, "ADD_negXi");
			
			vhdl << tab << "Zi <= Xr;" << endl;
			
			completeExecutionPath = false;
		}
	
		if(!reducedMultiplications && completeExecutionPath){
			mpz_class twRe, twIm;
			int wOutIm, wOutRe;
			
			twRe = getTwiddleConstant(TWIDDLERE);
			twIm = getTwiddleConstant(TWIDDLEIM);
			
			FixRealKCM *multiplyOperatorRe, *multiplyOperatorIm;
			
			multiplyOperatorRe = new FixRealKCM(target, -wF, wI-1, 1, -2*wF, getTwiddleConstantString(TWIDDLERE));
			oplist.push_back(multiplyOperatorRe);
			
			if((n/getGCD(n, 2*twiddleExponent))%4 == 0){
				multiplyOperatorIm = multiplyOperatorRe;
			}else{
				multiplyOperatorIm = new FixRealKCM(target, -wF, wI-1, 1, -2*wF, getTwiddleConstantString(TWIDDLEIM));
				oplist.push_back(multiplyOperatorIm);
			}
			IntAdder* addOperator =  new IntAdder(target, 2*w, inDelayMap("X",getCriticalPath()));
			oplist.push_back(addOperator);
			
			inPortMap (multiplyOperatorIm, "X", "Xi");
			outPortMap(multiplyOperatorIm, "R", "intXiYi");
			vhdl << instance(multiplyOperatorIm, "MUL_XiYi");
			
			inPortMap (multiplyOperatorRe, "X", "Xi");
			outPortMap(multiplyOperatorRe, "R", "intXiYr");
			vhdl << instance(multiplyOperatorRe, "MUL_XiYr");
						
			inPortMap (multiplyOperatorRe, "X", "Xr");
			outPortMap(multiplyOperatorRe, "R", "intXrYr");
			vhdl << instance(multiplyOperatorRe, "MUL_XrYr");
						
			inPortMap (multiplyOperatorIm, "X", "Xr");
			outPortMap(multiplyOperatorIm, "R", "intXrYi");
			vhdl << instance(multiplyOperatorIm, "MUL_XrYi");
			
			wOutIm = 1 + wI + 2*wF + ceil(log2(abs(getTwiddleConstant(TWIDDLEIM).get_si())));
			wOutRe = 1 + wI + 2*wF + ceil(log2(abs(getTwiddleConstant(TWIDDLERE).get_si())));
			
			vhdl << tab << declare("XrYr", 2*w) << " <= (" << 2*w-1-wOutRe << " downto 0 => intXrYr(" << wOutRe-1 << ")) & intXrYr;" << endl;
			vhdl << tab << declare("XiYi", 2*w) << " <= (" << 2*w-1-wOutIm << " downto 0 => intXiYi(" << wOutIm-1 << ")) & intXiYi;" << endl;
			vhdl << tab << declare("XrYi", 2*w) << " <= (" << 2*w-1-wOutIm << " downto 0 => intXrYi(" << wOutIm-1 << ")) & intXrYi;" << endl;
			vhdl << tab << declare("XiYr", 2*w) << " <= (" << 2*w-1-wOutRe << " downto 0 => intXiYr(" << wOutRe-1 << ")) & intXiYr;" << endl;
			
			syncCycleFromSignal("XiYr", false);
			
			std::string strXrYr, strNegXiYi, strXrYi, strXiYr, strCinRe, strCinIm;
			bool negateZr = false, negateZi = false;
			
			if((twIm<0) && (twRe<0)){
				vhdl << tab << declare("neg_XrYi", 2*w) << " <= XrYi xor (" << 2*w-1 << " downto 0 => \'1\');" << endl;
				
				strXrYr 	= "neg_XrYr";
				strNegXiYi 	= "XiYi";
				strCinRe	= "\'1\'";
				negateZr	= false;
				strXrYi		= "XrYi";
				strXiYr		= "XiYr";
				strCinIm	= "\'0\'";
				negateZi	= true;
				
				syncCycleFromSignal("neg_XrYi", false);
				nextCycle();
			}else if((twIm<0) && (twRe>=0)){
				vhdl << tab << declare("neg_XrYi", 2*w) << " <= XrYi xor (" << 2*w-1 << " downto 0 => \'1\');" << endl;
				
				strXrYr 	= "XrYr";
				strNegXiYi 	= "XiYi";
				strCinRe	= "\'0\'";
				negateZr	= false;
				strXrYi		= "neg_XrYi";
				strXiYr		= "XiYr";
				strCinIm	= "\'0\'";
				negateZi	= false;
				
				syncCycleFromSignal("neg_XrYi", false);
				nextCycle();
			}else if((twIm>=0) && (twRe<0)){
				vhdl << tab << declare("neg_XiYr", 2*w) << " <= XiYr xor (" << 2*w-1 << " downto 0 => \'1\');" << endl;
				
				strXrYr 	= "XrYr";
				strNegXiYi 	= "XiYi";
				strCinRe	= "\'0\'";
				negateZr	= true;
				strXrYi		= "XrYi";
				strXiYr		= "neg_XiYr";
				strCinIm	= "\'1\'";
				negateZi	= false;
				
				syncCycleFromSignal("neg_XiYr", false);
				nextCycle();
			}else if((twIm>=0) && (twRe>=0)){
				vhdl << tab << declare("neg_XiYi", 2*w) << " <= XiYi xor (" << 2*w-1 << " downto 0 => \'1\');" << endl;
				
				strXrYr 	= "XrYr";
				strNegXiYi 	= "neg_XiYi";
				strCinRe	= "\'1\'";
				negateZr	= false;
				strXrYi		= "XrYi";
				strXiYr		= "XiYr";
				strCinIm	= "\'0\'";
				negateZi	= false;
				
				syncCycleFromSignal("neg_XiYi", false);
				nextCycle();
			}
			
			inPortMap 	(addOperator, "X", 	 strXrYr);
			inPortMap 	(addOperator, "Y", 	 strNegXiYi);
			inPortMapCst(addOperator, "Cin", strCinRe);
			if(negateZr)
				outPortMap	(addOperator, "R", 	 "intZr");
			else
				outPortMap	(addOperator, "R", 	 "Zr", false);
			vhdl << instance(addOperator, "ADD_XrYrSubXiYi");
			
			inPortMap 	(addOperator, "X", 	 strXrYi);
			inPortMap 	(addOperator, "Y", 	 strXiYr);
			inPortMapCst(addOperator, "Cin", strCinIm);
			if(negateZi)
				outPortMap	(addOperator, "R", 	 "intZi");
			else
				outPortMap	(addOperator, "R", 	 "Zi", false);
			vhdl << instance(addOperator, "ADD_XrYiAddXiYr");
			
			if(negateZr){
				inPortMap 	(addOperator, "X", 	 "intZr");
				inPortMapCst(addOperator, "Y", 	 zg(2*w, 0));
				inPortMapCst(addOperator, "Cin", "\'1\'");
				outPortMap	(addOperator, "R", 	 "Zr", false);
				vhdl << instance(addOperator, "ADD_negZr");
			}
			if(negateZi){
				inPortMap 	(addOperator, "X", 	 "intZi");
				inPortMapCst(addOperator, "Y", 	 zg(2*w, 0));
				inPortMapCst(addOperator, "Cin", "\'1\'");
				outPortMap	(addOperator, "R", 	 "Zi", false);
				vhdl << instance(addOperator, "ADD_negZi");
			}
		}else if(reducedMultiplications && completeExecutionPath){
			try{
			mpz_class twRe, twImAddRe, twImSubRe;
			int wOutRe, wOutImAddRe, wOutImSubRe;
			
			twRe = getTwiddleConstant(TWIDDLERE);
			twImAddRe = getTwiddleConstant(TWIDDLEIMADDRE);
			twImSubRe = getTwiddleConstant(TWIDDLEIMSUBRE);
			
			wOutRe = 1 + wI + 2*wF + ceil(log2(abs(getTwiddleConstant(TWIDDLERE).get_si())));
			wOutImAddRe = 1 + wI + 2*wF + ceil(log2(abs(getTwiddleConstant(TWIDDLEIMADDRE).get_si())));
			wOutImSubRe = 1 + wI + 2*wF + ceil(log2(abs(getTwiddleConstant(TWIDDLEIMSUBRE).get_si())));
			
			FixRealKCM *multiplyOperatorRe, *multiplyOperatorImAddRe, *multiplyOperatorImSubRe;
			
			multiplyOperatorRe = new FixRealKCM(target, -wF, wI-1, 1, -2*wF, getTwiddleConstantString(TWIDDLERE));
			oplist.push_back(multiplyOperatorRe);
			if(abs(twImAddRe.get_si()) == abs(twImSubRe.get_si())){
				multiplyOperatorImAddRe = multiplyOperatorRe;
				multiplyOperatorImSubRe = multiplyOperatorRe;
			}else if(twImAddRe.get_si() == 0){
				multiplyOperatorImSubRe = new FixRealKCM(target, -wF, wI-1, 1, -2*wF, getTwiddleConstantString(TWIDDLEIMSUBRE));
				oplist.push_back(multiplyOperatorImSubRe);
			}else if(twImSubRe.get_si() == 0){
				multiplyOperatorImAddRe = new FixRealKCM(target, -wF, wI-1, 1, -2*wF, getTwiddleConstantString(TWIDDLEIMADDRE));
				oplist.push_back(multiplyOperatorImAddRe);
			}else{
				multiplyOperatorImAddRe = new FixRealKCM(target, -wF, wI-1, 1, -2*wF, getTwiddleConstantString(TWIDDLEIMADDRE));
				oplist.push_back(multiplyOperatorImAddRe);
				multiplyOperatorImSubRe = new FixRealKCM(target, -wF, wI-1, 1, -2*wF, getTwiddleConstantString(TWIDDLEIMSUBRE));
				oplist.push_back(multiplyOperatorImSubRe);
			}
			
			IntAdder* addOperator =  new IntAdder(target, w, inDelayMap("X",getCriticalPath()));
			oplist.push_back(addOperator);	
				
			inPortMap 	(addOperator, "X", 	 "Xr");
			inPortMap 	(addOperator, "Y",   "Xi");
			inPortMapCst(addOperator, "Cin", "\'0\'");
			outPortMap	(addOperator, "R", 	 "XrAddXi");
			vhdl << instance(addOperator, "ADD_XrXi");
			
			syncCycleFromSignal("YrAddYi", false);
			//nextCycle(); 
			
			inPortMap (multiplyOperatorRe, "X", "XrAddXi");
			outPortMap(multiplyOperatorRe, "R", "intK1");
			vhdl << instance(multiplyOperatorRe, "MUL_K1");
			
			if(twImSubRe == 0){
				vhdl << tab << declare("intK2", wOutImSubRe) << " <= " << zg(wOutImSubRe, 0) << ";" << endl;
				
				inPortMap (multiplyOperatorImAddRe, "X", "Xi");
				outPortMap(multiplyOperatorImAddRe, "R", "intintK3");
				vhdl << instance(multiplyOperatorImAddRe, "MUL_K3");
				
				vhdl << tab << declare("intK3", wOutImAddRe+1) << " <= intintK3 & \'0\';" << endl;
			}else if(twImAddRe == 0){
				vhdl << tab << declare("intK3", wOutImAddRe) << " <= " << zg(wOutImAddRe, 0) << ";" << endl;
				
				inPortMap (multiplyOperatorImSubRe, "X", "Xr");
				outPortMap(multiplyOperatorImSubRe, "R", "intintK2");
				vhdl << instance(multiplyOperatorImSubRe, "MUL_K2");
				
				vhdl << tab << declare("intK2", wOutImSubRe+1) << " <= intintK2 & \'0\';" << endl;
			}else{
				inPortMap (multiplyOperatorImSubRe, "X", "Xr");
				outPortMap(multiplyOperatorImSubRe, "R", "intK2");
				vhdl << instance(multiplyOperatorImSubRe, "MUL_K2");
				
				inPortMap (multiplyOperatorImAddRe, "X", "Xi");
				outPortMap(multiplyOperatorImAddRe, "R", "intK3");
				vhdl << instance(multiplyOperatorImAddRe, "MUL_K3");
			}
			
			vhdl << tab << declare("K1", 2*w) << " <= (" << 2*w-1-wOutRe << " downto 0 => intK1(" << wOutRe-1 << ")) & intK1;" << endl;
			if(twImAddRe == 0)
				vhdl << tab << declare("K2", 2*w) << " <= (" << 2*w-2-wOutImSubRe << " downto 0 => intK2(" << wOutImSubRe << ")) & intK2;" << endl;
			else
				vhdl << tab << declare("K2", 2*w) << " <= (" << 2*w-1-wOutImSubRe << " downto 0 => intK2(" << wOutImSubRe-1 << ")) & intK2;" << endl;
			if(twImSubRe == 0)
				vhdl << tab << declare("K3", 2*w) << " <= (" << 2*w-2-wOutImAddRe << " downto 0 => intK3(" << wOutImAddRe << ")) & intK2;" << endl;
			else
				vhdl << tab << declare("K3", 2*w) << " <= (" << 2*w-1-wOutImAddRe << " downto 0 => intK3(" << wOutImAddRe-1 << ")) & intK3;" << endl;
			
			syncCycleFromSignal("K1", false);
			//nextCycle(); 
			
			if(twRe<0){
				vhdl << tab << declare("neg_K1", 2*w) << " <= K1 xor (" << 2*w-1 << " downto 0 => \'1\');" << endl;
			}
			if(twImSubRe<0){
				vhdl << tab << declare("neg_K2", 2*w) << " <= K2 xor (" << 2*w-1 << " downto 0 => \'1\');" << endl;
			}
			vhdl << tab << declare("neg_K3", 2*w) << " <= K3 xor (" << 2*w-1 << " downto 0 => \'1\');" << endl;
			
			syncCycleFromSignal("neg_K3", false);
			//nextCycle();
			
			std::string strK1Im, strK1Re, strK2, strNegK3, strCinRe, strCinIm;
			bool negateZr = false, negateZi = false;
			
			if((twRe<0) && (twImAddRe<0) && (twImSubRe<0)){
				strK1Im = "K1";
				strK2 = "K2";
				strCinIm = "\'0\'";
				negateZi = true;
				strK1Re = "neg_K1";
				strNegK3 = "K3";
				strCinRe = "\'1\'";
				negateZr = false;
			}else if((twRe<0) && (twImAddRe<0) && (twImSubRe>0)){
				strK1Im = "neg_K1";
				strK2 = "K2";
				strCinIm = "\'1\'";
				negateZi = false;
				strK1Re = "neg_K1";
				strNegK3 = "K3";
				strCinRe = "\'1\'";
				negateZr = false;
			}else if((twRe<0) && (twImAddRe>0) && (twImSubRe<0)){
				strK1Im = "K1";
				strK2 = "K2";
				strCinIm = "\'0\'";
				negateZi = true;
				strK1Re = "K1";
				strNegK3 = "K3";
				strCinRe = "\'0\'";
				negateZr = true;
			}else if((twRe<0) && (twImAddRe>0) && (twImSubRe>0)){
				strK1Im = "neg_K1";
				strK2 = "K2";
				strCinIm = "\'1\'";
				negateZi = false;
				strK1Re = "K1";
				strNegK3 = "K3";
				strCinRe = "\'0\'";
				negateZr = true;
			}else if((twRe>0) && (twImAddRe<0) && (twImSubRe<0)){
				strK1Im = "K1";
				strK2 = "neg_K2";
				strCinIm = "\'1\'";
				negateZi = false;
				strK1Re = "K1";
				strNegK3 = "K3";
				strCinRe = "\'0\'";
				negateZr = false;
			}else if((twRe>0) && (twImAddRe<0) && (twImSubRe>0)){
				strK1Im = "K1";
				strK2 = "K2";
				strCinIm = "\'0\'";
				negateZi = false;
				strK1Re = "K1";
				strNegK3 = "K3";
				strCinRe = "\'0\'";
				negateZr = false;
			}else if((twRe>0) && (twImAddRe>0) && (twImSubRe<0)){
				strK1Im = "K1";
				strK2 = "neg_K2";
				strCinIm = "\'1\'";
				negateZi = false;
				strK1Re = "K1";
				strNegK3 = "neg_K3";
				strCinRe = "\'1\'";
				negateZr = false;
			}else if((twRe>0) && (twImAddRe>0) && (twImSubRe>0)){
				strK1Im = "K1";
				strK2 = "K2";
				strCinIm = "\'0\'";
				negateZi = false;
				strK1Re = "K1";
				strNegK3 = "neg_K3";
				strCinRe = "\'1\'";
				negateZr = false;
			}
			
			IntAdder *addOperator2 =  new IntAdder(target, 2*w, inDelayMap("X",getCriticalPath()));
			oplist.push_back(addOperator2);
			
			inPortMap 	(addOperator2, "X",   strK1Re);
			inPortMap 	(addOperator2, "Y",   strNegK3);
			inPortMapCst(addOperator2, "Cin", strCinRe);
			if(negateZr)
				outPortMap	(addOperator2, "R", 	 "intZr");
			else
				outPortMap	(addOperator2, "R", 	 "Zr", false);
			vhdl << instance(addOperator2, "ADD_K1MinK3");
			
			inPortMap 	(addOperator2, "X",   strK1Im);
			inPortMap 	(addOperator2, "Y",   strK2);
			inPortMapCst(addOperator2, "Cin", strCinIm);
			if(negateZr)
				outPortMap	(addOperator2, "R", 	 "intZi");
			else
				outPortMap	(addOperator2, "R", 	 "Zi", false);
			vhdl << instance(addOperator2, "ADD_K1AddK2");
			
			if(negateZr){
				inPortMap 	(addOperator2, "X", 	 "intZr");
				inPortMapCst(addOperator2, "Y", 	 zg(2*w, 0));
				inPortMapCst(addOperator2, "Cin", 	 "\'1\'");
				outPortMap	(addOperator2, "R", 	 "Zr", false);
				vhdl << instance(addOperator2, "ADD_negZr");
			}
			if(negateZi){
				inPortMap 	(addOperator2, "X", 	 "intZi");
				inPortMapCst(addOperator2, "Y", 	 zg(2*w, 0));
				inPortMapCst(addOperator2, "Cin", 	 "\'1\'");
				outPortMap	(addOperator2, "R", 	 "Zi", false);
				vhdl << instance(addOperator2, "ADD_negZi");
			}
			
			}catch(std::string str){
				cout << "execution interrupted: " << str << endl;
				exit(1);
			}
		}
	
	}	


	IntTwiddleMultiplier::~IntTwiddleMultiplier()
	{
	}
	
	//FIXME: correct the emulate function
	void IntTwiddleMultiplier::emulate ( TestCase* tc ) {
		mpz_class svXi = tc->getInputValue("Xi");
		mpz_class svXr = tc->getInputValue("Xr");
		mpz_class svYi = getTwiddleConstant(TWIDDLEIM);
		mpz_class svYr = getTwiddleConstant(TWIDDLERE);
		
		
		if (! signedOperator){

			mpz_class svZi = svXr*svYi + svXi*svYr;
			mpz_class svZr = svXr*svYr - svXi*svYi;
			
			// Don't allow overflow
			mpz_clrbit ( svZi.get_mpz_t(), 2*w );
			mpz_clrbit ( svZr.get_mpz_t(), 2*w );

			tc->addExpectedOutput("Zi", svZi);
			tc->addExpectedOutput("Zr", svZr);
		}else{
			mpz_class big1 = (mpz_class(1) << (w));
			mpz_class big1P = (mpz_class(1) << (w-1));
			mpz_class big2 = (mpz_class(1) << (w));
			mpz_class big2P = (mpz_class(1) << (w-1));

			if ( svXi >= big1P)
				svXi = svXi - big1;
			if ( svXr >= big1P)
				svXr = svXi - big1;

			if ( svYi >= big2P)
				svYi = svYi - big2;
			if ( svYr >= big2P)
				svYr = svYr - big2;
			
			mpz_class svXrYr = svXr*svYr;
			mpz_class svXiYi = svXi*svYi;
			mpz_class svXrYi = svXr*svYi;
			mpz_class svXiYr = svXi*svYr;
			
			if ( svXrYr < 0){
				mpz_class tmpSUB = (mpz_class(1) << (2*w));
				svXrYr = tmpSUB + svXrYr; 
			}
			if ( svXiYi < 0){
				mpz_class tmpSUB = (mpz_class(1) << (2*w));
				svXiYi = tmpSUB + svXiYi; 
			}
			if ( svXrYi < 0){
				mpz_class tmpSUB = (mpz_class(1) << (2*w));
				svXrYi = tmpSUB + svXrYi; 
			}
			if ( svXiYr < 0){
				mpz_class tmpSUB = (mpz_class(1) << (2*w));
				svXiYr = tmpSUB + svXiYr; 
			}
			
			mpz_class svZi = svXrYi + svXiYr;
			mpz_class svZr = svXrYr - svXiYi;
			 
			// Don't allow overflow
			mpz_clrbit ( svZi.get_mpz_t(), 2*w );
			mpz_clrbit ( svZr.get_mpz_t(), 2*w );
			
			tc->addExpectedOutput("Zi", svZi);
			tc->addExpectedOutput("Zr", svZr);
		}
	}
	
	
	mpz_class IntTwiddleMultiplier::getTwiddleConstant(int constantType){
		mpfr_t twiddleExp, twiddleIm, twiddleRe, constPi, temp;
		mpz_class intTemp;
		
		mpfr_init2(twiddleIm, 	10*w);
		mpfr_init2(twiddleRe, 	10*w);
		mpfr_init2(twiddleExp, 	10*w);
		mpfr_init2(constPi, 	10*w);
		mpfr_init2(temp, 		10*w);
		
		mpfr_const_pi(	constPi, GMP_RNDN);
		
		mpfr_set_d(		twiddleExp, twiddleExponent, 			GMP_RNDN);
		mpfr_mul_2si(	twiddleExp, twiddleExp, 	1, 			GMP_RNDN);
		mpfr_mul(		twiddleExp, twiddleExp, 	constPi, 	GMP_RNDN);
		mpfr_div_d(		twiddleExp, twiddleExp, 	n, 			GMP_RNDN);
		
		mpfr_sin_cos(	twiddleIm, 	twiddleRe, 		twiddleExp, GMP_RNDN);
		
		switch(constantType){
			case TWIDDLERE:
							mpfr_set(temp, twiddleRe, GMP_RNDN);
							break;
			case TWIDDLEIM:
							mpfr_set(temp, twiddleIm, GMP_RNDN);
							break;
			case TWIDDLEIMADDRE:
							mpfr_set(temp, twiddleIm, GMP_RNDN);
							mpfr_add(temp, temp, twiddleRe, GMP_RNDN);
							break;
			case TWIDDLEIMSUBRE:
							mpfr_set(temp, twiddleIm, GMP_RNDN);
							mpfr_sub(temp, temp, twiddleRe, GMP_RNDN);
							break;
		}
		
		mpfr_mul_2si(temp, temp, wF, GMP_RNDN);
		mpfr_get_z(intTemp.get_mpz_t(), temp,  GMP_RNDN);
		
		mpfr_free_cache();
		mpfr_clears (twiddleExp, twiddleIm, twiddleRe, constPi, NULL);
		
		return intTemp;
	}
	
	std::string IntTwiddleMultiplier::getTwiddleConstantString(int constantType){
		std::ostringstream result;
		mpz_class temp;
		
		temp = getTwiddleConstant(constantType);
		if(temp<0)
			temp *= (-1);
		result << temp;
		
		return result.str();
	}
	
	int IntTwiddleMultiplier::getGCD(int x, int y){
		int temp;
		
		if(x==0)
			return y;
		else if(y==0)
			return x;
		else if(x<y){
			temp = x;
			x = y;
			y = temp;
		}
		
		return getGCD(y, x % y);
	}

}




















