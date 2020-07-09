/*
  An FP logarithm for FloPoCo

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon

  Author : Florent de Dinechin, Florent.de.Dinechin@ens-lyon.fr

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,
  2008-2010.
  All rights reserved.

*/

#include <sstream>
#include "FPLogPolynomial.hpp"
#include "Table.hpp"
#include "utils.hpp"


using namespace std;


namespace flopoco{

	FPLogPolynomial::FPLogPolynomial(OperatorPtr parentOp, Target* target,
		int wE, int wF, int inTableSize):
		FPLog(parentOp, target, wE, wF)
	{
		setCopyrightString("F. de Dinechin, C. Klein, Q. Corradi  (2008-2020)");

		ostringstream o;
		srcFileName = "FPLogPolynomial";

		o << "FPLogPolynomial_" << wE << "_" << wF << "_" << inTableSize << "_";
		if(getTarget()->isPipelined())
			o << getTarget()->frequencyMHz();
		else
			o << "comb";
		setNameWithFreqAndUID(o.str());

		addFPInput("X", wE, wF);
		addFPOutput("R", wE, wF, 2); // 2 because faithfully rounded

		if (inTableSize==0 || inTableSize < 4) {
			// find default value
			throw string("FPLogPolynomial error: tables need at least 4 input bits");
		}

		constexpr int g = 10; // guard bits
		// see ???
		addComment("Centered mantissa");
		// more complex center condition?
		// getTarget()->lutInputs();
		vhdl << tab << declare("CenterCond") << " <= X"<<of(wF-1) << ";" << endl
			// << tab << declare("CenterCond") << " <= X > threshold;" << endl
			<< tab << declareFixPoint("Y", false, 0, -wF-1) << " <= ('1' & X"<<range(wF-1, 0)
			<< " & '0') when CenterCond = '0' else (\"01\" & X"<<range(wF-1, 0) << ");" << endl;

		addComment("Compute log(2^exponent)");
		// addComment("Exponent is E = Ex - E0 + CenterCond = Ex - ((1 << (wE - 1)) - 1) + CenterCond");
		// vhdl << tab << declare("E", wE + 1) << " <= X"<<range(wF+wE-1, wF) << " - ('0' & " << og(wE-1) << ") + CenterCond;" << endl;
		addComment("Exponent is E = Ex - E0 + CenterCond = EX + ((1 << wE) + 1) + CenterCond");
		vhdl << tab << declare("E", wE + 1) << " <= X"<<range(wF+wE-1, wF) << " + ('1' & " << zg(wE-2) << " & '1') + CenterCond;" << endl;
		// E needs one more bit compared to Ex because of one value out of range, it would be unnecessary using IEEE floats

		addComment("Compute log(2^exponent) as E*log(2)");
		// lsb is -wF-3 because of 2 guard bits and 1/2 < log(2) < 1
		cout << join("signedIn=1 msbIn=",wE+1," lsbIn=0 lsbOut=",-wF-3," constant=\"log(2)\"") << endl;
		newInstance("FixRealKCM", "MulLog2",
			join("signedIn=1 msbIn=",wE+1," lsbIn=0 lsbOut=",-wF-3," constant=\"log(2)\""),
			"X=>E", "R=>Elog2");
		vhdl << tab << declareFixPoint("L0", true, wE, -wF-3) << " <= Elog2;" << endl;

		addComment("Range reduction: Table lookup for log and reciprocal, indexed by A");
		vhdl << tab << declare("A", inTableSize) << " <= X"<<range(wF-1, wF-inTableSize) << ";" << endl;
		const int wL = wF+inTableSize+g;
		{ // Reciprocal & log table
			const int threshold = 1 << (inTableSize - 1);
			const int max_val = 1 << inTableSize;
			vector<mpz_class> it0Content(max_val+1);
			vector<mpz_class> lt0Content(max_val+1);
			mpfr_t r, l;
			mpfr_init2(r, inTableSize+2); // will contain the reciprocal
			mpfr_init2(l, wL); // will contain the log of that reciprocal
			// A = 0, reciprocal is 1
			it0Content[0] = max_val << 1;
			lt0Content[0] = 0;
			for (int A = 1; A < threshold; ++A) { // below theshold, no centering
				mpfr_set_d(r, (max_val << 1)/(((max_val | A) << 1) | 1), MPFR_RNDN); // reciprocal of 1/(1 + (A + 1/2)2**-inTableSize)
				mpfr_log(l, r, MPFR_RNDN);
				mpfr_mul_2si(r, r, inTableSize+1, MPFR_RNDN);
				mpfr_mul_2si(l, l, wL, MPFR_RNDN);
				mpfr_get_z(it0Content[A].get_mpz_t(), r, MPFR_RNDN);
				mpfr_get_z(lt0Content[A].get_mpz_t(), l, MPFR_RNDN);
				cout << it0Content[A] << endl;
				cout << lt0Content[A] << endl;
			}
			for (int A = threshold; A < max_val - 1; ++A) { // above threshold, centering occured
				mpfr_set_d(r, max_val/(((max_val | A) << 1) | 1), MPFR_RNDN); // reciprocal of 1/(1 + (A + 1/2)2**-inTableSize)/2
				mpfr_log(l, r, MPFR_RNDN);
				mpfr_mul_2si(r, r, inTableSize+1, MPFR_RNDN);
				mpfr_mul_2si(l, l, wL, MPFR_RNDN);
				mpfr_get_z(it0Content[A].get_mpz_t(), r, MPFR_RNDN);
				mpfr_get_z(lt0Content[A].get_mpz_t(), l, MPFR_RNDN);
				cout << it0Content[A] << endl;
				cout << lt0Content[A] << endl;
			}
			// A = -1
			it0Content[max_val] = max_val << 1;
			lt0Content[max_val] = 0;
			// finalize
			mpfr_clear(r);
			mpfr_clear(l);
			Table::newUniqueInstance(this, "A", "Yhat", it0Content, "RecY0Table", inTableSize, inTableSize+2);
			Table::newUniqueInstance(this, "A", "LYhat", lt0Content, "LogY0Table", inTableSize, wL);
		}
		/*
		vhdl << tab << declareFixPoint("L1", true, 0, -wL) << " <= LYhat;" << endl;
		newInstance("IntMultiplier", "Mul0", join("wX=",wF+2," wY=",inTableSize+g), "X=>Y Y=>Yhat", "R=>YYhat");
		vhdl << tab << declareFixPoint("Z", true, -inTableSize, -inTableSize-wF-2) << " <= YYhat"<<range(wF+2, 0) << ";" << endl; // exact computation atm
		addComment("Get sign, it will be discarded by the LZOC so it has to be reinroduced afterwards");
		vhdl << tab << declare("signBit") << "YYhat"<<of(wF+2) << ";" << endl;

		addComment("Final exponent prediction when E=L1=0");
		addComment("Count how many singBit in the fraction");
		vhdl << tab << declare("Zns", wF+2) << " <= YYhat"<<range(wF+1, 0) << ";" << endl;
		newInstance("LZOCShifterSticky", "ExponentPredictorShifter",
			join("wIn=",wF+2," wOut=",wF+2," countType=-1 wCount=",intlog2(wF-inTableSize)),
			"I=>Zns OZB=>signBit", "Count=>Zshift O=>Zs");
		// in the case there is more than wF-inTableSize signBit it means that not(E=L1=0), so the result is irrelevant

		addComment("Last step, hardware for E=L1=0 is merged with the other case");
		vhdl << tab << declare("EL1") << " <= E = 0 and (A = 0 or A = -1);" << endl
			tab << declare("Zu", wF+3) << " <= (singBit & Zs) when EL1 else Z;" << endl;
		if (false) { // more range reduction: piecewise polynomial
			// ostringstream o;
			// o << "f=\"log1p(x/2^" << inTableSize << ")\" lsbIn=-?? lsbOut="-wF-g;
			// newInstance("FixFunctionByPiecewisePoly", "LastLog", , "X=>Zu", "Y=>PZu")
		} else { // simple polynomial
			addComment("Simple polynomial approximation");
			newInstance("FixFunctionBySimplePoly", "P",
				join("f=\"log1p(x/2^",inTableSize,")/x - 1\" signedIn=1 lsbIn=",-wF-2," lsbOut=",-wF-g),
				"X=>Z", "Y=>PZ");
		}
		if (true || plainVHDL) { // restart HERE
			newInstance("IntMultiplier", "Mul1", join("wX=",wF+g," wY=",wF+3," wOut=",wF+g), "X=>PZ Y=>Zu", "R=>PZZu");
			vhdl << tab << declare("PZum1", 2*wF+g+4) << " <= PZZu"<<of(2*wF+g+2) << " & PZZu;" << endl;
			resizeFixPoint();
			// TODO: Sign extend for overflow
			newInstance("IntAdder", "PolyAdd", join("wIn=",wF+g), "X=>Zu Y=>PZum1", "R=>PZu");
			// TODO: Check output bit size!
			// TODO: declarefixpoint PZu
			vhdl << tab << declareFixPoint("L2", true, -inTableSize, -inTableSize-wF-g) << " <= PZu;" << endl;
			addComment("Fixed point computation then floating point conversion when not(E=L1=0)");
			// TODO: adjust bit pos and sign extend
			newInstance("IntAdder", "FixedAdd1", join("wIn=",wF+wE+inTableSize+g), "X=>L0 Y=>L1", "R=>ELlog");
			// TODO: adjust bit pos and sign extend
			newInstance("IntAdder", "FixedAdd2", join("wIn=",wF+wE+inTableSize+g), "X=>ELlog Y=>PZu", "R=>FixL");
		} else {
			throw "use without plainVHDL needs FixMultAdd to be revived";
			// newInstance("FixMultAdd", , "X=>PZ Y=>Zu A=>Zu", "R=>L1");
			addComment("Fixed point computation then floating point conversion when not(E=L1=0)");
			// BitHeap bh(this, msb, lsb);
			// bh.addSignal("L0", shift_from_lsb);
			// bh.addSignal("L1", shift_from_lsb);
			// bh.addSignal("L2", shift_from_lsb);
			// bh.startCompression();
		}
		addComment("Almost final exponent");
		vhdl << tab << declare("Ez", intlog2(wF)) << " <= -" << inTableSize << " - Zshift;" << endl;
		// TODO: get signs

		// TODO: OZB=>signbit
		newInstance("LZOCShifterSticky", "Normalizer",
			join("wIn=",wF+wE+inTableSize+g," wOut=",wF+g," countType=-1 wCount=",intlog2(wF)+intlog2(wE)),
			"I=>FixL", "Count=>Er O=>Mr");

		addComment("Result selection");
		vhdl << tab << declare("Er") << " <= Ez when EL1 else Er;" << endl
			<< tab << declare("Fr") << " <= PZ when EL1 else Mr;" << endl
			<< tab << declare("Sr") << " <= ";//TODO: the sign

		addComment(join("Final rounding with wE=",wE," and wF=",wF));
		// Carry propagate using IntAdder?
		// also need to xor bits with sign bit
		// newInstance("IntAdder", "FinalRounding", join("wIn="), "X=> Cin=>", );
		addComment("Conversion to FloPoCo floating point format");
		*/
		vhdl << tab << "R <= 0;" << endl;
	}

	FPLogPolynomial::~FPLogPolynomial()
	{
	}

}
