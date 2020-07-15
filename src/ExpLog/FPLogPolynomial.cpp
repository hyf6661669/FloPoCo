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

	string signed_const(int val, int sz)
	{
		ostringstream o;
		o << '"';
		for (int i = 1 << (sz - 1); i > 0; i>>=1) o << ((val&i)!=0 ? '1' : '0');
		o << '"';
		return o.str();
	}

	FPLogPolynomial::FPLogPolynomial(OperatorPtr parentOp, Target* target,
		int wE, int wF, int inTableSize, int maxDegree):
		FPLog(parentOp, target, wE, wF)
	{
		setCopyrightString("F. de Dinechin, Q. Corradi  (2008-2020)");

		ostringstream o;
		srcFileName = "FPLogPolynomial";

		o << "FPLogPolynomial_" << wE << "_" << wF << "_" << inTableSize << "_";
		if (getTarget()->isPipelined()) o << getTarget()->frequencyMHz();
		else o << "comb";
		setNameWithFreqAndUID(o.str());

		addFPInput("X", wE, wF);
		addFPOutput("R", wE, wF, 2); // 2 because faithfully rounded

		if (inTableSize==0) inTableSize = 6;
		if (inTableSize < 4) throw string("FPLogPolynomial error: tables need at least 4 input bits");

		// see ???
		addComment("Centered mantissa");
		// more complex center condition?
		// getTarget()->lutInputs();
		vhdl << tab << declare("CenterCond") << " <= X"<<of(wF-1) << ";" << endl
			// << tab << declare("CenterCond") << " <= X > threshold;" << endl
			<< tab << declare("Y", wF+2) // u(0, -wF-1)
			<< " <= ('1' & X"<<range(wF-1, 0) << " & '0') when CenterCond = '0' else (\"01\" & X"<<range(wF-1, 0) << ");" << endl;

		addComment("Compute log(2^exponent)");
		addComment("Exponent is E = Ex - E0 + CenterCond = EX + -((1 << wE) - 1) + CenterCond");
		vhdl << tab << declare("E", wE+1) // s(wE, 0)
			<< " <= X"<<range(wF+wE-1, wF) << " + "<<signed_const(1 - (1 << wE), wE+1)<<" + CenterCond;" << endl;
		// E needs one more bit compared to Ex because of one value out of range, it would be unnecessary using IEEE floats

		addComment("Compute log(2^exponent) as E*log(2)");
		// lsb is -wF-3 because only 2 guard bits are needed and 1/2 < log(2) < 1
		newInstance("FixRealKCM", "MulLog2",
			"signedIn=0 msbIn=" + to_string(wE)
			+ " lsbIn=0"
			+ " lsbOut=" + to_string(-wF-3)
			+ " constant=log(2)",
			"X=>E", "R=>Elog2"); // Elog2 is s(wE, -wF-3)

		addComment("Range reduction: Table lookup for log and reciprocal, indexed by A");
		vhdl << tab << declare("A", inTableSize) << " <= X"<<range(wF-1, wF-inTableSize) << ";" << endl;
		const int wL = wF+inTableSize+3;
		{ // Reciprocal & log table
			const int threshold = 1 << (inTableSize - 1);
			const int max_val = 1 << inTableSize;
			vector<mpz_class> it0Content(max_val);
			vector<mpz_class> lt0Content(max_val);
			mpz_class t;
			mpfr_t r, l;
			mpfr_init2(r, inTableSize+2); // will contain the reciprocal
			mpfr_init2(l, wL); // will contain the log of that reciprocal
			// A = 0, reciprocal is 1
			it0Content[0] = max_val << 1;
			lt0Content[0] = 0;
			for (int A = 1; A < threshold; ++A) { // below theshold, no centering
				// reciprocal of 1/(1 + (A + 1/2)2**-inTableSize)
				mpfr_set_d(r, double(max_val << 1)/(((max_val | A) << 1) | 1), MPFR_RNDN);
				mpfr_log(l, r, MPFR_RNDN);
				mpfr_mul_2si(r, r, inTableSize+1, MPFR_RNDN);
				mpfr_mul_2si(l, l, wL, MPFR_RNDN);
				mpfr_get_z(it0Content[A].get_mpz_t(), r, MPFR_RNDN);
				mpfr_get_z(t.get_mpz_t(), l, MPFR_RNDN);
				lt0Content[A] = signedToBitVector(t, wL);
			}
			for (int A = threshold; A < max_val - 1; ++A) { // above threshold, centering occured
				// reciprocal of 2/(1 + (A + 1/2)2**-inTableSize)
				mpfr_set_d(r, (double)(max_val << 2)/(((max_val | A) << 1) | 1), MPFR_RNDN);
				mpfr_log(l, r, MPFR_RNDN);
				mpfr_mul_2si(r, r, inTableSize+1, MPFR_RNDN);
				mpfr_mul_2si(l, l, wL, MPFR_RNDN);
				mpfr_get_z(it0Content[A].get_mpz_t(), r, MPFR_RNDN);
				mpfr_get_z(t.get_mpz_t(), l, MPFR_RNDN);
				lt0Content[A] = signedToBitVector(t, wL);
			}
			// A = -1
			it0Content[max_val-1] = max_val << 1;
			lt0Content[max_val-1] = 0;
			// finalize
			mpfr_clear(r);
			mpfr_clear(l);
			Table::newUniqueInstance(this, "A", "Yhat", it0Content, "RecY0Table", inTableSize, inTableSize+2);
			Table::newUniqueInstance(this, "A", "L1", lt0Content, "LogY0Table", inTableSize, wL); // L1 is s(0, -wL+1)
		}
		newInstance("IntMultiplier", "Mul0",
			"wX=" + to_string(wF+2)
			+ " wY=" + to_string(inTableSize+2),
			"X=>Y, Y=>Yhat", "R=>YYhat"); // exact computation atm, YYhat is s(0, -inTableSize-wF-2)
		vhdl << tab << declare("Z", wF+3) // s(0, -wF-2) but there is a exponent of -inTableSize associated
			<< " <= YYhat"<<range(wF+2, 0) << ";" << endl;
		addComment("Get sign, it will be discarded by the LZOC so it has to be reinroduced afterwards");
		vhdl << tab << declare("SZ") << " <= YYhat"<<of(wF+2) << ";" << endl;

		addComment("Final exponent prediction when E=L1=0");
		addComment("Count how many singBit in the fraction");
		vhdl << tab << declare("Zns", wF+2) << " <= Z"<<range(wF+1, 0) << ";" << endl;
		newInstance("LZOCShifterSticky", "ExponentPredictorShifter",
			"wIn=" + to_string(wF+2)
			+ " wOut=" + to_string(wF+2)
			+ " countType=-1"
			+ " wCount=" + to_string(intlog2(wF-inTableSize)),
			"I=>Zns, OZb=>SZ", "Count=>Zshift, O=>mZ"); // Zshift is u(intlog2(wF-inTableSize)-1, 0)
		// (SZ & mZ) is s(0, -wF-2) but there is an exponent of -inTableSize-Zshift associated
		// in the case there is more than wF-inTableSize sign bits it means that not(E=L1=0), so the result is irrelevant anyway

		addComment("Last step, hardware for E=L1=0 is merged with the other case");
		vhdl << tab << declare("EL1")
			<< " <= '1' when E = "<<zg(wE+1)<<" and (A = "<<zg(inTableSize)<<" or A = "<<og(inTableSize)<<") else '0';" << endl
			<< tab << declare("Zu", wF+3) << " <= (SZ & mZ) when EL1 = '1' else Z;" << endl; // s(0, -wF-2)
		constexpr int g = 3;
		if (maxDegree >= 0) { // more range reduction
			addComment("Piecewise polynomial approximation");
			newInstance("FixFunctionByPiecewisePoly", "P",
				"f=(log1p(x*1b" + to_string(-inTableSize-1) + ")/(x*1b" + to_string(-inTableSize-1) + ")-1)"
				+ " signedIn=1"
				+ " lsbIn=" + to_string(-wF-2)
				+ " lsbOut=" + to_string(-wF-g)
				+ " d=" + to_string(maxDegree), // should depend on wF
				"X=>Z", "Y=>PZ"); // PZ is s(-inTableSize-1, -inTableSize-wF-g)
		} else { // unlimited degree
			addComment("Simple polynomial approximation");
			newInstance("FixFunctionBySimplePoly", "P",
				"f=(log1p(x*1b"+to_string(-inTableSize-1)+")/(x*1b"+to_string(-inTableSize-1)+")-1)*1b"+to_string(inTableSize)
				+ " signedIn=1"
				+ " lsbIn=" + to_string(-wF-2)
				+ " lsbOut=" + to_string(-wF-g),
				"X=>Z", "Y=>PZ"); // PZ is s(-inTableSize-1, -inTableSize-wF-g)
		}
		if (true || getTarget()->plainVHDL()) {
			newInstance("IntMultiplier", "Mul1",
				"wX=" + to_string(wF+g)
				+ " wY=" + to_string(wF+3)
				+ " wOut=" + to_string(wF+2-inTableSize),
				"X=>PZ, Y=>Zu", "R=>PZZu"); // PZZu is s(-inTableSize-1, -wF-2)
			// We need guard bits on Zu for cancellation, and PZZu has to be aligned
			vhdl << tab << declare("PZumZu", wF+3)
				<< " <= "<<rangeAssign(wF+2, wF+2-inTableSize, "PZZu"+of(wF+1-inTableSize))<<" & PZZu;" << endl;
			newInstance("IntAdder", "PolyAdd",
				"wIn=" + to_string(wF+3),
				"X=>Zu, Y=>PZumZu", "R=>PZu", "Cin=>'0'"); // PZu is s(0, -wF-2)

			addComment("Fixed point computation then floating point conversion when not(E=L1=0)");
			// Adjust bit pos and sign extend
			vhdl << tab << declare("L0", wE+wL) << " <= Elog2 & "<<zg(inTableSize-1)<<";" << endl // s(wE, -wL+1)
				<< tab << declare("L2", wL) << " <= "<<zg(inTableSize)<<" & PZu;" << endl; // s(0, -wL+1)
			newInstance("IntAdder", "FixedAdd1",
				"wIn=" + to_string(wL),
				"X=>L1, Y=>L2", "R=>LPlog", "Cin=>'0'"); // LPlog is s(0, -wL+1)
			vhdl << tab << declare("L12", wE+wL) << " <= "<<rangeAssign(wE+wL-1, wL, "LPlog"+of(wL-1))<<" & LPlog;" << endl;
			newInstance("IntAdder", "FixedAdd2",
				"wIn=" + to_string(wF+wE+inTableSize+g),
				"X=>L12, Y=>L0", "R=>FixL", "Cin=>'0'"); // FixL is s(wE,-wL+1)
		} else {
			throw "use without plainVHDL needs FixMultAdd to be revived";
			// newInstance("FixMultAdd", , "X=>PZ, Y=>Zu A=>Zu", "R=>L1");
			addComment("Fixed point computation then floating point conversion when not(E=L1=0)");
			// BitHeap bh(this, msb, lsb);
			// bh.addSignal("L0", shift_from_lsb);
			// bh.addSignal("L1", shift_from_lsb);
			// bh.addSignal("L2", shift_from_lsb);
			// bh.startCompression();
		}

		addComment("Almost final exponent");
		const int Esize = max(intlog2(wE+inTableSize), intlog2(wF-inTableSize));
		// Sign of FixL is always the right sign, so it is the final sign
		vhdl << tab << declare("Sr") << " <= FixL"<<of(wE+wL-1) << ";" << endl
		// We also need to correct the value of Ez by the table input size and 1 more if there is a cancellation on PZu
			<< tab << declare("PZuLT1") << " <= '1' when Sr = FixL"<<of(wE+wL-2) << " else '0';" << endl
			<< tab << declare("Ez", Esize) << " <= "<<signed_const(-inTableSize, Esize)
			<<" - ("<<zg(intlog2(wE+inTableSize)-intlog2(wF-inTableSize))<<" & Zshift) - PZuLT1;" << endl;

		// Exponent for the path not(E=L0=0), bounded by wE+inTableSize
		newInstance("LZOCShifterSticky", "Normalizer",
			"wIn=" + to_string(wE+wL)
			+ " wOut=" + to_string(wF+3)
			+ " countType=-1"
			+ " wCount=" + to_string(intlog2(wE+inTableSize)),
			"I=>FixL, OZb=>Sr", "Count=>Efix, O=>Mr");

		addComment("Result selection");
		vhdl << tab << declare("Er", Esize)
			<< " <= Ez when EL1 = '1' else ("<<zg(intlog2(wF-inTableSize)-intlog2(wE+inTableSize))<<" & Efix);" << endl
			<< tab << declare("Fr", wF+3) << " <= PZu when EL1 = '1' else Mr;" << endl;

		addComment("Final rounding with wE=" + to_string(wE) + " and wF=" + to_string(wF));
		vhdl << tab << declare("xFr", wF+3) << " <= "<<rangeAssign(wF+2, 0, "Sr")<<" xor Fr;" << endl
			<< tab << declare("Lfp", wE+wF+3) << " <= (Er + "<<signed_const((1 << wE) - 1, wE)<<") & xFr;" << endl
			<< tab << declare("halfulp", wE+wF+3) << " <= "<<zg(wE+wF)<<" & \"100\";" << endl;
		newInstance("IntAdder", "FinalRounding",
			"wIn=" + to_string(wE+wF+3),
			"X=>Lfp, Y=>halfulp, Cin=>Sr", "R=>REF");

		addComment("Conversion to FloPoCo floating point format");
		vhdl << tab << "with X"<<range(wE+wF+2, wE+wF) << " select " << declare("RES", 3) << " <=\n"
			<< tab << tab << "\"101\" when \"000\",\n"
			<< tab << tab << "\"101\" when \"001\",\n"
			<< tab << tab << "\"01\" & Sr when \"010\",\n"
			<< tab << tab << "\"110\" when \"011\",\n"
			<< tab << tab << "\"100\" when \"100\",\n"
			<< tab << tab << "\"110\" when \"101\",\n"
			<< tab << tab << "\"110\" when \"110\",\n"
			<< tab << tab << "\"110\" when \"111\",\n"
			<< tab << tab << "\"---\" when others;" << endl
			<< tab << "R <= RES & REF"<<range(wE+wF+2, 3) << ";" << endl;
	}

	FPLogPolynomial::~FPLogPolynomial()
	{
	}

}
