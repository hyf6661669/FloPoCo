#include <iostream>
#include <sstream>

#include "gmp.h"
#include "mpfr.h"

#include "CordicSinCosClassic.hpp"

using namespace std;

namespace flopoco{

	CordicSinCosClassic::CordicSinCosClassic(Target* target, int w_, map<string, double> inputDelays) 
		: Operator(target), w(w_)
	{
		int wcs = 1+w, wx = 1+w;
		
		//TODO: verify the validity of the necessary guard bits
		guard = ceil(log2(wcs));
		
		mpfr_init2(scale, 10*w);
		mpfr_set_d(scale, -1.0, GMP_RNDN);
		mpfr_mul_2si(scale, scale, -w, GMP_RNDN);
		mpfr_add_d(scale, scale, 1.0, GMP_RNDN);
		
		srcFileName="CordicSinCosClassic";
		ostringstream name;

		setCopyrightString ( "Istoan Matei, Florent de Dinechin (2008-2012)" );
		if(target->isPipelined())
			name << "CordicSinCosClassic_" << 1+w <<"_f"<< target->frequencyMHz() << "_uid" << getNewUId();
		else
			name << "CordicSinCosClassic_" << 1+w << "_uid" << getNewUId();
		setName( name.str() );
		
		//RE
		initResourceEstimation();
		//RE
		
		//FLP
		initFloorplanning(0.45);
		//FLP

		// declaring inputs
		addInput  ( "X"  , wcs, true );

		// declaring output
		addOutput  ( "C"  , wcs, 2 );
		addOutput  ( "S"  , wcs, 2 );
		
		setCriticalPath(getMaxInputDelays(inputDelays));
		
		manageCriticalPath(target->localWireDelay(wcs + guard) + target->adderDelay(1+w) + target->lutDelay());
		
		//reduce the argument X to [0, 1/2)
		vhdl << tab << declare("absX", wx ) << "<= (X xor (" << wx-1 << " downto 0 => X(" << wx-1 << ")))"
											<< " + " 
											<< "(" << zg(wx-1, 0) << " & X(" << wx-1 << "));"<< endl;
		vhdl << tab << declare("reducedX", wx + guard) 
						<< "<= (absX(" << wx-1 << " downto " << wx-2 << ") - (\'0\' & (absX(" << wx-1 << ") xor absX(" << wx-2 << "))))" 
						<< " & absX(" << wx-3 << " downto 0) & " << zg(guard, 0) << ";" << endl;
		vhdl << tab << declare("quadrantX", 2) << " <= X(" << wx-1 << " downto " << wx-2 << ");" << endl;
		
		syncCycleFromSignal("absX");
		
		//RE
		resourceEstimate << addArithOp(wx, 2, 1);
		resourceEstimate << addAdderSubtracter(wx, wx, 0, 1);
		resourceEstimate << addAdderSubtracter(2, 2, 0, 1);
		resourceEstimate << addArithOp(1, 2, 1);
		//RE
		
		//create the C0, S0, X0 and D0 signals for the first stage
		//compute the scale factor
		mpfr_t kfactor, temp;
		
		mpfr_init2(kfactor, 10*wcs);
		mpfr_init2(temp, 10*wcs);
		mpfr_set_d(kfactor, 1, GMP_RNDN);
		for(int i=0; i<wcs-1; i++){
			mpfr_set_d(temp, 1, GMP_RNDN);
			mpfr_mul_2si(temp, temp, -2*i, GMP_RNDN);
			mpfr_add_d(temp, temp, 1.0, GMP_RNDN);
			
			mpfr_mul(kfactor, kfactor, temp, GMP_RNDN);
		}
		mpfr_sqrt(kfactor, kfactor, GMP_RNDN);
		mpfr_d_div(kfactor, 1.0, kfactor, GMP_RNDN);

		mpfr_mul(kfactor, kfactor, scale, GMP_RNDN);
		
		
		manageCriticalPath(target->localWireDelay(wcs + guard) + target->lutDelay());
		
		vhdl << tab << declare("C0", wcs + guard) << "<= " << zg(1, 0) << " & \"" << generateFixPointNumber(kfactor, 0, wcs-1+guard) << "\";" << endl;
		vhdl << tab << declare("S0", wcs + guard) << "<= " << zg(1, 0) << " & " << zg(wcs-1+guard) << ";" << endl;
		vhdl << tab << declare("X0", wx  + guard) << "<= reducedX;" << endl;
		vhdl << tab << declare("D0") << "<= reducedX(" << wx+guard-1 << ");" << endl;
		
		//RE
		resourceEstimate << addReg(wcs+guard-1);
		resourceEstimate << addReg(wcs+guard-1);
		//RE
				
		//create the stages of micro-rotations
		int stage;
		
		wcs += guard;
		wx += guard;
		
		//build the cordic stages
		for(stage=0; stage<wcs-1; stage++){
			CordicIteration* cordicIteration = new CordicIteration(target, wx, wcs, stage);
			oplist.push_back(cordicIteration);
			
			inPortMap(cordicIteration, "Xin", join("X", stage));
			inPortMap(cordicIteration, "Cin", join("C", stage));
			inPortMap(cordicIteration, "Sin", join("S", stage));
			inPortMap(cordicIteration, "Din", join("D", stage));
			outPortMap(cordicIteration, "Xout", join("X", stage+1));
			outPortMap(cordicIteration, "Cout", join("C", stage+1));
			outPortMap(cordicIteration, "Sout", join("S", stage+1));
			outPortMap(cordicIteration, "Dout", join("D", stage+1));
			vhdl << instance(cordicIteration, join("cordicIteration", stage)) << endl;

			syncCycleFromSignal(join("X", stage+1));
			
			//decrement the size of Z
			wx--;
		}
		
		manageCriticalPath(target->localWireDelay(1+w+1) + target->adderDelay(1+w+1));
		syncCycleFromSignal(join("D", stage));
		
		vhdl << tab << declare("preRoundedIntCout", 1+w+1) << "<= " << join("C", stage) << "(" << wcs-1 << " downto " << guard-1 << ");" << endl;
		vhdl << tab << declare("preRoundedIntSout", 1+w+1) << "<= " << join("S", stage) << "(" << wcs-1 << " downto " << guard-1 << ");" << endl;
		
		vhdl << tab << declare("roundedIntCout", 1+w+1) << "<= preRoundedIntCout "
															<< "+"
															<< " (" << zg(1+w, 0) << " & \'1\');" << endl;
		vhdl << tab << declare("roundedIntSout", 1+w+1) << "<= preRoundedIntSout "
															<< "+"
															<< " (" << zg(1+w, 0) << " & \'1\');" << endl;
															
		//RE
		resourceEstimate << addAdderSubtracter(w+2, w+2, 0, 2);
		//RE
															
		//assign output
		manageCriticalPath(target->localWireDelay(1+w) + target->adderDelay(1+w));
		
		vhdl << tab << declare("reducedC", 1+w) << "<= roundedIntCout(" << 1+w << " downto 1);" << endl;
		vhdl << tab << declare("reducedS", 1+w) << "<= roundedIntSout(" << 1+w << " downto 1);" << endl;
		
		vhdl << tab << declare("negReducedC", 1+w) << "<= (reducedC xor (" << w << " downto 0 => \'1\'))"
												   << " + " 
												   << "(" << zg(w, 0) << " & \'1\');"<< endl;
		vhdl << tab << declare("negReducedS", 1+w) << "<= (reducedS xor (" << w << " downto 0 => \'1\'))"
												   << " + " 
												   << "(" << zg(w, 0) << " & \'1\');"<< endl;
												   
		//RE
		resourceEstimate << addAdderSubtracter(w+1, w+1, 0, 2);
		//RE
												   
		manageCriticalPath(target->localWireDelay(1+w) + target->lutDelay());
		
		vhdl << tab << "with quadrantX select" << endl;
		vhdl << tab << tab << "C <= reducedC when \"00\"," << endl;
		vhdl << tab << tab << tab << " negReducedS when \"01\"," << endl;
		vhdl << tab << tab << tab << " negReducedS when \"10\"," << endl;
		vhdl << tab << tab << tab << " reducedC when \"11\"," << endl;
		vhdl << tab << tab << tab << " reducedC when others;" << endl;	//edit to signal error
		
		vhdl << tab << "with quadrantX select" << endl;
		vhdl << tab << tab << "S <= reducedS when \"00\"," << endl;
		vhdl << tab << tab << tab << " reducedC when \"01\"," << endl;
		vhdl << tab << tab << tab << " negReducedC when \"10\"," << endl;
		vhdl << tab << tab << tab << " negReducedS when \"11\"," << endl;
		vhdl << tab << tab << tab << " reducedS when others;" << endl;	//edit to signal error
		
		//RE
		resourceEstimate << addMux(w+1, 4, 2);
		
		resourceEstimate << addPipelineFF();
		resourceEstimate << addWireCount();
		resourceEstimate << addPortCount();
		resourceEstimate << addComponentResourceCount();
		//RE
	};
	
	
	
	void CordicSinCosClassic::emulate(TestCase * tc) 
	{
		/* Get I/O values */
		mpz_class svZ = tc->getInputValue("X");
		mpfr_t z, constPi, rsin, rcos;
		mpz_t rsin_z, rcos_z;
		
		/* Compute correct value */
		mpfr_init2(z, 10*w);
		
		mpfr_init2(constPi, 10*w);
		
		
		mpfr_set_z (z, svZ.get_mpz_t(), GMP_RNDN); // this rounding is exact
		mpfr_div_2si (z, z, w, GMP_RNDN); // this rounding is acually exact
		
		mpfr_const_pi( constPi, GMP_RNDN);
		mpfr_mul(z, z, constPi, GMP_RNDN);
		
		mpfr_init2(rsin, 10*w); 
		mpfr_init2(rcos, 10*w); 


		mpz_init2 (rsin_z, 2*w);
		mpz_init2 (rcos_z, 2*w);
		mpfr_sin(rsin, z, GMP_RNDN); 
		mpfr_cos(rcos, z, GMP_RNDN);
		mpfr_mul(rsin, rsin, scale, GMP_RNDN);
		mpfr_mul(rcos, rcos, scale, GMP_RNDN);

		mpfr_add_d(rsin, rsin, 6.0, GMP_RNDN);
		mpfr_add_d(rcos, rcos, 6.0, GMP_RNDN);
		mpfr_mul_2si (rsin, rsin, w, GMP_RNDN); // exact rnd here
		mpfr_mul_2si (rcos, rcos, w, GMP_RNDN); // exact rnd here

		// Rounding down
		mpfr_get_z (rsin_z, rsin, GMP_RNDD); // there can be a real rounding here
		mpfr_get_z (rcos_z, rcos, GMP_RNDD); // there can be a real rounding here
		mpz_class sin_zc (rsin_z), cos_zc (rcos_z);
		sin_zc -= mpz_class(6)<<w;
		cos_zc -= mpz_class(6)<<w;

		tc->addExpectedOutput ("S", sin_zc);
		tc->addExpectedOutput ("C", cos_zc);
		
		// Rounding up
		mpfr_get_z (rsin_z, rsin, GMP_RNDU); // there can be a real rounding here
		mpfr_get_z (rcos_z, rcos, GMP_RNDU); // there can be a real rounding here
		mpz_class sin_zcu (rsin_z), cos_zcu (rcos_z);
		sin_zcu -= mpz_class(6)<<w;
		cos_zcu -= mpz_class(6)<<w;

		tc->addExpectedOutput ("S", sin_zcu);
		tc->addExpectedOutput ("C", cos_zcu);
		


		// clean up
		mpfr_clears (z, rsin, rcos, NULL);		
		mpfr_free_cache();
	}


	void CordicSinCosClassic::buildStandardTestCases(TestCaseList * tcl) 
	{
		TestCase* tc;
		mpf_t zinit;
		mpfr_t z;
		mpz_t z_z;
		
		//mpf_set_default_prec (1+wI+wF+guard);
		
		mpfr_init2(z, 1+w+ceil(log2(1 + w)));
		mpz_init2 (z_z, 1+w+ceil(log2(1 + w)));
		
		//z=0
		tc = new TestCase (this);
		tc -> addInput ("X",mpz_class(0));
		emulate(tc);
		tcl->add(tc);
		
		//z=pi/2
		tc = new TestCase (this);
		
		mpf_init2   (zinit, 1+w+ceil(log2(1 + w)));
		//mpf_set_str (zinit, "1.5707963267949e0", 10);
		mpf_set_str (zinit, "0.5e0", 10);
		mpfr_set_f (z, zinit, GMP_RNDD); 
		
		mpfr_mul_2si (z, z, w, GMP_RNDD); 
		mpfr_get_z (z_z, z, GMP_RNDD);  
		tc -> addInput ("X",mpz_class(z_z));
		emulate(tc);
		tcl->add(tc);
		
		//z=pi/6
		tc = new TestCase (this); 
		
		mpf_init2   (zinit, 1+w+ceil(log2(1 + w)));
		//mpf_set_str (zinit, "0.5235987755983e0", 10);
		mpf_set_str (zinit, "0.16666666666666e0", 10);
		mpfr_set_f (z, zinit, GMP_RNDD); 
		
		mpfr_mul_2si (z, z, w, GMP_RNDD); 
		mpfr_get_z (z_z, z, GMP_RNDD);  
		tc -> addInput ("X",mpz_class(z_z));
		emulate(tc);
		tcl->add(tc);
		
		//z=pi/4
		tc = new TestCase (this);
		
		mpf_init2   (zinit, 1+w+ceil(log2(1 + w)));
		//mpf_set_str (zinit, "0.78539816339745e0", 10);
		mpf_set_str (zinit, "0.25e0", 10);
		mpfr_set_f (z, zinit, GMP_RNDD); 
		
		mpfr_mul_2si (z, z, w, GMP_RNDD); 
		mpfr_get_z (z_z, z, GMP_RNDD);  
		tc -> addInput ("X",mpz_class(z_z));
		emulate(tc);
		tcl->add(tc);
		
		//z=pi/3
		tc = new TestCase (this);
		
		mpf_init2   (zinit, 1+w+ceil(log2(1 + w)));
		//mpf_set_str (zinit, "1.0471975511966e0", 10);
		mpf_set_str (zinit, "0.33333333333333e0", 10);
		mpfr_set_f (z, zinit, GMP_RNDD);
		
		mpfr_mul_2si (z, z, w, GMP_RNDD); 
		mpfr_get_z (z_z, z, GMP_RNDD);  
		
		tc -> addInput ("X",mpz_class(z_z));
		emulate(tc);
		tcl->add(tc);
		
		mpfr_clears (z, NULL);
	}



	std::string CordicSinCosClassic::generateFixPointNumber(mpfr_t x, int wI, int wF)
	{
		std::string result;
		int size = wI+wF;
		mpz_class h;
		
		if(x<0){
			mpfr_neg (x, x,  GMP_RNDN);
		}
		
		mpfr_mul_2si(x, x, wF, GMP_RNDN);
		
		mpfr_get_z(h.get_mpz_t(), x,  GMP_RNDN);         
		result = unsignedBinary(h, size);
		
		return result;
	}

}











