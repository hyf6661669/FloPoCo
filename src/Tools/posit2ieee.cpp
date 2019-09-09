/*
 * Utility for creating the conversion table of a given posit format in a given
 * IEEE format
 *
 * Author : Luc Forget
 *
 * This file is part of the FloPoCo project developed by the Arenaire
 * team at Ecole Normale Superieure de Lyon
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
*/

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <gmpxx.h>
#include <cstdio>
#include <mpfr.h>
#include <cstdlib>

#include"../utils.hpp"
#include "../TestBenches/PositNumber.hpp"
#include "../TestBenches/IEEENumber.hpp"

using namespace std;
using namespace flopoco;

static void usage(char *name) 
{
	cerr << endl << "Usage: " << name << " N WES WE WF" << endl;
	cerr << "N WES : format parameters of the posit number source format" << endl;
	cerr << "WE WF : format parameters of the destination IEEE format" << endl;

	exit(EXIT_FAILURE);
}

int check_positive(char* s, char* cmd, bool strict=true) {
  int n=atoi(s);
  string strstr = (strict) ? "strict" : "";
  if (n<0 || (n==0 and strict)){
    cerr<<"ERROR: got "<<s<<", expected " << strstr  << " positive number."<<endl;
    usage(cmd);
  }
  return n;
}

int main(int argc, char* argv[] )
{
	if(argc != 5) usage(argv[0]);
	int N = check_positive(argv[1], argv[0]);
	int WES = check_positive(argv[2], argv[0], false);
	int WE = check_positive(argv[3], argv[0]);
	int WF = check_positive(argv[4], argv[0]);

	mpz_class val;
	mpz_class limit = mpz_class{1} << N;

	mpfr_t posit_mpfr;
	mpfr_init2(posit_mpfr, N);

	mpfr_t ieee_mpfr;
	mpfr_init2(ieee_mpfr, WF+1);

	ofstream output_file;

	stringstream name;
	name << "posit_" << N << "_" << WES << "_to_ieee_" << WE << "_" << WF << "_conversion_table.txt";
	output_file.open(name.str());

	if (! output_file.is_open()) {
		cerr << "Error while opening file " << name.str() << " for writing output" << endl;
		exit(EXIT_FAILURE);
	}

	char* mpfr_str_repr = 0;

	output_file << "posit_" << N << "_" << WES << 
		", float_" << WE << "_" << WF << "_val, " << 
		"float_" << WE << "_" << WF << "_bin" << endl ;

	for (val = 0 ; val < limit ; val += 1) {
		output_file << "0x" << val.get_str(16) << ", ";

		auto pos_val = PositNumber{N, WES, val};
		pos_val.getMPFR(posit_mpfr);

		mpfr_set(ieee_mpfr, posit_mpfr, MPFR_RNDN);

		mpfr_asprintf(&mpfr_str_repr, "%Rf", ieee_mpfr);
		string ieee_val_repr{mpfr_str_repr};
		output_file << ieee_val_repr << ", ";		

		mpfr_free_str(mpfr_str_repr);
		
		output_file << ieee2bin(ieee_mpfr, WE, WF) << endl;
	}

	return 0;
}


