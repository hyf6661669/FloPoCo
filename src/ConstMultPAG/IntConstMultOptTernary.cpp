/*
  Integer constant multiplication using minimum number of adders due to

  Gustafsson, O., Dempster, A., Johansson, K., Macleod, M., & Wanhammar, L. (2006).
  Simplified Design of Constant Coefficient Multipliers. Circuits, Systems, and Signal Processing, 25(2), 225–251.

  All constants up to 19 bit will be realized optimal using precomputed tables provided by the SPIRAL project (http://spiral.ece.cmu.edu/mcm/).
  
  Author : Martin Kumm kumm@uni-kassel.de, (emulate adapted from Florent de Dinechin, Florent.de.Dinechin@ens-lyon.fr)

  All rights reserved.

*/

#include <iostream>
#include <sstream>
#include <vector>
#include <gmp.h>
#include <mpfr.h>
#include <gmpxx.h>
#include "../utils.hpp"
#include "../Operator.hpp"

#include "IntConstMultOptTernary.hpp"
#include "../ConstMultPAG/tscm_solutions.hpp"
#include "../ConstMultPAG/pagexponents.hpp"

#include "rpag/types.h"
#include "rpag/compute_successor_set.h"
#include "rpag/log2_64.h"
#include "fundamental_extended.hpp"
//#include "../ConstMultPAG/paglib/compute_successor_set.hpp"
//#include "../ConstMultPAG/paglib/log2_64.hpp"
//#include "../ConstMultPAG/paglib/fundamental.hpp"

#include <algorithm>

namespace rpag {} //required for clang to know the namespace

using namespace std;
using namespace rpag;

namespace flopoco{

    IntConstMultOptTernary::IntConstMultOptTernary(Target* target, int wIn, int coeff, bool syncInOut) : ConstMultPAG(target, wIn, "", false, syncInOut, 1000, false)
    {
		int maxCoefficient = 4194303; //=2^22-1

        int shift;
        int coeffOdd = fundamental(coeff, &shift);

		if(coeffOdd <= maxCoefficient)
        {
            cout << "nofs[" << coeff << "]=" << nofs[(coeff-1)>>1][0] << " " << nofs[(coeff-1)>>1][1] << endl;

            int nof1 = nofs[(coeffOdd-1)>>1][0];
            int nof2 = nofs[(coeffOdd-1)>>1][1];

//            cout << "coeff=" << coeff << ", nof1=" << nof1 << ", nof2=" << nof2 << endl;

            set<int> coefficient_set;
            coefficient_set.insert(coeff);

            set<int> nof_set;
            if(nof1 != 0) nof_set.insert(nof1);
            if(nof2 != 0) nof_set.insert(nof2);

            string adderGraphString;

            computeAdderGraphTernary(nof_set, coefficient_set, adderGraphString);

            stringstream adderGraphStringStream;
            adderGraphStringStream << "{" << adderGraphString;

            if(coeff != coeffOdd)
            {
                if(adderGraphString.length() > 1) adderGraphStringStream << ",";
				adderGraphStringStream << "{" << coeff << ",2," << coeffOdd << ",1," << shift << ",0,0,0}";
            }
            adderGraphStringStream << "}";

            cout << "adder_graph=" << adderGraphStringStream.str() << endl;

            ProcessConstMultPAG(target,adderGraphStringStream.str());

            ostringstream name;
            name << "IntConstMultOptTernary_" << coeffOdd << "_" << wIn;
            setName(name.str());

        }
        else
        {
			cerr << "Error: Coefficient too large, max. odd coefficient is " << maxCoefficient << endl;
            exit(-1);
        }

	}

    OperatorPtr flopoco::IntConstMultOptTernary::parseArguments( Target *target, vector<string> &args ) {
        int wIn, constant;

        UserInterface::parseInt( args, "wIn", &wIn );
        UserInterface::parseInt( args, "constant", &constant );

        return new IntConstMultOptTernary( target, wIn, constant, false);
    }

    void flopoco::IntConstMultOptTernary::registerFactory() {

        UserInterface::add( "IntConstMultOptTernary", // name
                            "Integer constant multiplication using shift and ternary additions in an optimal way (i.e., with minimum number of ternary adders). Works for coefficients up to 4194303 (22 bit)", // description, string
                            "ConstMultDiv", // category, from the list defined in UserInterface.cpp
                            "", //seeAlso
                            "wIn(int): Input word size; \
                            constant(int): constant;",
                            "Nope.",
                            IntConstMultOptTernary::parseArguments
                          ) ;
    }

}
