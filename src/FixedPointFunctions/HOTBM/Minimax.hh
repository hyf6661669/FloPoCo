#ifndef _MINIMAX_HH_
#define _MINIMAX_HH_

#include "../Function.hpp"
#include "MPPolynomial.hh"

using namespace std;


namespace flopoco{

	class Minimax {
	public:
		Minimax(Function &f, double ia, double ib, int d);
		~Minimax();

		MPPolynomial &getMPP() const;
		void getMPErr(mpfr_t mpErr_) const;
	
		//! Get the number of points used in the infinite norm calc
		static int getInfNormPoints();
	
		//! Set the number of points used in the infinite norm calc, and return the previous value
		static int setInfNormPoints(int pp);

	private:
		MPPolynomial *mpP;
		mpfr_t mpErr;
	
		static int infNormPoints;
	};
}
#endif // _MINIMAX_HH_
