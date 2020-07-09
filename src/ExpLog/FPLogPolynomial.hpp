#ifndef __FPLOGPOLYNOMIAL_HPP
#define __FPLOGPOLYNOMIAL_HPP

#include "FPLog.hpp"

namespace flopoco{

	class FPLogPolynomial: public FPLog
	{
	public:
		FPLogPolynomial(OperatorPtr parentOp, Target* target, int wE, int wF, int inTableSize=0);
		~FPLogPolynomial();
	};
}

#endif
