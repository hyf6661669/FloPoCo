#ifndef _PIECEWISEFUNCTION_HH_
#define _PIECEWISEFUNCTION_HH_

#include <string>
#include <iostream>
#include <vector>
#include "HOTBM/Util.hh"
#include "Function.hpp"
#include <stdlib.h>
//using namespace std;


namespace flopoco{

	class PiecewiseFunction {
	public:
		PiecewiseFunction(std::string name_);
		virtual ~PiecewiseFunction();

		std::string getName() const;
		std::vector<Function*>getPiecewiseFunctionArray() const;
    Function* getPiecewiseFunctionArray(int i) const;
	private:
		std::string name;
		std::vector<Function*> fArray;
	}
		;
}
#endif // _PIECEWISEFUNCTION_HH_
