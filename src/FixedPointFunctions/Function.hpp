#ifndef _FUNCTION_HH_
#define _FUNCTION_HH_

#include <string>
#include <iostream>

#include "HOTBM/Util.hh"
#include "../sollya.h"	

using namespace std;


namespace flopoco{

	class Function {
	public:
		Function(string name_, double xmin = 0, double xmax = 1, double scale = 1);
		virtual ~Function();

		string getName() const;
		double eval(double x) const;
		void eval(mpfr_t y, mpfr_t x) const;
		sollya_node_t getSollyaNode() const;
	
		/** Find x such that eval(x)==y, where a<x<b
		**/
		void eval_inverse(mpfr_t x, mpfr_t y, mpfr_t a, mpfr_t b) const;
	private:
		string name;
		sollya_node_t node;
		sollya_node_t diff;
	
	};
}
#endif // _FUNCTION_HH_
