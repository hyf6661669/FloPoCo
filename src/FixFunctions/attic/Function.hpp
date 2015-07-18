#ifndef _FUNCTION_HH_
#define _FUNCTION_HH_

#include <string>
#include <iostream>

#include "HOTBM/Util.hh"

#include <sollya.h>

//using namespace std;


namespace flopoco{

	class Function {
	public:
		Function(std::string name_, double xmin = 0, double xmax = 1, double scale = 1);
		virtual ~Function();

		std::string getName() const;
		double eval(double x) const;
		void eval(mpfr_t r, mpfr_t x) const;
		sollya_obj_t getSollyaNode() const;

	private:
		std::string name;
		sollya_obj_t node;
	
	}
		;
}
#endif // _FUNCTION_HH_
