#ifndef BaseMultiplier_HPP
#define BaseMultiplier_HPP

#include <string>
#include <iostream>
#include <string>
#include "Target.hpp"

namespace flopoco {


    class BaseMultiplier {

	public:
        BaseMultiplier();
        ~BaseMultiplier();

		
		
    private:

	

        string srcFileName; //for debug outputs

        string uniqueName_; /**< useful only to enable same kind of reporting as for FloPoCo operators. */
	
	};
}
#endif
