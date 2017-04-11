#ifndef BaseMultiplierCollection_HPP
#define BaseMultiplierCollection_HPP

#include <string>
#include <iostream>
#include <string>
#include "Target.hpp"

namespace flopoco {


    class BaseMultiplierCollection {

	public:
        BaseMultiplierCollection(Target *target);
        ~BaseMultiplierCollection();

		//BaseMultiplier* getBaseMultiplier(int shape)
		
		
    private:

	
		Target* target_;

        string srcFileName; //for debug outputs

        string uniqueName_; /**< useful only to enable same kind of reporting as for FloPoCo operators. */
	
	};
}
#endif
