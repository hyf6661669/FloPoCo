#ifndef flopoco_random_utils_make_table_hpp
#define flopoco_random_utils_make_table_hpp

#include "Table.hpp"
#include "gmp.h"
#include "mpfr.h"
#include "FPNumber.hpp"

namespace flopoco
{
namespace random
{
	
inline Operator *MakeSinglePortTable(Target *target, std::string name, int wElts, const std::vector<mpz_class> &contents, bool hardRam=false, map<string, double> inputDelays = emptyDelayMap )
{
	class SinglePortTable
		: public Table
	{
	private:
		std::vector<mpz_class> m_elements;
	public:
		SinglePortTable(Target* target, std::string name, int wIn, int wOut, const std::vector<mpz_class> &elements, bool hardRam, map<string, double> inputDelays = emptyDelayMap )
			: Table(target, wIn, wOut, /*minIn*/ 0, /*maxIn*/elements.size()-1, /*logicTable*/hardRam?0:1,  inputDelays)
			, m_elements(elements)
		{
			setName(name);
			assert((1u<<wIn) >=elements.size());
		}
		
		virtual mpz_class function(int x)
		{
			return m_elements.at(x);
		}
	};
	
	int wIn=(int)ceil(log(contents.size())/log(2.0));
	
	mpz_class mm=mpz_class(1)<<wElts;
	for(unsigned i=0;i<contents.size();i++){
		if(contents[i] < 0)
			throw std::string("MakeSinglePortTable - Currently elements must be positive (need to be encoded for signed values).");
		if(contents[i] >=mm)
			throw std::string("MakeSinglePortTable - Element is not in range [0,2^wElts).");
	}
	
	return new SinglePortTable(target, name, wIn, wElts, contents, hardRam, inputDelays);
}

}; // random
}; // flopoco

#endif
