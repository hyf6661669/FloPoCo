#ifndef GENERICTABLE_HPP
#define GENERICTABLE_HPP

#include <vector>

#include "../Table.hpp"

namespace flopoco{

	class GenericTable : public Table
	{
	public:
		GenericTable(Target* target, int wIn, int wOut, std::vector<mpz_class> values,  std::map<std::string, double> inputDelays = emptyDelayMap);

		~GenericTable();

		mpz_class function(int x);

		int wIn;
		int wOut;
		std::vector<mpz_class> values; //the values to be stored in the table
	};

}
#endif //GENERICTABLE_HPP

