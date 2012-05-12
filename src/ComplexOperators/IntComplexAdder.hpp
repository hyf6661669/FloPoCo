#ifndef IntComplexAdder_HPP
#define IntComplexAdder_HPP
#include <vector>
#include <sstream>

#include "../Operator.hpp"
#include "../IntAdder.hpp"

namespace flopoco{
	
	/**
	 * Complex adder for fixed point numbers
	 */
	class IntComplexAdder : public Operator
	{
	public:
		IntComplexAdder(Target* target, int wI, int wF, bool signedOperator = true, map<string, double> inputDelays = emptyDelayMap);
		~IntComplexAdder();
		
		void emulate(TestCase * tc);

		int wI, wF, w;
		bool signedOperator;

	};
}
#endif
