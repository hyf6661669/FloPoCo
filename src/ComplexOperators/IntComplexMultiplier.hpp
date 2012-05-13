#ifndef IntComplexMultiplier_HPP
#define IntComplexMultiplier_HPP
#include <vector>
#include <sstream>

#include "../Operator.hpp"
#include "../IntMultiplier.hpp"
#include "../IntAdder.hpp"

namespace flopoco{
	
	/**
	 * Complex multiplier for complex numbers
	 * Depending on the value of hasLessMultiplications, the multiplication
	 * (a+jb)*(c+jd) can be either
	 * 		Re(z)=a*c+-b*d
	 * 		Im(z)=a*d+b*c, with 4 multiplications and 3 additions
	 * or
	 * 		m1=(a+b)*c
	 * 		m2=(d+c)*b
	 * 		m3=(d-c)*a
	 * 		Re(z)=m1-m2
	 * 		Im(z)=m1+m3, with 3 multiplications and 5 additions
	 */
	class IntComplexMultiplier : public Operator
	{
	public:
		IntComplexMultiplier(Target* target, int wI, int wF, bool signedOperator = true, bool hasLessMultiplications = false);
		~IntComplexMultiplier();
		
		void emulate(TestCase * tc);

		int wI, wF, w;
		bool signedOperator;

	};
}
#endif
