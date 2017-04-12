#ifndef BaseMultiplier_HPP
#define BaseMultiplier_HPP

#include <string>
#include <iostream>
#include <string>
#include "Target.hpp"
#include "Operator.hpp"

namespace flopoco {


    class BaseMultiplierLUT : public BaseMultiplier
    {

	public:
        BaseMultiplierLUT(bool isSignedX, bool isSignedY, int wX, int wY);

        /**
         * @brief generateOperator generates an instance of the corresponding Operator that realizes the given shape
         * @return the generated operator
         */
        virtual Operator *generateOperator();

        /**
         * @brief Returns true if x and y coordinates are at valid shape positions
         * @param x: x-coordinate (relative to multiplier coordinate)
         * @param y: y-coordinate (relative to multiplier coordinate)
         * @return true if x and y coordinates are at valid shape positions, false otherwise
         */
        virtual bool shapeValid(int x, int y);

    private:

	};

    class BaseMultiplierLUTOp : public Operator
    {
        BaseMultiplierLUTOp();
    };
}
#endif
