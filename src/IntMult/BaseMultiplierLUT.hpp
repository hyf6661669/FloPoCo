#ifndef BaseMultiplierLUT_HPP
#define BaseMultiplierLUT_HPP

#include <string>
#include <iostream>
#include <string>
#include <gmp.h>
#include <gmpxx.h>
#include "Target.hpp"
#include "Operator.hpp"
#include "Table.hpp"
#include "BaseMultiplier.hpp"

namespace flopoco {


    class BaseMultiplierLUT : public BaseMultiplier
    {

	public:
        BaseMultiplierLUT(bool isSignedX, bool isSignedY, int wX, int wY);

        /**
         * @brief generateOperator generates an instance of the corresponding Operator that realizes the given shape
         * @return the generated operator
         */
        virtual Operator *generateOperator(Target *target);

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
    public:
        BaseMultiplierLUTOp(Target* target, bool isSignedX, bool isSignedY, int wX, int wY);
    };

    class BaseMultiplierLUTTable : public Table
    {
    public:
        BaseMultiplierLUTTable(Target* target, int dx, int dy, int wO, bool negate=false, bool signedX=false, bool signedY=false );
        mpz_class function(int x);
    protected:
        int dx, dy, wO;
        bool negate, signedX, signedY;
    };
}
#endif
