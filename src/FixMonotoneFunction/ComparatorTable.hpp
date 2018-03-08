//
// Created by Viktor Schmidt.
//
#include "Table.hpp"
#include "utils.hpp"
#include <string>

#ifndef FLOPOCO_COMPARATORTABLE_HPP
#define FLOPOCO_COMPARATORTABLE_HPP


namespace flopoco {
    class ComparatorTable : public Table {
        int inputWidth;
        int outputWidth;
        std::vector<mpz_class> values;

    public:
        ComparatorTable(OperatorPtr parentOp, Target *target, int inputWidth_, int outputWidth_, std::vector<mpz_class> values_);

        mpz_class function(int x);
    };
}


#endif //FLOPOCO_COMPARATORTABLE_HPP
