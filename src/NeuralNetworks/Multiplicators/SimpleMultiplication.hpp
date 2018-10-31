//
// Created by nfiege on 28/06/18.
//

#ifndef FLOPOCO_SIMPLEMULTIPLICATION_HPP
#define FLOPOCO_SIMPLEMULTIPLICATION_HPP

#include "Operator.hpp"

#include "utils.hpp"
#include <string>

namespace flopoco {
    class SimpleMultiplication : public Operator {

    public:
        SimpleMultiplication(Target *target, int wIn1, int wIn2);
    };
}

#endif //FLOPOCO_SIMPLEMULTIPLICATION_HPP
