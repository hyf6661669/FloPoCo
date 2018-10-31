//
// Created by nfiege on 12/07/18.
//

#ifndef FLOPOCO_FINDINDEXOFMAXIMUM_HPP
#define FLOPOCO_FINDINDEXOFMAXIMUM_HPP

#include "Operator.hpp"

#include "utils.hpp"
#include <string>
#include <vector>

namespace flopoco {

    class FindIndexOfMaximum : public Operator {

    public:
        FindIndexOfMaximum(Target* target, unsigned int wIn, unsigned int numberOfInputs);

        ~FindIndexOfMaximum() {}
        unsigned int getIndexWordSize() const { return this->indexWordSize; }

        /** Factory method that parses arguments and calls the constructor */
        static OperatorPtr parseArguments(Target *target , vector<string> &args);

        /** Factory register method */
        static void registerFactory();

    private:
        unsigned int indexWordSize;
        vector <string> inputNames;
        vector <string> indexNames;
    };
} // namespace flopoco

#endif //FLOPOCO_FINDINDEXOFMAXIMUM_HPP