//
// Created by Viktor Schmidt.
//

#include "ComparatorTable.hpp"

flopoco::ComparatorTable::ComparatorTable(OperatorPtr parentOp, flopoco::Target *target, int inputWidth_, int outputWidth_, std::vector<mpz_class> values_)
        : Table(parentOp, target, values_, inputWidth_, outputWidth_), inputWidth(inputWidth_), outputWidth(outputWidth_), values(values_)  {


    srcFileName="ComparatorTable";

    // definition of the name of the operator
    ostringstream name;
    name << "ComparatorTable" << inputWidth;
    setName(name.str());
    // Copyright
    //setCopyrightString("Viktor Schmidt 2018");

    REPORT(DEBUG,"Created ComparatorTable with " << inputWidth << " input bits and " << outputWidth << " outputBits.");
}

mpz_class flopoco::ComparatorTable::function(int x) {
    return values.at(x);
}
