#ifndef GenericLut_H
#define GenericLut_H

#include "Operator.hpp"
#include "utils.hpp"

#include <vector>
#include <string>
#include "bool_eq.hpp"

namespace flopoco {

	// new operator class declaration
    class GenericLut : public Operator {
        std::vector<bool_eq> equations_;
        unsigned int wIn_;
        unsigned int wOut_;
      public:
        GenericLut( Target *target, const std::string &name, const std::vector<bool_eq> &equations );

        GenericLut( Target *target, const std::string &name, const std::map<unsigned int, unsigned int> &groups, const unsigned int &wIn = 0, const unsigned int &wOut = 0 );

        void build();

        void build_select();

		// destructor
        ~GenericLut() {};

	};


}//namespace

#endif
