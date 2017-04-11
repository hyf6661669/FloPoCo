#ifndef MultiplierSolutionParser_HPP
#define MultiplierSolutionParser_HPP

#include <string>
#include <iostream>
#include <fstream>
#include "Operator.hpp" //for THROWERRROR, etc
#include "utils.hpp"
#include "BitHeap/BitHeap.hpp"
#include <list>
#include <utility>

namespace flopoco {


    class MultiplierSolutionParser {

	public:
        MultiplierSolutionParser(std::string fileName);
        ~MultiplierSolutionParser();

        list<pair< unsigned int, pair<unsigned int, unsigned int> > > getSolution();
        bool readSolution();
		
		list<pair< unsigned int, pair<unsigned int, unsigned int> > > solution; 
		//solution.first: type
		//solution.second.first: x-coordinate
		//solution.second.second: y-coordinate
		
    private:

		bool variableIsTrue(string line);
		void addVariableToSolution(string line);

        string solFileName;

        string srcFileName; //for debug outputs

        string uniqueName_; /**< useful only to enable same kind of reporting as for FloPoCo operators. */
	
	};
}
#endif
