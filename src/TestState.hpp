#ifndef TESTSTATE_HPP
#define TESTSTATE_HPP

#include <ostream>
#include <sstream>
#include <vector>
#include <string>
#include <gmpxx.h>

//using namespace std;
namespace flopoco {
	class TestState
	{
	public :
		TestState ( std::string param );
	
		/**
		 * Test the equality between two TestState
		 **/
		bool equality (  TestState * ts );

		std::vector < int > vectInt;
		std::vector < float > vectFloat;
		std::vector < std::string > vectString;
		std::vector < mpz_class > vectMpz;
		std::vector < bool > vectBool;

		std::string paramTypes; /** String representing type of parameters of the chosen operator int the right order */
		std::string toString(); /** get a textual representation of the state, to be passed to the flopoco command line for instance*/
		int counter;   /**< counter to be aware of how many test we have done in this instance*/
	};
}
#endif
