#ifdef HAVE_PAGLIB

#include "Operator.hpp"
#include "utils.hpp"

namespace flopoco {

	// new operator class declaration
	class FullyParallelFFT : public Operator {
	public:
        int wIn;
        string rotator_file_name,FFT_realization_file_name;


	public:

        FullyParallelFFT(Target* target, int wIn_ = 16,  string rotator_file_name_=0, string FFT_realization_file_name_=0);

		// destructor
		~FullyParallelFFT() {};
		void emulate(TestCase * tc);
		void buildStandardTestCases(TestCaseList* tcl);

		/** Factory method that parses arguments and calls the constructor */
		static OperatorPtr parseArguments(Target *target , vector<string> &args);
		
		/** Factory register method */ 
		static void registerFactory();

	};

}//namespace
#endif // HAVE_PAGLIB
