#ifdef HAVE_PAGLIB
// general c++ library for manipulating streams
#include <iostream>
#include <sstream>
#include <fstream>


#include "gmp.h"
#include "mpfr.h"

#include "FullyParallelFFT.hpp"

using namespace std;
namespace flopoco {




FullyParallelFFT::FullyParallelFFT(Target* target, int wIn_, string rotatorFileName_, string FFTRealizationFileName_)
    : Operator(target),
      wIn(wIn_),
      rotatorFileName(rotatorFileName_),
      FFTRealizationFileName(FFTRealizationFileName_)
{

    std::ifstream rotFile(rotatorFileName);
    std::ifstream FFTFile(FFTRealizationFileName);

    // Parse in rotators
    std::string line;
    vector<string> rotators;
    while (std::getline(rotFile, line))
    {
        std::istringstream newLine(line);
        int numRot;
        string rotRealization;
        newLine >> numRot >> rotRealization;
        if ((unsigned int) (numRot+1 )> rotators.size())
            rotators.resize(numRot+1,"");
        rotators[numRot]=rotRealization;
    }

    if (DEBUG<=UserInterface::verbose){
        REPORT(DEBUG,"rotators which were read in");
        for (unsigned n=0; n<rotators.size(); n++) {
            REPORT(DEBUG,"rotator "<< n <<" is implemented by "<<rotators.at( n ));
        }
    }

    // Parse in FFT realization
    string FFTRealizationString;
    vector<vector<int64_t> > FFTRealization;
    std::getline(FFTFile, FFTRealizationString);

    size_t pos=0;
    string tmp=FFTRealizationString;
    vector<int64_t> FFTRow;
    while(tmp.size()>0)
    {
        pos=tmp.find_first_of(",; ");

        if(pos != string::npos)
        {
            switch (tmp[pos]) {
            case ',':
                FFTRow.push_back( atol(tmp.substr(0,pos).c_str()) );
                break;
            case ' ':
            case ';':
                FFTRow.push_back( atol(tmp.substr(0,pos).c_str()) );
                FFTRealization.push_back(FFTRow);
                FFTRow.clear();
                break;
            default:
                break;
            }
            tmp = tmp.substr(pos+1);
        }
        else
        {
            FFTRow.push_back( atol(tmp.c_str()) );
            FFTRealization.push_back(FFTRow);
            tmp = "";
        }
    }

    if (DEBUG<=UserInterface::verbose){
        REPORT(DEBUG,"FFT implementation which was read in");
        for (unsigned r=0; r<FFTRealization.size(); r++) {
            REPORT(DEBUG,"Row " << r);
            for (unsigned c=0; c<FFTRealization[r].size(); c++) {
                REPORT(DEBUG,FFTRealization[r].at(c) << " ");
            }
        }
    }



    unsigned int N = FFTRealization.size();
    unsigned int n = FFTRealization[0].size()+1;
    srcFileName="FullyParallelFFT";
    ostringstream name;
    name << "FullyParallelFFT"<< N << "p" << wIn << "bit";

    useNumericStd();
    setName(name.str());
    // Copyright
    setCopyrightString("Konrad Möller, Martin Kumm, Mario Garrido");


    for (unsigned r=0; r<N; r++) {
        // declaring inputs
        {stringstream tmpSignal; tmpSignal << "In_"<< r << "_x";
            addInput(tmpSignal.str() , wIn);}
        {stringstream tmpSignal; tmpSignal << "In_"<< r << "_y";
            addInput(tmpSignal.str() , wIn);}

        //declare intermediate signals
        for (unsigned c=0; c<n-1; c++) {

            {stringstream tmpSignal; tmpSignal << "post_s"<<c+1<<"_r"<<r<<"_x";
                declare(tmpSignal.str(),wIn);}
            {stringstream tmpSignal; tmpSignal << "post_s"<<c+1<<"_r"<<r<<"_y";
                declare(tmpSignal.str(),wIn);}
        }

        // declaring outputs
        {stringstream tmpSignal; tmpSignal << "Out_"<< r << "_x";
            addOutput(tmpSignal.str() , wIn);}
        {stringstream tmpSignal; tmpSignal << "Out_"<< r << "_y";
            addOutput(tmpSignal.str() , wIn);}
    }

    //Building the FFTadders
    for (unsigned stage=1; stage<n+1; stage++) {
        for (unsigned row=0; row<N; row++) {
            if (stage==1){ //first stage requires the input
                stringstream tmpSignalx; tmpSignalx << "pre_s"<<stage<<"_r"<<row<<"_x";
                stringstream tmpSignaly; tmpSignaly << "pre_s"<<stage<<"_r"<<row<<"_y";
                if(row<N/2){
                    vhdl << tab << declare(tmpSignalx.str(),wIn+1) << " <= " << "std_logic_vector(signed(resize(signed(In_"<<row<<"_x)+signed(In_"<<row+N/(2*stage)<<"_x),"<<wIn+1<<")));"<<endl;
                    vhdl << tab << declare(tmpSignaly.str(),wIn+1) << " <= " << "std_logic_vector(signed(resize(signed(In_"<<row<<"_y)+signed(In_"<<row+N/(2*stage)<<"_y),"<<wIn+1<<")));"<<endl;
                }else{
                    vhdl << tab << declare(tmpSignalx.str(),wIn+1) << " <= " << "std_logic_vector(signed(resize(signed(In_"<<row-N/(2*stage)<<"_x)-signed(In_"<<row<<"_x),"<<wIn+1<<")));"<<endl;
                    vhdl << tab << declare(tmpSignaly.str(),wIn+1) << " <= " << "std_logic_vector(signed(resize(signed(In_"<<row-N/(2*stage)<<"_y)-signed(In_"<<row<<"_y),"<<wIn+1<<")));"<<endl;
                }
            }
            else if (stage<=n) {
                stringstream tmpSignalx; tmpSignalx << "pre_s"<<stage<<"_r"<<row<<"_x";
                stringstream tmpSignaly; tmpSignaly << "pre_s"<<stage<<"_r"<<row<<"_y";
                if(row<N/2){
                    vhdl << tab << declare(tmpSignalx.str(),wIn+1) << " <= " << "std_logic_vector(signed(resize(signed(post_s"<<stage-1<<"_r"<<row<<"_x)+signed(post_s"<<stage-1<<"_r"<<row+N/(2*stage)<<"_x),"<<wIn+1<<")));"<<endl;
                    vhdl << tab << declare(tmpSignaly.str(),wIn+1) << " <= " << "std_logic_vector(signed(resize(signed(post_s"<<stage-1<<"_r"<<row<<"_y)+signed(post_s"<<stage-1<<"_r"<<row+N/(2*stage)<<"_y),"<<wIn+1<<")));"<<endl;
                }else{
                    vhdl << tab << declare(tmpSignalx.str(),wIn+1) << " <= " << "std_logic_vector(signed(resize(signed(post_s"<<stage-1<<"_r"<<row-N/(2*stage)<<"_x)-signed(post_s"<<stage-1<<"_r"<<row<<"_x),"<<wIn+1<<")));"<<endl;
                    vhdl << tab << declare(tmpSignaly.str(),wIn+1) << " <= " << "std_logic_vector(signed(resize(signed(post_s"<<stage-1<<"_r"<<row-N/(2*stage)<<"_y)-signed(post_s"<<stage-1<<"_r"<<row<<"_y),"<<wIn+1<<")));"<<endl;
                }

            }
            else{ //last stage feeds outputs and no rotation
                stringstream tmpSignalx; tmpSignalx << "pre_s"<<stage<<"_r"<<row<<"_x";
                stringstream tmpSignaly; tmpSignaly << "pre_s"<<stage<<"_r"<<row<<"_y";
                if(row<N/2){
                    vhdl << tab << declare(tmpSignalx.str(),wIn+1) << " <= " << "std_logic_vector(signed(resize(signed(post_s"<<stage-1<<"_r"<<row<<"_x)+signed(post_s"<<stage-1<<"_r"<<row+N/(2*stage)<<"_x),"<<wIn+1<<")));"<<endl;
                    vhdl << tab << declare(tmpSignaly.str(),wIn+1) << " <= " << "std_logic_vector(signed(resize(signed(post_s"<<stage-1<<"_r"<<row<<"_y)+signed(post_s"<<stage-1<<"_r"<<row+N/(2*stage)<<"_y),"<<wIn+1<<")));"<<endl;
                }else{
                    vhdl << tab << declare(tmpSignalx.str(),wIn+1) << " <= " << "std_logic_vector(signed(resize(signed(post_s"<<stage-1<<"_r"<<row-N/(2*stage)<<"_x)-signed(post_s"<<stage-1<<"_r"<<row<<"_x),"<<wIn+1<<")));"<<endl;
                    vhdl << tab << declare(tmpSignaly.str(),wIn+1) << " <= " << "std_logic_vector(signed(resize(signed(post_s"<<stage-1<<"_r"<<row-N/(2*stage)<<"_y)-signed(post_s"<<stage-1<<"_r"<<row<<"_y),"<<wIn+1<<")));"<<endl;
                }
            }

        }
    }

}


void FullyParallelFFT::emulate(TestCase * tc) {
    /* This function will be used when the TestBench command is used in the command line
           we have to provide a complete and correct emulation of the operator, in order to compare correct output generated by this function with the test input generated by the vhdl code */
}


void FullyParallelFFT::buildStandardTestCases(TestCaseList * tcl) {
    // please fill me with regression tests or corner case tests!
}


OperatorPtr FullyParallelFFT::parseArguments(Target *target, vector<string> &args) {
    int wIn;
    string rotatorFileName, FFTRealizationFileName;
    UserInterface::parseInt(args, "wIn", &wIn); // param0 has a default value, this method will recover it if it doesnt't find it in args,
    UserInterface::parseString(args, "rotatorFileName", &rotatorFileName);
    UserInterface::parseString(args, "FFTRealizationFileName", &FFTRealizationFileName);
    return new FullyParallelFFT(target, wIn, rotatorFileName, FFTRealizationFileName);
}

void FullyParallelFFT::registerFactory(){
    UserInterface::add("FullyParallelFFT", // name
                       "Generator for a fully parallel multiplier-less FFT implementation using shift and add", // description, string
                       "Filters and FFTs", // category, from the list defined in UserInterface.cpp
                       "", //seeAlso
                       // Now comes the parameter description string.
                       // Respect its syntax because it will be used to generate the parser and the docs
                       // Syntax is: a semicolon-separated list of parameterDescription;
                       // where parameterDescription is parameterName (parameterType)[=defaultValue]: parameterDescriptionString
                       "wIn(int)=16: A first parameter, here used as the input size; \
                       rotatorFileName(string): A second parameter, here used as the output size; \
            FFTRealizationFileName(string): A third parameter, here used as the output size",
                                            // More documentation for the HTML pages. If you want to link to your blog, it is here.
                                            "Feel free to experiment with its code, it will not break anything in FloPoCo. <br> Also see the developper manual in the doc/ directory of FloPoCo.",
                                            FullyParallelFFT::parseArguments
                                            ) ;
}

}//namespace
#endif // HAVE_PAGLIB
