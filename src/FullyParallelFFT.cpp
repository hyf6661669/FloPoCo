#ifdef HAVE_PAGLIB
// general c++ library for manipulating streams
#include <iostream>
#include <sstream>
#include <fstream>


#include "gmp.h"
#include "mpfr.h"
#include "FloPoCo.hpp"
#include "ConstMultPAG/ConstMultPAG.hpp"
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
    vector<pair<int,int> > rotatorVal;
    while (std::getline(rotFile, line))
    {
        std::istringstream newLine(line);
        int numRot,real,imag;
        string rotRealization;
        newLine >> numRot >> real >> imag >> rotRealization;
        if ((unsigned int) (numRot+1 )> rotators.size()){
            rotators.resize(numRot+1,"");
            rotatorVal.resize(numRot+1,make_pair(0,0));
        }
        rotators[numRot]=rotRealization;
        rotatorVal[numRot]=make_pair(real,imag);
    }

    if (DEBUG<=UserInterface::verbose){
        REPORT(DEBUG,"rotators which were read in");
        for (unsigned n=0; n<rotators.size(); n++) {
            REPORT(DEBUG,"rotator "<< n << " " << rotatorVal[n].first << " " << rotatorVal[n].second <<" is implemented by "<<rotators.at( n ));
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

        /*  //declare intermediate signals
        for (unsigned c=0; c<n-1; c++) {

            {stringstream tmpSignal; tmpSignal << "post_s"<<c+1<<"_r"<<r<<"_x";
                declare(tmpSignal.str(),wIn);}
            {stringstream tmpSignal; tmpSignal << "post_s"<<c+1<<"_r"<<r<<"_y";
                declare(tmpSignal.str(),wIn);}
        }*/

        // declaring outputs
        {stringstream tmpSignal; tmpSignal << "Out_"<< r << "_x";
            addOutput(tmpSignal.str() , wIn);}
        {stringstream tmpSignal; tmpSignal << "Out_"<< r << "_y";
            addOutput(tmpSignal.str() , wIn);}
    }

    vector<vector<bool> > negX, negY;
    negX.resize(N);
    negY.resize(N);
    for (unsigned stage=1; stage<n+1; stage++) {
        for (unsigned row=0; row<N; row++) {
            negX[row].resize(n+1,false);
            negY[row].resize(n+1,false);
        }
    }

    //Building the FFT
    for (unsigned stage=1; stage<n+1; stage++) {
        for (unsigned row=0; row<N; row++) {
            bool swap=false;
            stringstream tmpSignalx,tmpSignaly;

            tmpSignalx << "pre_s"<<stage<<"_r"<<row<<"_x";
            tmpSignaly << "pre_s"<<stage<<"_r"<<row<<"_y";

            // get output signs and connect rotator
            if (stage<n) {
                // get rotator for current row
                string rotatorString = rotators.at(FFTRealization[row].at(stage-1));
                //put rotators here
                stringstream tmpOutSignalx; tmpOutSignalx << "post_s"<<stage<<"_r"<<row<<"_x";
                stringstream tmpOutSignaly; tmpOutSignaly << "post_s"<<stage<<"_r"<<row<<"_y";

                Operator* cMult = new flopoco::ConstMultPAG(target,wIn+1,rotatorString,false,false,1000,false);
                cMult->setName("rotator",to_string(FFTRealization[row].at(stage-1)));
                inPortMap(cMult,"x_in0",declare(tmpSignalx.str(),wIn+1));
                inPortMap(cMult,"x_in1",declare(tmpSignaly.str(),wIn+1));
                list<flopoco::ConstMultPAG::output_signal_info> tmp_output_list = ((flopoco::ConstMultPAG*)(cMult))->GetOutputList();

                list<flopoco::ConstMultPAG::output_signal_info>::iterator it=tmp_output_list.begin();

                if((*it).output_factors[0][0]==rotatorVal[FFTRealization[row].at(stage-1)].first && (*it).output_factors[0][1]==-rotatorVal[FFTRealization[row].at(stage-1)].second)
                {
                    swap=false; negX[row].at(stage)=false;
                }
                else if((*it).output_factors[0][0]==-rotatorVal[FFTRealization[row].at(stage-1)].first && (*it).output_factors[0][1]==rotatorVal[FFTRealization[row].at(stage-1)].second)
                {
                    swap=false; negX[row].at(stage)=true;
                }
                else if((*it).output_factors[0][0]==rotatorVal[FFTRealization[row].at(stage-1)].second && (*it).output_factors[0][1]==rotatorVal[FFTRealization[row].at(stage-1)].first)
                {
                    swap=true; negY[row].at(stage)=false;
                }
                else if((*it).output_factors[0][0]==-rotatorVal[FFTRealization[row].at(stage-1)].second && (*it).output_factors[0][1]==-rotatorVal[FFTRealization[row].at(stage-1)].first)
                {
                    swap=true; negY[row].at(stage)=true;
                }


                if (!swap)
                    outPortMap(cMult,(*it).signal_name,tmpOutSignalx.str(),true);
                else
                    outPortMap(cMult,(*it).signal_name,tmpOutSignaly.str(),true);

                it++;

                if((*it).output_factors[0][0]==rotatorVal[FFTRealization[row].at(stage-1)].first && (*it).output_factors[0][1]==-rotatorVal[FFTRealization[row].at(stage-1)].second)
                {
                    swap=true; negX[row].at(stage)=false;
                }
                else if((*it).output_factors[0][0]==-rotatorVal[FFTRealization[row].at(stage-1)].first && (*it).output_factors[0][1]==rotatorVal[FFTRealization[row].at(stage-1)].second)
                {
                    swap=true; negX[row].at(stage)=true;
                }
                else if((*it).output_factors[0][0]==rotatorVal[FFTRealization[row].at(stage-1)].second && (*it).output_factors[0][1]==rotatorVal[FFTRealization[row].at(stage-1)].first)
                {
                    swap=false; negY[row].at(stage)=false;
                }
                else if((*it).output_factors[0][0]==-rotatorVal[FFTRealization[row].at(stage-1)].second && (*it).output_factors[0][1]==-rotatorVal[FFTRealization[row].at(stage-1)].first)
                {
                    swap=false; negY[row].at(stage)=true;
                }

                if (swap){
                    outPortMap(cMult,(*it).signal_name,tmpOutSignalx.str(),true);
                }
                else
                    outPortMap(cMult,(*it).signal_name,tmpOutSignaly.str(),true);
                addSubComponent(cMult);
                stringstream rotatorName;
                rotatorName << "rotator_" << FFTRealization[row].at(stage-1) << "at_s" << stage << "_r" << row;
                cerr << rotatorName.str() << "\t swap " << swap << " negX " << negX[row].at(stage) <<  " negY " << negY[row].at(stage) << endl;
                vhdl << instance(cMult, rotatorName.str());
            }
//            //put butterflies here
//            if (stage==1){ //first stage requires the input
//                if(!((row) & (1<<((int)(log2(N)-stage))))){
//                    vhdl << tab << tmpSignalx.str() << " <= " << "std_logic_vector(resize(signed(In_"<<row<<"_x),"<<wIn+1<<")+resize(signed(In_"<<row+N/(pow(2,stage))<<"_x),"<<wIn+1<<"));"<<endl;
//                    vhdl << tab << tmpSignaly.str() << " <= " << "std_logic_vector(resize(signed(In_"<<row<<"_y),"<<wIn+1<<")+resize(signed(In_"<<row+N/(pow(2,stage))<<"_y),"<<wIn+1<<"));"<<endl;
//                }else{
//                    vhdl << tab << tmpSignalx.str() << " <= " << "std_logic_vector(resize(signed(In_"<<row-N/(pow(2,stage))<<"_x),"<<wIn+1<<")-resize(signed(In_"<<row<<"_x),"<<wIn+1<<"));"<<endl;
//                    vhdl << tab << tmpSignaly.str() << " <= " << "std_logic_vector(resize(signed(In_"<<row-N/(pow(2,stage))<<"_y),"<<wIn+1<<")-resize(signed(In_"<<row<<"_y),"<<wIn+1<<"));"<<endl;
//                }
//            }
//            else if (stage<n) { // intermediate stages

//                if(!((row) & (1<<((int)(log2(N)-stage))))){
//                    vhdl << tab << tmpSignalx.str() << " <= " << "std_logic_vector(resize("<<(negX[row].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row<<"_x(post_s"<<stage-1<<"_r"<<row<<"_x'length-1 downto post_s"<<stage-1<<"_r"<<row<<"_x'length-"<<wIn<<")),"<<wIn+1<<")+resize("<<(negX[row+N/(pow(2,stage))].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_x(post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_x'length-1 downto post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_x'length-"<<wIn<<")),"<<wIn+1<<"));"<<endl;
//                    vhdl << tab << tmpSignaly.str() << " <= " << "std_logic_vector(resize("<<(negY[row].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row<<"_y(post_s"<<stage-1<<"_r"<<row<<"_y'length-1 downto post_s"<<stage-1<<"_r"<<row<<"_y'length-"<<wIn<<")),"<<wIn+1<<")+resize("<<(negY[row+N/(pow(2,stage))].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_y(post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_y'length-1 downto post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_y'length-"<<wIn<<")),"<<wIn+1<<"));"<<endl;
//                }else{
//                    vhdl << tab << tmpSignalx.str() << " <= " << "std_logic_vector(resize("<<(negX[row-N/(pow(2,stage))].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row-N/(pow(2,stage))<<"_x(post_s"<<stage-1<<"_r"<<row-N/(pow(2,stage))<<"_x'length-1 downto post_s"<<stage-1<<"_r"<<row-N/(pow(2,stage))<<"_x'length-"<<wIn<<")),"<<wIn+1<<")-resize("<<(negX[row].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row<<"_x(post_s"<<stage-1<<"_r"<<row<<"_x'length-1 downto post_s"<<stage-1<<"_r"<<row<<"_x'length-"<<wIn<<")),"<<wIn+1<<"));"<<endl;
//                    vhdl << tab << tmpSignaly.str() << " <= " << "std_logic_vector(resize("<<(negY[row-N/(pow(2,stage))].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row-N/(pow(2,stage))<<"_y(post_s"<<stage-1<<"_r"<<row-N/(pow(2,stage))<<"_y'length-1 downto post_s"<<stage-1<<"_r"<<row-N/(pow(2,stage))<<"_y'length-"<<wIn<<")),"<<wIn+1<<")-resize("<<(negY[row].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row<<"_y(post_s"<<stage-1<<"_r"<<row<<"_y'length-1 downto post_s"<<stage-1<<"_r"<<row<<"_y'length-"<<wIn<<")),"<<wIn+1<<"));"<<endl;
//                }
//            }
//            else{ //last stage feeds outputs and doesn't include rotations
//                if(!((row) & (1<<((int)(log2(N)-stage))))){
//                    vhdl << tab << declare(tmpSignalx.str(),wIn+1) << " <= " << "std_logic_vector(resize("<<(negX[row].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row<<"_x(post_s"<<stage-1<<"_r"<<row<<"_x'length-1 downto post_s"<<stage-1<<"_r"<<row<<"_x'length-"<<wIn<<")),"<<wIn+1<<")+resize("<<(negX[row+N/(pow(2,stage))].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_x(post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_x'length-1 downto post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_x'length-"<<wIn<<")),"<<wIn+1<<"));"<<endl;
//                    vhdl << tab << declare(tmpSignaly.str(),wIn+1) << " <= " << "std_logic_vector(resize("<<(negY[row].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row<<"_y(post_s"<<stage-1<<"_r"<<row<<"_y'length-1 downto post_s"<<stage-1<<"_r"<<row<<"_y'length-"<<wIn<<")),"<<wIn+1<<")+resize("<<(negY[row+N/(pow(2,stage))].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_y(post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_y'length-1 downto post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_y'length-"<<wIn<<")),"<<wIn+1<<"));"<<endl;
//                }else{
//                    vhdl << tab << declare(tmpSignalx.str(),wIn+1) << " <= " << "std_logic_vector(resize("<<(negX[row-N/(pow(2,stage))].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row-N/(pow(2,stage))<<"_x(post_s"<<stage-1<<"_r"<<row-N/(pow(2,stage))<<"_x'length-1 downto post_s"<<stage-1<<"_r"<<row-N/(pow(2,stage))<<"_x'length-"<<wIn<<")),"<<wIn+1<<")-resize("<<(negX[row].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row<<"_x(post_s"<<stage-1<<"_r"<<row<<"_x'length-1 downto post_s"<<stage-1<<"_r"<<row<<"_x'length-"<<wIn<<")),"<<wIn+1<<"));"<<endl;
//                    vhdl << tab << declare(tmpSignaly.str(),wIn+1) << " <= " << "std_logic_vector(resize("<<(negY[row-N/(pow(2,stage))].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row-N/(pow(2,stage))<<"_y(post_s"<<stage-1<<"_r"<<row-N/(pow(2,stage))<<"_y'length-1 downto post_s"<<stage-1<<"_r"<<row-N/(pow(2,stage))<<"_y'length-"<<wIn<<")),"<<wIn+1<<")-resize("<<(negY[row].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row<<"_y(post_s"<<stage-1<<"_r"<<row<<"_y'length-1 downto post_s"<<stage-1<<"_r"<<row<<"_y'length-"<<wIn<<")),"<<wIn+1<<"));"<<endl;
//                }
//                //Assigning the output
//                vhdl << tab  <<  "Out_"<< getBitReverse(row,log2(N)) << "_x <= " << tmpSignalx.str() << "(" << wIn << " downto 1);" << endl;
//                vhdl << tab  <<  "Out_"<< getBitReverse(row,log2(N)) << "_y <= " << tmpSignaly.str() << "(" << wIn << " downto 1);" << endl;
//            }
            //put butterflies here
            if (stage==1){ //first stage requires the input
                if(!((row) & (1<<((int)(log2(N)-stage))))){
                    vhdl << tab << tmpSignalx.str() << " <= " << "std_logic_vector(resize(signed(In_"<<row<<"_x),"<<wIn+1<<")+resize(signed(In_"<<row+N/(pow(2,stage))<<"_x),"<<wIn+1<<"));"<<endl;
                    vhdl << tab << tmpSignaly.str() << " <= " << "std_logic_vector(resize(signed(In_"<<row<<"_y),"<<wIn+1<<")+resize(signed(In_"<<row+N/(pow(2,stage))<<"_y),"<<wIn+1<<"));"<<endl;
                }else{
                    vhdl << tab << tmpSignalx.str() << " <= " << "std_logic_vector(resize(signed(In_"<<row-N/(pow(2,stage))<<"_x),"<<wIn+1<<")-resize(signed(In_"<<row<<"_x),"<<wIn+1<<"));"<<endl;
                    vhdl << tab << tmpSignaly.str() << " <= " << "std_logic_vector(resize(signed(In_"<<row-N/(pow(2,stage))<<"_y),"<<wIn+1<<")-resize(signed(In_"<<row<<"_y),"<<wIn+1<<"));"<<endl;
                }
            }
            else if (stage<n) { // intermediate stages

                if(!((row) & (1<<((int)(log2(N)-stage))))){
                    vhdl << tab << tmpSignalx.str() << " <= " << "std_logic_vector(resize("<<(negX[row].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row<<"_x("<<wIn+12-1<<" downto "<<wIn+12-16<<")),"<<wIn+1<<")+resize("<<(negX[row+N/(pow(2,stage))].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_x("<<wIn+12-1<<" downto "<<wIn+12-16<<")),"<<wIn+1<<"));"<<endl;
                    vhdl << tab << tmpSignaly.str() << " <= " << "std_logic_vector(resize("<<(negY[row].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row<<"_y("<<wIn+12-1<<" downto "<<wIn+12-16<<")),"<<wIn+1<<")+resize("<<(negY[row+N/(pow(2,stage))].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_y("<<wIn+12-1<<" downto "<<wIn+12-16<<")),"<<wIn+1<<"));"<<endl;
                }else{
                    vhdl << tab << tmpSignalx.str() << " <= " << "std_logic_vector(resize("<<(negX[row-N/(pow(2,stage))].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row-N/(pow(2,stage))<<"_x("<<wIn+12-1<<" downto "<<wIn+12-16<<")),"<<wIn+1<<")-resize("<<(negX[row].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row<<"_x("<<wIn+12-1<<" downto "<<wIn+12-16<<")),"<<wIn+1<<"));"<<endl;
                    vhdl << tab << tmpSignaly.str() << " <= " << "std_logic_vector(resize("<<(negY[row-N/(pow(2,stage))].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row-N/(pow(2,stage))<<"_y("<<wIn+12-1<<" downto "<<wIn+12-16<<")),"<<wIn+1<<")-resize("<<(negY[row].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row<<"_y("<<wIn+12-1<<" downto "<<wIn+12-16<<")),"<<wIn+1<<"));"<<endl;
                }
            }
            else{ //last stage feeds outputs and doesn't include rotations
                if(!((row) & (1<<((int)(log2(N)-stage))))){
                    vhdl << tab << declare(tmpSignalx.str(),wIn+1) << " <= " << "std_logic_vector(resize("<<(negX[row].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row<<"_x("<<wIn+12-1<<" downto "<<wIn+12-16<<")),"<<wIn+1<<")+resize("<<(negX[row+N/(pow(2,stage))].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_x("<<wIn+12-1<<" downto "<<wIn+12-16<<")),"<<wIn+1<<"));"<<endl;
                    vhdl << tab << declare(tmpSignaly.str(),wIn+1) << " <= " << "std_logic_vector(resize("<<(negY[row].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row<<"_y("<<wIn+12-1<<" downto "<<wIn+12-16<<")),"<<wIn+1<<")+resize("<<(negY[row+N/(pow(2,stage))].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_y("<<wIn+12-1<<" downto "<<wIn+12-16<<")),"<<wIn+1<<"));"<<endl;
                }else{
                    vhdl << tab << declare(tmpSignalx.str(),wIn+1) << " <= " << "std_logic_vector(resize("<<(negX[row-N/(pow(2,stage))].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row-N/(pow(2,stage))<<"_x("<<wIn+12-1<<" downto "<<wIn+12-16<<")),"<<wIn+1<<")-resize("<<(negX[row].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row<<"_x("<<wIn+12-1<<" downto "<<wIn+12-16<<")),"<<wIn+1<<"));"<<endl;
                    vhdl << tab << declare(tmpSignaly.str(),wIn+1) << " <= " << "std_logic_vector(resize("<<(negY[row-N/(pow(2,stage))].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row-N/(pow(2,stage))<<"_y("<<wIn+12-1<<" downto "<<wIn+12-16<<")),"<<wIn+1<<")-resize("<<(negY[row].at(stage-1)?"-":"")<<"signed(post_s"<<stage-1<<"_r"<<row<<"_y("<<wIn+12-1<<" downto "<<wIn+12-16<<")),"<<wIn+1<<"));"<<endl;
                }
                //Assigning the output
                vhdl << tab  <<  "Out_"<< getBitReverse(row,log2(N)) << "_x <= " << tmpSignalx.str() << "(" << wIn << " downto 1);" << endl;
                vhdl << tab  <<  "Out_"<< getBitReverse(row,log2(N)) << "_y <= " << tmpSignaly.str() << "(" << wIn << " downto 1);" << endl;
            }
        }
    }

}


void FullyParallelFFT::emulate(TestCase * tc) {
    /* This function will be used when the TestBench command is used in the command line
           we have to provide a complete and correct emulation of the operator, in order to compare correct output generated by this function with the test input generated by the vhdl code */
}


int FullyParallelFFT::getBitReverse(int value, int bits) {
    int revValue=0;
    for (int i=0;i<bits;i++){
        if (((value) & (1<<i))) revValue+= pow(2,bits-i-1);
    }
    return revValue;
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
