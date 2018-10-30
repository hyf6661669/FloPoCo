#ifdef HAVE_PAGLIB
// general c++ library for manipulating streams
#include <iostream>
#include <sstream>
#include <fstream>


#include "gmp.h"
#include "mpfr.h"
#include "FloPoCo.hpp"
#include "PrimitiveComponents/Xilinx/Xilinx_GenericAddSub.hpp"
#include "ConstMultPAG/ConstMultPAG.hpp"
#include "FullyParallelFFT.hpp"


using namespace std;
namespace flopoco {




FullyParallelFFT::FullyParallelFFT(Target* target, int wIn_, int fftSize_, int bC_, int b_, int wLE_, string alg_, string rotatorFileName_, string FFTAlgorithmFileName_, bool intPip_)
    : Operator(target),
      wIn(wIn_),
      fftSize(fftSize_),
      bC(bC_),
      b(b_),
      wLE(wLE_),
      alg(alg_),
      rotatorFileName(rotatorFileName_),
      FFTAlgorithmFileName(FFTAlgorithmFileName_),
      intPip(intPip_)
{
  cout << "rotatorFileName=" << rotatorFileName << endl;
  cout << "FFTAlgorithmFileName=" << FFTAlgorithmFileName << endl;
  cout << "alg=" << alg << endl;

    if(alg == "-" && (FFTAlgorithmFileName == "-" || rotatorFileName == "-"))
    {
      THROWERROR("Either alg or rotatorFileName and FFTAlgorithmFileName have to be specified");
    }

    if(alg != "-")
    {
      if((alg != "SR") && (alg != "BT"))
      {
        THROWERROR("Parameter alg has to be either SR or BT");
      }
    }

    if(rotatorFileName == "-")
    {
      rotatorFileName = "data/FFTRotators/rotators_" + alg + "_bestWLe" + to_string(wLE) + "b" + to_string(b) + ".txt";
    }
    std::ifstream rotFile(rotatorFileName);

		if (!rotFile.is_open())
			THROWERROR("Failed to open rotator file " + rotatorFileName);

    if(FFTAlgorithmFileName == "-")
    {
      FFTAlgorithmFileName = "data/FFTAlgorithms/phi" + to_string(fftSize) + "_" + alg + "_WLe" + to_string(wLE) + "b" + to_string(b) + ".txt";
    }
    std::ifstream FFTFile(FFTAlgorithmFileName);

		if (!FFTFile.is_open())
			THROWERROR("Failed to open FFT algorithm file " + FFTAlgorithmFileName);

    // Parse in rotators
    std::string line;
    vector<string> rotators;
    vector<pair<int64_t,int64_t> > rotatorVal;
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
    vector<vector<int> > cWordsize;
    negX.resize(N);
    negY.resize(N);
    cWordsize.resize(N);
    for (unsigned stage=1; stage<n+1; stage++) {
        for (unsigned row=0; row<N; row++) {
            negX[row].resize(n+1,false);
            negY[row].resize(n+1,false);
            cWordsize[row].resize(n+1,0);
        }
    }
    nextCycle();
    //Building the FFT
    vector<bool> realized;
    realized.resize(300,false);
    for (unsigned stage=1; stage<n+1; stage++) {
        int cC = this->getCurrentCycle();
        int maxDepth = 0;
        for (unsigned row=0; row<N; row++) {
            bool swap_outputs = false;
            bool swap_inputs = false;

            stringstream tmpSignalx,tmpSignaly;

            tmpSignalx << "pre_s"<<stage<<"_r"<<row<<"_x";
            tmpSignaly << "pre_s"<<stage<<"_r"<<row<<"_y";

            declare(tmpSignalx.str(),wIn+1);
            declare(tmpSignaly.str(),wIn+1);

            nextCycle();
            // get right output signs and connect rotator
            if (stage<n) {
                // get rotator for current row
                string rotatorString = rotators.at(FFTRealization[row].at(stage-1));
                //put rotators here
                stringstream tmpOutSignalx; tmpOutSignalx << "post_s"<<stage<<"_r"<<row<<"_x";
                stringstream tmpOutSignaly; tmpOutSignaly << "post_s"<<stage<<"_r"<<row<<"_y";

				        REPORT(DEBUG,"implementing rotator " << FFTRealization[row].at(stage-1) << " as " << rotatorString);

                Operator* cMult = new flopoco::ConstMultPAG(target,wIn+1,rotatorString,intPip,false,1,false);
								cMult->setName("rotator",to_string(FFTRealization[row].at(stage-1)));

                list<flopoco::ConstMultPAG::output_signal_info> tmp_output_list = ((flopoco::ConstMultPAG*)(cMult))->GetOutputList();

                list<flopoco::ConstMultPAG::output_signal_info>::iterator first_output=tmp_output_list.begin();
                list<flopoco::ConstMultPAG::output_signal_info>::iterator second_output= first_output;second_output++;

								if ((*first_output).wordsize !=  (*second_output).wordsize) {THROWERROR("Wordsize missmatch at rotator : " << FFTRealization[row].at(stage-1));}
                else {cWordsize[row].at(stage)=(*first_output).wordsize;}


								if (        (*first_output).output_factors[0][0]==rotatorVal[FFTRealization[row].at(stage-1)].second
                        &&  (*first_output).output_factors[0][1] == rotatorVal[FFTRealization[row].at(stage-1)].first
                        &&  (*second_output).output_factors[0][0] == rotatorVal[FFTRealization[row].at(stage-1)].first
                        &&  (*second_output).output_factors[0][1]== -rotatorVal[FFTRealization[row].at(stage-1)].second)
                {
                    swap_outputs = true;
                }
								else if (        (*first_output).output_factors[0][0]==rotatorVal[FFTRealization[row].at(stage-1)].first
                             &&  (*first_output).output_factors[0][1] == -rotatorVal[FFTRealization[row].at(stage-1)].second
                             &&  (*second_output).output_factors[0][0] == -rotatorVal[FFTRealization[row].at(stage-1)].second
                             &&  (*second_output).output_factors[0][1]== -rotatorVal[FFTRealization[row].at(stage-1)].first)
                {
                    negY[row].at(stage)=true;
                }
                else if (    (*first_output).output_factors[0][0]==-rotatorVal[FFTRealization[row].at(stage-1)].first
                             &&  (*first_output).output_factors[0][1] == rotatorVal[FFTRealization[row].at(stage-1)].second
                             &&  (*second_output).output_factors[0][0] == rotatorVal[FFTRealization[row].at(stage-1)].second
                             &&  (*second_output).output_factors[0][1]== rotatorVal[FFTRealization[row].at(stage-1)].first)
                {
                    negX[row].at(stage)=true;
                }
                else if (    (*first_output).output_factors[0][0]==-rotatorVal[FFTRealization[row].at(stage-1)].first
                             &&  (*first_output).output_factors[0][1] == rotatorVal[FFTRealization[row].at(stage-1)].second
                             &&  (*second_output).output_factors[0][0] == -rotatorVal[FFTRealization[row].at(stage-1)].second
                             &&  (*second_output).output_factors[0][1]== -rotatorVal[FFTRealization[row].at(stage-1)].first)
                {
                    negY[row].at(stage)=true;
                    negX[row].at(stage)=true;
                }
                else if (    (*first_output).output_factors[0][0]==rotatorVal[FFTRealization[row].at(stage-1)].first
                             &&  (*first_output).output_factors[0][1] == rotatorVal[FFTRealization[row].at(stage-1)].second
                             &&  (*second_output).output_factors[0][0] == -rotatorVal[FFTRealization[row].at(stage-1)].second
                             &&  (*second_output).output_factors[0][1]== rotatorVal[FFTRealization[row].at(stage-1)].first)
                {
                    swap_inputs = true;
                    swap_outputs = true;
                }

                else if (    (*first_output).output_factors[0][0]==-rotatorVal[FFTRealization[row].at(stage-1)].first
                             &&  (*first_output).output_factors[0][1] == -rotatorVal[FFTRealization[row].at(stage-1)].second
                             &&  (*second_output).output_factors[0][0] == rotatorVal[FFTRealization[row].at(stage-1)].second
                             &&  (*second_output).output_factors[0][1]== -rotatorVal[FFTRealization[row].at(stage-1)].first)
                {
                    swap_inputs = true;
                    swap_outputs = true;
                    negX[row].at(stage)=true;
                    negY[row].at(stage)=true;
                }
                else if (    (*first_output).output_factors[0][0]== -rotatorVal[FFTRealization[row].at(stage-1)].second
                             &&  (*first_output).output_factors[0][1] == -rotatorVal[FFTRealization[row].at(stage-1)].first
                             &&  (*second_output).output_factors[0][0] == rotatorVal[FFTRealization[row].at(stage-1)].first
                             &&  (*second_output).output_factors[0][1]== -rotatorVal[FFTRealization[row].at(stage-1)].second)
                {
                    swap_outputs = true;
                    negY[row].at(stage)=true;
                }
                else if (    (*first_output).output_factors[0][0]== -rotatorVal[FFTRealization[row].at(stage-1)].second
                             &&  (*first_output).output_factors[0][1] == -rotatorVal[FFTRealization[row].at(stage-1)].first
                             &&  (*second_output).output_factors[0][0] == -rotatorVal[FFTRealization[row].at(stage-1)].first
                             &&  (*second_output).output_factors[0][1]== rotatorVal[FFTRealization[row].at(stage-1)].second)
                {
                    swap_outputs = true;
                    negX[row].at(stage)=true;
                    negY[row].at(stage)=true;

                }
                else if (    (*first_output).output_factors[0][0]==rotatorVal[FFTRealization[row].at(stage-1)].first
                             &&  (*first_output).output_factors[0][1] == -rotatorVal[FFTRealization[row].at(stage-1)].second
                             &&  (*second_output).output_factors[0][0] == rotatorVal[FFTRealization[row].at(stage-1)].second
                             &&  (*second_output).output_factors[0][1]== rotatorVal[FFTRealization[row].at(stage-1)].first)
                {
                    //cerr << "ok for " << FFTRealization[row].at(stage-1) << endl;
                }
                else {
										THROWERROR("!!! no rule to build rotator " << FFTRealization[row].at(stage-1));
								}


                if (!swap_inputs)
                {
                    inPortMap(cMult,"x_in0",tmpSignalx.str());
                    inPortMap(cMult,"x_in1",tmpSignaly.str());
                }
                else
                {
                    inPortMap(cMult,"x_in0",tmpSignaly.str());
                    inPortMap(cMult,"x_in1",tmpSignalx.str());
                }

                if (!swap_outputs)
                {
                    outPortMap(cMult,(*first_output).signal_name,tmpOutSignalx.str(),true);
                    outPortMap(cMult,(*second_output).signal_name,tmpOutSignaly.str(),true);
                }
                else
                {
                    outPortMap(cMult,(*first_output).signal_name,tmpOutSignaly.str(),true);
                    outPortMap(cMult,(*second_output).signal_name,tmpOutSignalx.str(),true);
                }

                addSubComponent(cMult);

                if (cMult->getPipelineDepth()>maxDepth) maxDepth=cMult->getPipelineDepth();

                stringstream rotatorName;
                rotatorName << "rotator_" << FFTRealization[row].at(stage-1) << "at_s" << stage << "_r" << row;

                vhdl << instance(cMult, rotatorName.str());
            }
            this->setCycle(cC);
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
            else { // intermediate stages

								if(!((row) & (1<<((int)(log2(N)-stage)))))
								{

                    if (cWordsize[row].at(stage-1) > bC+wIn)
                    {
												vhdl << tab<< declare(join("post_s",stage-1,"_r",row,"_x_l_v"),wIn) << "<= std_logic_vector(signed(post_s"<<stage-1<<"_r"<<row<<"_x(post_s"<<stage-1<<"_r"<<row<<"_x'length-1 downto post_s"<<stage-1<<"_r"<<row<<"_x'length-" << wIn << ")));" << endl;
                        vhdl << tab<< declare(join("post_s",stage-1,"_r",row,"_x_l"),wIn+1) << "<= std_logic_vector(signed(shift_left(resize(signed(post_s"<<stage-1<<"_r"<<row<<"_x_l_v),"<<wIn+1<<"),1))) ;" << endl;
                    }
                    else
                    {
												vhdl << tab<< declare(join("post_s",stage-1,"_r",row,"_x_l"),wIn+1) << "<= std_logic_vector(resize(signed(post_s"<<stage-1<<"_r"<<row<<"_x(post_s"<<stage-1<<"_r"<<row<<"_x'length-1 downto post_s"<<stage-1<<"_r"<<row<<"_x'length-" << wIn << ")),"<<wIn+1<<"));" << endl;
                    }

                    if (cWordsize[row+N/(pow(2,stage))].at(stage-1) > bC+wIn)
                    {
												vhdl << tab<< declare(join("post_s",stage-1,"_r",row+N/(pow(2,stage)),"_x_l_v"),wIn) <<"<= std_logic_vector(signed(post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_x(post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_x'length-1 downto post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_x'length-" << wIn << ")));"<<endl;
                        vhdl << tab<< declare(join("post_s",stage-1,"_r",row+N/(pow(2,stage)),"_x_l"),wIn+1) << "<= std_logic_vector(signed(shift_left(resize(signed(post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_x_l_v),"<<wIn+1<<"),1))) ;" << endl;
                    }
                    else
                    {
												vhdl << tab<< declare(join("post_s",stage-1,"_r",row+N/(pow(2,stage)),"_x_l"),wIn+1) <<"<= std_logic_vector(resize(signed(post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_x(post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_x'length-1 downto post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_x'length-" << wIn << ")),"<<wIn+1<<"));"<<endl;
                    }

                    if (negX[row].at(stage-1) && negX[row+N/(pow(2,stage))].at(stage-1) && target->getID()=="Virtex6") //sub with 2 neg. inputs
                    {
                        Xilinx_GenericAddSub *addsub = new Xilinx_GenericAddSub( target, wIn+1 , true );
                        addSubComponent( addsub );
                        inPortMap( addsub, "x_i", join( "post_s",stage-1,"_r",row,"_x_l" ) );
                        inPortMap( addsub, "y_i", join( "post_s",stage-1,"_r",row+N/(pow(2,stage)),"_x_l") );
                        inPortMapCst( addsub, "neg_x_i", "'1'" );
                        inPortMapCst( addsub, "neg_y_i", "'1'" );
                        outPortMap( addsub, "sum_o", tmpSignalx.str() , false );
                        vhdl << instance( addsub, join( "add_pre_s" , stage , "_r" , row , "_x") );
                    }
                    else {
                        vhdl << tab << tmpSignalx.str() << " <= " << "std_logic_vector("<<(negX[row].at(stage-1)?"-":"")<<"resize(signed(post_s"<<stage-1<<"_r"<<row<<"_x_l),"<<wIn+1<<")"<<(negX[row+N/(pow(2,stage))].at(stage-1)?"-":"+")<<"resize(signed(post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_x_l),"<<wIn+1<<"));"<<endl;
                    }

                    if (cWordsize[row].at(stage-1) > bC+wIn)
                    {
												vhdl << tab<< declare(join("post_s",stage-1,"_r",row,"_y_l_v"),wIn) << "<= std_logic_vector(signed(post_s"<<stage-1<<"_r"<<row<<"_y(post_s"<<stage-1<<"_r"<<row<<"_y'length-1 downto post_s"<<stage-1<<"_r"<<row<<"_y'length-" << wIn << ")));" << endl;
                        vhdl << tab<< declare(join("post_s",stage-1,"_r",row,"_y_l"),wIn+1) << "<= std_logic_vector(signed(shift_left(resize(signed(post_s"<<stage-1<<"_r"<<row<<"_y_l_v),"<<wIn+1<<"),1))) ;" << endl;
                    }
                    else
                    {
												vhdl << tab<< declare(join("post_s",stage-1,"_r",row,"_y_l"),wIn+1) << "<= std_logic_vector(resize(signed(post_s"<<stage-1<<"_r"<<row<<"_y(post_s"<<stage-1<<"_r"<<row<<"_y'length-1 downto post_s"<<stage-1<<"_r"<<row<<"_y'length-" << wIn << ")),"<<wIn+1<<"));" << endl;
                    }

                    if (cWordsize[row+N/(pow(2,stage))].at(stage-1) > bC+wIn)
                    {
												vhdl << tab<< declare(join("post_s",stage-1,"_r",row+N/(pow(2,stage)),"_y_l_v"),wIn) <<"<= std_logic_vector(signed(post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_y(post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_y'length-1 downto post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_y'length-" << wIn << ")));"<<endl;
                        vhdl << tab<< declare(join("post_s",stage-1,"_r",row+N/(pow(2,stage)),"_y_l"),wIn+1) << "<= std_logic_vector(signed(shift_left(resize(signed(post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_y_l_v),"<<wIn+1<<"),1))) ;" << endl;
                    }
                    else
                    {
												vhdl << tab<< declare(join("post_s",stage-1,"_r",row+N/(pow(2,stage)),"_y_l"),wIn+1) <<"<= std_logic_vector(resize(signed(post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_y(post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_y'length-1 downto post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_y'length-" << wIn << ")),"<<wIn+1<<"));"<<endl;
                    }

                    if (negY[row].at(stage-1) && negY[row+N/(pow(2,stage))].at(stage-1) && target->getID()=="Virtex6" ) //sub with 2 neg. inputs
                    {
                        Xilinx_GenericAddSub *addsub = new Xilinx_GenericAddSub( target, wIn+1, true );
                        addSubComponent( addsub );
                        inPortMap( addsub, "x_i", join( "post_s",stage-1,"_r",row,"_y_l" ) );
                        inPortMap( addsub, "y_i", join( "post_s",stage-1,"_r",row+N/(pow(2,stage)),"_y_l") );
                        inPortMapCst( addsub, "neg_x_i", "'1'" );
                        inPortMapCst( addsub, "neg_y_i", "'1'" );
                        outPortMap( addsub, "sum_o", tmpSignaly.str() , false );
                        vhdl << instance( addsub, join( "add_pre_s" , stage , "_r" , row , "_y") );
                    }
                    else {
                        vhdl << tab << tmpSignaly.str() << " <= " << "std_logic_vector("<<(negY[row].at(stage-1)?"-":"")<<"resize(signed(post_s"<<stage-1<<"_r"<<row<<"_y_l),"<<wIn+1<<")"<<(negY[row+N/(pow(2,stage))].at(stage-1)?"-":"+")<<"resize(signed(post_s"<<stage-1<<"_r"<<row+N/(pow(2,stage))<<"_y_l),"<<wIn+1<<"));"<<endl;
                    }
                }
								else
								{

                    if (negX[row-N/(pow(2,stage))].at(stage-1) && !negX[row].at(stage-1) && target->getID()=="Virtex6") //sub with 2 neg. inputs
                    {
                        Xilinx_GenericAddSub *addsub = new Xilinx_GenericAddSub( target, wIn+1 , true );
                        addSubComponent( addsub );
                        inPortMap( addsub, "x_i", join( "post_s",stage-1,"_r",row-N/(pow(2,stage)),"_x_l" ) );
                        inPortMap( addsub, "y_i", join( "post_s",stage-1,"_r",row,"_x_l") );
                        inPortMapCst( addsub, "neg_x_i", "'1'" );
                        inPortMapCst( addsub, "neg_y_i", "'1'" );
                        outPortMap( addsub, "sum_o", tmpSignalx.str() , false );
                        vhdl << instance( addsub, join( "add_pre_s" , stage , "_r" , row , "_x") );
                    }
                    else
                    {
                        vhdl << tab << tmpSignalx.str() << " <= " << "std_logic_vector("<<(negX[row-N/(pow(2,stage))].at(stage-1)?"-":"")<<"resize(signed(post_s"<<stage-1<<"_r"<<row-N/(pow(2,stage))<<"_x_l),"<<wIn+1<<")"<<(negX[row].at(stage-1)?"+":"-")<<"resize(signed(post_s"<<stage-1<<"_r"<<row<<"_x_l),"<<wIn+1<<"));"<<endl;
                    }
                    if (negY[row-N/(pow(2,stage))].at(stage-1) && !negY[row].at(stage-1) && target->getID()=="Virtex6") //sub with 2 neg. inputs
                    {
                        Xilinx_GenericAddSub *addsub = new Xilinx_GenericAddSub( target, wIn+1 , true );
                        addSubComponent( addsub );
                        inPortMap( addsub, "x_i", join( "post_s",stage-1,"_r",row-N/(pow(2,stage)),"_y_l" ) );
                        inPortMap( addsub, "y_i", join( "post_s",stage-1,"_r",row,"_y_l") );
                        inPortMapCst( addsub, "neg_x_i", "'1'" );
                        inPortMapCst( addsub, "neg_y_i", "'1'" );
                        outPortMap( addsub, "sum_o", tmpSignaly.str() , false );
                        vhdl << instance( addsub, join( "add_pre_s" , stage , "_r" , row , "_y") );
                    }
                    else
                    {
                        vhdl << tab << tmpSignaly.str() << " <= " << "std_logic_vector("<<(negY[row-N/(pow(2,stage))].at(stage-1)?"-":"")<<"resize(signed(post_s"<<stage-1<<"_r"<<row-N/(pow(2,stage))<<"_y_l),"<<wIn+1<<")"<<(negY[row].at(stage-1)?"+":"-")<<"resize(signed(post_s"<<stage-1<<"_r"<<row<<"_y_l),"<<wIn+1<<"));"<<endl;
                    }
                }


                if (stage==n) { // last stage -> assign outputs

                    //Assigning the output
                    vhdl << tab  <<  "Out_"<< getBitReverse(row,log2(N)) << "_x <= " << tmpSignalx.str() << "(" << wIn << " downto 1);" << endl;
                    vhdl << tab  <<  "Out_"<< getBitReverse(row,log2(N)) << "_y <= " << tmpSignaly.str() << "(" << wIn << " downto 1);" << endl;
                }
            }
        }
        if(stage<n)
        {
            this->setCycle(cC+maxDepth+2);
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
    int wIn, bC, wLE, b, fftSize;
    string rotatorFileName, FFTAlgorithmFileName, alg;
    bool intPip;
    UserInterface::parseInt(args, "wIn", &wIn); // param0 has a default value, this method will recover it if it doesnt't find it in args,
    UserInterface::parseInt(args, "fftSize", &fftSize);
    UserInterface::parseInt(args, "bC", &bC);
    UserInterface::parseInt(args, "wLE", &wLE);
    UserInterface::parseInt(args, "b", &b);
    UserInterface::parseString(args, "alg", &alg);
    UserInterface::parseString(args, "rotatorFileName", &rotatorFileName);
    UserInterface::parseString(args, "FFTAlgorithmFileName", &FFTAlgorithmFileName);
    UserInterface::parseBoolean(args, "intPip", &intPip);
    return new FullyParallelFFT(target, wIn, fftSize, bC, b, wLE, alg, rotatorFileName, FFTAlgorithmFileName, intPip);
}

void FullyParallelFFT::registerFactory(){
    UserInterface::add("FullyParallelFFT", // name
                       "Generator for a fully parallel multiplier-less FFT implementation using shift and add", // description, string
                       "FiltersEtc", // category, from the list defined in UserInterface.cpp
                       "", //seeAlso
                       // Now comes the parameter description string.
                       // Respect its syntax because it will be used to generate the parser and the docs
                       // Syntax is: a semicolon-separated list of parameterDescription;
                       // where parameterDescription is parameterName (parameterType)[=defaultValue]: parameterDescriptionString
                       "wIn(int): input word size;  \
                       alg(string): algorithm, may be BT for binary tree or SR for the split-radix algorithm; \
                       fftSize(int)=32: size of FFT;\
                       bC(int)=12: constant word length of rotators;\
                       wLE(int)=16: Effective word length;\
                       b(int)=20: coefficient word length;\
                       rotatorFileName(string)=-: rotator file name (see data/FFTRotators), overwrites the path determined by alg; \
                       FFTAlgorithmFileName(string)=-: FFT algorithm file name (see data/FFTAlgorithms), overwrites the path determined by alg; \
                       intPip(bool)=true: activate interal rotator pipelining",
                       // More documentation for the HTML pages. If you want to link to your blog, it is here.
                       "A fully parallel multiplierless FFT implementation",
                       FullyParallelFFT::parseArguments
            ) ;
}

}//namespace
#endif // HAVE_PAGLIB
