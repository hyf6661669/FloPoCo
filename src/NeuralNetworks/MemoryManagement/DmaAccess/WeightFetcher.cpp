// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "WeightFetcher.hpp"
#include "NeuralNetworks/Utility/ModuloCounter.hpp"
#include "NeuralNetworks/Utility/DataGuard.hpp"
#include "NeuralNetworks/Utility/Register.hpp"
#include "NeuralNetworks/Layers/Layer.hpp"
#include "NeuralNetworks/MemoryManagement/DmaAccess/AddressCalculation.hpp"

#include "PrimitiveComponents/GenericLut.hpp"

using namespace std;
namespace flopoco {

    WeightFetcher::WeightFetcher(Target* target, unsigned int weightsPerAccess_, unsigned int weightWidth, unsigned int startAddress_, unsigned int innerCounterMax, unsigned int outerCounterMax, unsigned int weightsPerOuterCounterStep, bool lutBasedAddressCalculation) :
            Operator(target), weightsPerAccess(weightsPerAccess_), startAddress(startAddress_) {

        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="WeightFetcher";

        // use numeric_std
        useNumericStd();

        // author
        setCopyrightString("Nicolai Fiege, 2018");

        // definition of the name of the operator
        ostringstream name;
        name << "WeightFetcher_" << outerCounterMax << "_" << innerCounterMax << "_" << weightsPerOuterCounterStep << "_" << weightsPerAccess << "_" << weightWidth << "_" << hex << startAddress << dec << (lutBasedAddressCalculation==true?"_lut":"_arith");
        setName(name.str());

        if(weightWidth>32)
        {
            THROWERROR("Width of weights must be smaller than 32!");
        }

        ////////////////////
        // inputs/outputs //
        ////////////////////////////////////
        // connection with DMA Controller //
        ////////////////////////////////////
        addInput(WeightFetcher::getDataInPortName(),32);
        addInput(WeightFetcher::getValidInPortName(),1);
        addInput(WeightFetcher::getLastInPortName(),1);
        addOutput(WeightFetcher::getNewReadAccessPortName(),1);
        addOutput(WeightFetcher::getNumberOfBytesPortName(),23);
        addOutput(WeightFetcher::getStartAddressPortName(),32);

        ///////////////////////////////////////
        // connection with rest of the layer //
        ///////////////////////////////////////
        addOutput(WeightFetcher::getNextWeightsValidPortName(),1);
        unsigned int innerCounterWidth = ceil(log2(innerCounterMax+1));
        if(innerCounterWidth>0)
        {
            addInput(WeightFetcher::getInnerCounterPortName(),innerCounterWidth);
        }
        unsigned int outerCounterWidth = ceil(log2(outerCounterMax+1));
        if(outerCounterWidth>0)
        {
            addInput(WeightFetcher::getOuterCounterPortName(),outerCounterWidth);
        }
        addInput(WeightFetcher::getStartNewReadPortName(),1);
        addInput(WeightFetcher::getWeightsArrivedAtConvCorePortName(),1);
        for(unsigned int i=0; i<this->weightsPerAccess; i++)
        {
            addOutput(this->getNextWeightsPortName(i), weightWidth);
        }

        ///////////////////////////
        // start new read access //
        ///////////////////////////
        this->vhdl << "process(clk)" << endl;
        this->vhdl << "begin" << endl;
        this->vhdl << tab << "if(rising_edge(clk)) then" << endl;
        this->vhdl << tab << tab << "if(rst = '1') then" << endl;
        this->vhdl << tab << tab << tab << declare(WeightFetcher::getStartNewReadPortName()+"_reg",1) << " <= \"0\";" << endl;
        this->vhdl << tab << tab << "else" << endl;
        this->vhdl << tab << tab << tab << WeightFetcher::getStartNewReadPortName() << "_reg <= " << WeightFetcher::getStartNewReadPortName() << ";" << endl;
        this->vhdl << tab << tab << "end if;" << endl;
        this->vhdl << tab << "end if;" << endl;
        this->vhdl << "end process;" << endl;

        this->vhdl << declare(WeightFetcher::getStartNewReadPortName()+"_rise",1) << " <= " << WeightFetcher::getStartNewReadPortName() << " and (not " << WeightFetcher::getStartNewReadPortName() << "_reg);" << endl;

        this->vhdl << "process(clk)" << endl;
        this->vhdl << "begin" << endl;
        this->vhdl << tab << "if(rising_edge(clk)) then" << endl;
        this->vhdl << tab << tab << "if(rst='1') then" << endl;
        this->vhdl << tab << tab << tab << WeightFetcher::getNewReadAccessPortName() << "_temp <= \"0\";" << endl;
        this->vhdl << tab << tab << "else" << endl;
        this->vhdl << tab << tab << tab << "if(" << WeightFetcher::getStartNewReadPortName() << "_rise = \"1\") then" << endl;
        this->vhdl << tab << tab << tab << tab << WeightFetcher::getNewReadAccessPortName() << "_temp <= \"1\";" << endl;
        this->vhdl << tab << tab << tab << "elsif(" << WeightFetcher::getLastInPortName() << " = \"1\") then" << endl;
        this->vhdl << tab << tab << tab << tab << WeightFetcher::getNewReadAccessPortName() << "_temp <= \"0\";" << endl;
        this->vhdl << tab << tab << tab << "end if;" << endl;
        this->vhdl << tab << tab << "end if;" << endl;
        this->vhdl << tab << "end if;" << endl;
        this->vhdl << "end process;" << endl;

        unsigned int weightsPerDataBeat = (unsigned int)floor(32.0/((double)weightWidth)); // number of weights, that fit into a 32-bit signal
        unsigned int dataBeatsPerAccess = (unsigned int)(ceil(((double)this->weightsPerAccess)/((double)weightsPerDataBeat))); // must be 32-bit-aligned (= 4-Byte-Aligned)
        unsigned int dataBeatsPerLastAccess;
        if(weightsPerOuterCounterStep==0)
        {
            // the number of weights on the last access is the same as on every other access
            dataBeatsPerLastAccess = dataBeatsPerAccess;
        }
        else
        {
            unsigned int modValue = weightsPerOuterCounterStep % this->weightsPerAccess;
            if(modValue == 0)
            {
                dataBeatsPerLastAccess = dataBeatsPerAccess;
            }
            else
            {
                dataBeatsPerLastAccess = (unsigned int)(ceil(((double)modValue)/((double)weightsPerDataBeat))); // must be 32-bit-aligned (= 4-Byte-Aligned)
            }
        }

        if(dataBeatsPerLastAccess != dataBeatsPerAccess && innerCounterWidth>0)
            this->vhdl << WeightFetcher::getNumberOfBytesPortName() << " <= std_logic_vector(to_unsigned(" << dataBeatsPerLastAccess*4 << ",23)) when unsigned(" << WeightFetcher::getInnerCounterPortName() << ") = " << innerCounterMax << " else std_logic_vector(to_unsigned(" << dataBeatsPerAccess*4 << ",23));" << endl;
        else if(dataBeatsPerLastAccess != dataBeatsPerAccess && innerCounterWidth==0)
            this->vhdl << WeightFetcher::getNumberOfBytesPortName() << " <= std_logic_vector(to_unsigned(" << dataBeatsPerLastAccess*4 << ",23));" << endl;
        else
            this->vhdl << WeightFetcher::getNumberOfBytesPortName() << " <= std_logic_vector(to_unsigned(" << dataBeatsPerAccess*4 << ",23));" << endl;

        ///////////////////////////////////////////////////////
        // Determine next memory address (LUT or arithmetic) //
        ///////////////////////////////////////////////////////
        int numberOfDelays = 0;
		if(outerCounterMax>0 || innerCounterMax>0)
		{
			if(lutBasedAddressCalculation)
			{

				map<unsigned int, unsigned int> LUTMap = this->getLUTData(outerCounterMax, innerCounterMax, dataBeatsPerAccess, dataBeatsPerLastAccess);
				stringstream lutName;
				lutName << "Start_Address_Lut_" << this->startAddress;
				GenericLut* startAddressLut = new GenericLut(target,lutName.str(),LUTMap,innerCounterWidth+outerCounterWidth,32);
				addSubComponent(startAddressLut);
				for(unsigned int i=0; i<innerCounterWidth; i++)
				{
					string signalForPortName = WeightFetcher::getInnerCounterPortName()+"_"+to_string(i);
					this->vhdl << declare(signalForPortName,1,false) << " <= " << WeightFetcher::getInnerCounterPortName() << "(" << i << ");" << endl;
					inPortMap(startAddressLut,join("i",i),signalForPortName);
				}
				for(unsigned int i=0; i<outerCounterWidth; i++)
				{
					string signalForPortName = WeightFetcher::getOuterCounterPortName()+"_"+to_string(i);
					this->vhdl << declare(signalForPortName,1,false) << " <= " << WeightFetcher::getOuterCounterPortName() << "(" << i << ");" << endl;
					inPortMap(startAddressLut,join("i",innerCounterWidth+i),signalForPortName);
				}
				for(unsigned int i=0; i<32; i++)
				{
					string signalForPortName = "Start_Address_Lut_out_"+to_string(i);
					outPortMap(startAddressLut,join("o",i),signalForPortName,true);
					this->vhdl << WeightFetcher::getStartAddressPortName() << "(" << i << ") <= " << signalForPortName << ";" << endl;
				}
				this->vhdl << instance(startAddressLut, "Start_Address_Lut_instance");
			}
        	else
			{
				int bytesPerAccess = dataBeatsPerAccess*4;
				int bytesPerFeature = dataBeatsPerLastAccess*4;
				if(innerCounterMax > 0)
				{
					bytesPerFeature += innerCounterMax*bytesPerAccess;
				}
				auto addrCalc = new AddressCalculation(target,this->startAddress,innerCounterWidth,outerCounterWidth,bytesPerAccess,bytesPerFeature);
				addSubComponent(addrCalc);
				if(innerCounterWidth>0)
					inPortMap(addrCalc,"InnerCounter",WeightFetcher::getInnerCounterPortName());
				if(outerCounterWidth>0)
					inPortMap(addrCalc,"OuterCounter",WeightFetcher::getOuterCounterPortName());
				outPortMap(addrCalc,"NextAddress",WeightFetcher::getStartAddressPortName(),false);
				this->vhdl << instance(addrCalc,"Address_Calc");
				numberOfDelays = addrCalc->getPipelineDepth();
			}
		}
		else
		{
			this->vhdl << WeightFetcher::getStartAddressPortName() << " <= std_logic_vector(to_unsigned(" << this->startAddress << ",32));" << endl;
		}

        //////////////////////////////////////////////////
        // Read the weights and store them in Registers //
        ///////////////////////////////////////////////////////////////////////////////////////////
        // Create Modulo counter to keep track, which set of weights is being read at the moment //
        ///////////////////////////////////////////////////////////////////////////////////////////
        // Reset Modulo Counter at every new read access
        // Enable Modulo Counter when a valid-signal comes from the DMA AXI Stream
        if(dataBeatsPerAccess>1)
        {
            ModuloCounter* weightSetCounter = new ModuloCounter(target,dataBeatsPerAccess,true);
            inPortMap(weightSetCounter,"enable",WeightFetcher::getValidInPortName());
            inPortMap(weightSetCounter,"manualReset",WeightFetcher::getLastInPortName());
            outPortMap(weightSetCounter,"counter","Modulo_Counter_out",true);
            addSubComponent(weightSetCounter);
            this->vhdl << instance(weightSetCounter, "Weight_Set_Counter_instance");
        }
        else
        {
            this->vhdl << declare("Modulo_Counter_out",1) << " <= \"0\";" << endl;
        }


        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Create Data Guards as enable for the Weight-Registers and set the Guard number on the specific Modulo Counter value //
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        unsigned int positionInDataStream = 0;
        unsigned int moduloCounterValueCounter = 0;
        for(unsigned int i=0; i<weightsPerAccess; i++)
        {
            moduloCounterValueCounter = (unsigned int)floor(((double)i)/((double)weightsPerDataBeat));
            positionInDataStream = i % weightsPerDataBeat;
            unsigned int downto = positionInDataStream*weightWidth;
            unsigned int from = downto+weightWidth-1;
            unsigned int guardWidth=1;
            if(dataBeatsPerAccess>1)
            {
                guardWidth=ceil(log2(dataBeatsPerAccess));
            }
            DataGuard* dat = new DataGuard(target,1,guardWidth,moduloCounterValueCounter);
            addSubComponent(dat);
            this->vhdl << declare("Data_Guard_"+to_string(i)+"_in",1) << " <= " << WeightFetcher::getValidInPortName() << ";" << endl;
            inPortMap(dat,"Data_in","Data_Guard_"+to_string(i)+"_in");
            inPortMap(dat,"Guard_in","Modulo_Counter_out");
            outPortMap(dat,"Data_out",this->getNextWeightsPortName(i)+"_enable",true);
            this->vhdl << instance(dat,"Data_Guard_"+to_string(i)+"_instance");

            this->vhdl << declare("Data_Register_"+to_string(i)+"_in",weightWidth) << " <= " << WeightFetcher::getDataInPortName() << "(" << from << " downto " << downto << ");" << endl;
            Register* reg = new Register(target,weightWidth,1,true);
            addSubComponent(reg);
            inPortMap(reg,"X","Data_Register_"+to_string(i)+"_in");
            inPortMap(reg,"E",this->getNextWeightsPortName(i)+"_enable");
            outPortMap(reg,"R",this->getNextWeightsPortName(i),false);
            this->vhdl << instance(reg,"Weight_Register_"+to_string(i));
        }

        //////////////////////////////////
        // Set weights-are-valid-signal //
        //////////////////////////////////
        this->vhdl << "process(clk)" << endl;
        this->vhdl << "begin" << endl;
        this->vhdl << tab << "if(rising_edge(clk)) then" << endl;
        this->vhdl << tab << tab << "if(rst='1' or " << WeightFetcher::getStartNewReadPortName() << "_rise" << "=\"1\" or " << WeightFetcher::getWeightsArrivedAtConvCorePortName() << "=\"1\") then" << endl;
        this->vhdl << tab << tab << tab << WeightFetcher::getNextWeightsValidPortName() << " <= \"0\";" << endl;
        this->vhdl << tab << tab << "elsif(" << WeightFetcher::getLastInPortName() << "=\"1\") then" << endl;
        this->vhdl << tab << tab << tab << WeightFetcher::getNextWeightsValidPortName() << " <= \"1\";" << endl;
        this->vhdl << tab << tab << "end if;" << endl;
        this->vhdl << tab << "end if;" << endl;
        this->vhdl << "end process;" << endl;

        ///////////////////////////////////////
		// Set output to request new weights //
		///////////////////////////////////////
		declare(WeightFetcher::getNewReadAccessPortName()+"_temp",1);
		if(lutBasedAddressCalculation)
		{
			// lut is not pipelined (obviously...)
			this->vhdl << WeightFetcher::getNewReadAccessPortName() << " <= " << WeightFetcher::getNewReadAccessPortName() << "_temp;" << endl;
		}
		else if(numberOfDelays > 0)
		{
			// set output when address calculation finishes
			auto reg = new Register(target,1,numberOfDelays);
			addSubComponent(reg);
			inPortMap(reg,"X",WeightFetcher::getNewReadAccessPortName()+"_temp");
			outPortMap(reg,"R",WeightFetcher::getNewReadAccessPortName(),false);
			this->vhdl << instance(reg,"NewWeightRequestRegister");
		}
    }

    string WeightFetcher::getNextWeightsPortName(unsigned int weightNumber) const
    {
        if(weightNumber>=this->weightsPerAccess)
        {
            THROWERROR("WeightFetcher.nextWeightsPortName: requested port doesn't exist");
        }
        return "Next_Weights_"+to_string(weightNumber);
    }

    map<unsigned int, unsigned  int> WeightFetcher::getLUTData(unsigned int outerCounterMax, unsigned int innerCounterMax, unsigned int dataBeatsPerAccess, unsigned int dataBeatsPerLastAccess) const
    {
        unsigned int wordSizeInner = ceil(log2(innerCounterMax+1));
        unsigned int wordSizeOuter = ceil(log2(outerCounterMax+1));

        if(wordSizeInner+wordSizeOuter > 32)
        {
            THROWERROR("The weight fetcher can't be implemented, GenericLut-Constructor only supports 32 bits");
        }
        map<unsigned int, unsigned int> returnMe;

        double innerMax = pow(2,wordSizeInner)-1;
        double outerMax = pow(2,wordSizeOuter)-1;
        double LUTMax_d = (outerMax*pow(2,wordSizeInner))+innerMax;
        unsigned int LUTMax = (unsigned int)LUTMax_d;
        unsigned int addressCounter = this->startAddress;
        for(unsigned int i=0; i<=LUTMax; i++)
        {
            unsigned int outerC = i >> (wordSizeInner);
            unsigned int innerC = i-(outerC << wordSizeInner);
            if(innerC <= innerCounterMax && outerC <= outerCounterMax)
            {
                    returnMe[i] = addressCounter;

                    if(innerC == innerCounterMax) addressCounter += (dataBeatsPerLastAccess*4);
                    else addressCounter += (dataBeatsPerAccess*4);
            }
            else
            {
                returnMe[i] = 0;
            }
        }
        return returnMe;
    }

}//namespace flopoco