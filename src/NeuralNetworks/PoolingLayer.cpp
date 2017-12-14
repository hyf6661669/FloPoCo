// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "PoolingLayer.hpp"

//include other Operators
#include "WindowShiftRegister.hpp"
#include "PaddingGenerator.hpp"
#include "PoolingCore.hpp"


using namespace std;
namespace flopoco {




    PoolingLayer::PoolingLayer(Target* target, unsigned int wordSize_, unsigned int horizontalSize_, unsigned int verticalSize_, unsigned int numberOfFeatures_, bool parallel_, int paddingTop_, unsigned int windowSize_, unsigned int stride_, int paddingBot_, int paddingLeft_, int paddingRight_, string paddingType_) :
        Operator(target), wordSize(wordSize_), horizontalSize(horizontalSize_), verticalSize(verticalSize_), numberOfFeatures(numberOfFeatures_), parallel(parallel_), paddingTop(paddingTop_), windowSize(windowSize_), stride(stride_), paddingBot(paddingBot_), paddingLeft(paddingLeft_), paddingRight(paddingRight_), paddingType(paddingType_) {
cout << "### IN POOLINGLAYER CONSTRUCTOR" << endl;
        this->useNumericStd();

        // throw error
        if(stride<1)
        {
            stringstream e;
            e << "The stride must be > 0";
            THROWERROR(e);
        }
		
        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="PoolingLayer";

        // author
        setCopyrightString("Nicolai Fiege, 2017");

        // name of the operator
        ostringstream name;
        name << "PoolingLayer_wordSize_" << wordSize << "_horizontalSize_" << horizontalSize << "_verticalSize_" << verticalSize << "_numberOfFeatures_" << numberOfFeatures << "_parallel_" << (parallel==true?"1":"0") << "_padTop_" << paddingTop << "_windowSize_" << windowSize << "_stride_" << stride << "_padBot_" << paddingBot << "_paddingLeft_" << paddingLeft << "_paddingRight_" << paddingRight << "_paddingType_" << paddingType;
        setName(name.str());
cout << "###    1" << endl;
        // add inputs and outputs
        addInput("newStep",1);
        addOutput("finished",1);
cout << "###    2" << endl;
        unsigned int forLimit=(this->parallel==true?this->numberOfFeatures:1);
        unsigned int maxLatency=0;
        unsigned int poolingCoreMaxLatencyNumber=0;
        for(unsigned int outFC=0; outFC<forLimit; outFC++)
        {
cout << "###   3.1:" << outFC << endl;
            addInput("X_"+to_string(outFC),wordSize);
            addOutput("R_"+to_string(outFC),wordSize);
            addOutput("getNewData"+to_string(outFC),1);
            addOutput("validData_o"+to_string(outFC),1);
            addInput("validData_i"+to_string(outFC),1);
cout << "###    3.2:" << outFC << endl;
            //shift reg
            WindowShiftRegister* winOp = new WindowShiftRegister(target,wordSize,windowSize,horizontalSize);
            addSubComponent(winOp);
            inPortMap(winOp,"X","X_"+to_string(outFC));
            inPortMap(winOp,"enable","validData_i"+to_string(outFC));
            inPortMap(winOp,"newStep","newStep");
            for(unsigned int i=0; i<windowSize*windowSize; i++)
            {
                outPortMap(winOp,winOp->outputNames[i],"winShiftOut_"+to_string(i)+"_"+to_string(outFC),true);
            }
            vhdl << instance(winOp,"winShiftInstance_"+to_string(outFC));
cout << "###    3.3:" << outFC  << endl;
            //padding
            PaddingGenerator* padOp = new PaddingGenerator(target,wordSize,windowSize,horizontalSize,verticalSize,paddingTop,stride,paddingType,paddingBot,paddingLeft,paddingRight,(outFC==0?true:false));
            addSubComponent(padOp);
            for(unsigned int i=0; i<windowSize*windowSize; i++)
            {
                inPortMap(padOp,"X"+to_string(i),"winShiftOut_"+to_string(i)+"_"+to_string(outFC));
                outPortMap(padOp,"R"+to_string(i),"paddingOut_"+to_string(i)+"_"+to_string(outFC),true);
            }
            inPortMap(padOp,"validData_i","validData_i"+to_string(outFC));
            inPortMap(padOp,"newStep","newStep");
            outPortMap(padOp,"getNewData","getNewData"+to_string(outFC),false);
cout << "###    3.4:" << outFC  << endl;
            if(outFC==0)
            {
                outPortMap(padOp,"finished","finished",false);
                outPortMap(padOp,"validData_o","validData_temp",true);
            }
            vhdl << instance(padOp,"paddingInstance_"+to_string(outFC));
cout << "###    3.5:" << outFC  << endl;
            //pooling core
            PoolingCore* pooOp = new PoolingCore(target,wordSize,windowSize);
cout << "###    3.6:" << outFC  << endl;
            addSubComponent(pooOp);
            for(unsigned int i=0; i<windowSize*windowSize; i++)
            {
                inPortMap(pooOp,"X"+to_string(i),"paddingOut_"+to_string(i)+"_"+to_string(outFC));
            }
            outPortMap(pooOp,"R","poolOut_"+to_string(outFC));
//            if(pooOp->getPipelineDepth()>maxLatency)
//            {
//                poolingCoreMaxLatencyNumber=outFC;
//                maxLatency=pooOp->getPipelineDepth();
//            }
cout << "###    3.7:" << outFC  << endl;
            vhdl << instance(pooOp,"poolingInstance_"+to_string(outFC));
        }
cout << "###    4" << endl;
//        if(maxLatency>0)
//        {
//            syncCycleFromSignal("poolOut_"+to_string(poolingCoreMaxLatencyNumber));
//        }
cout << "###    5" << endl;
        for(unsigned int outFC=0; outFC<forLimit; outFC++)
        {
            vhdl << "R_" << outFC << " <= poolOut_" << outFC << ";" << endl;
            vhdl << "validData_o" << outFC << " <= validData_temp;" << endl;
        }
cout << "###    6" << endl;
    }

    unsigned int PoolingLayer::getNumberOfInstances()
    {
        if(this->parallel==false)
        {
            return 1;
        }
        return this->numberOfFeatures;
    }

}//namespace flopoco
