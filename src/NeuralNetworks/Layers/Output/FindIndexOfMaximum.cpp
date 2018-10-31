//
// Created by nfiege on 12/07/18.
//

#include "FindIndexOfMaximum.hpp"

#include "NeuralNetworks/Layers/Pooling/Max2Input.hpp"

namespace flopoco {
    FindIndexOfMaximum::FindIndexOfMaximum(Target* target, unsigned int wIn, unsigned int numberOfInputs)
            : Operator(target), inputNames(), indexNames()
    {
        useNumericStd();
        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="FindIndexOfMaximum";

        // author
        setCopyrightString("Nicolai Fiege, 2018");

        // definition of the name of the operator
        ostringstream name;
        name << "FindIndexOfMaximum_" << wIn << "_" << numberOfInputs;
        setName(name.str());

        // add inputs and output
        for(unsigned int i=0; i<numberOfInputs; i++)
        {
            inputNames.push_back("X"+to_string(i));
            addInput(inputNames[i], wIn);
        }
        this->indexWordSize=ceil(log2(numberOfInputs));
        addOutput("R",this->indexWordSize);
        for(unsigned int i=0; i<numberOfInputs; i++)
        {
            indexNames.push_back("I"+to_string(i));
            this->vhdl << declare("I"+to_string(i),this->indexWordSize) << " <= std_logic_vector(to_unsigned(" << i << "," << this->indexWordSize << "));" << endl;
        }

        // vector for comparator tree
        vector < vector < Max2Input* > > comparators;
        vector < Max2Input* > firstStage;
        comparators.push_back(firstStage);

        // fill vector and declare an output signal for every comparator
        unsigned int stageCounter=0;
        for(unsigned int i=1; i<=numberOfInputs-1; i++)
        {
            if(i>=((unsigned int)1<<(stageCounter+1)))
            {
                stageCounter++;
                vector < Max2Input* > tmp;
                comparators.push_back(tmp);
            }
            Max2Input* op = new Max2Input(target, wIn, true, true, this->indexWordSize);
            comparators[stageCounter].push_back(op);
            this->addSubComponent(op);
        }

        // connect comparators
        unsigned int inputCounter=0;
        unsigned int comparatorCounter=0;
        map < Max2Input*,string > outputSignalNames;
        map < Max2Input*,string > indexSignalNames;
        for(unsigned int i=0; i<comparators.size(); i++)
        {
            unsigned int i0=comparators.size()-1-i;
            for(unsigned int j=0; j<comparators[i0].size(); j++)
            {
                unsigned int j0 = j;
                //outportmap
                outputSignalNames[comparators[i0][j0]]=declare("s"+to_string(comparatorCounter),wIn);
                this->outPortMap(comparators[i0][j0],"R","s"+to_string(comparatorCounter),false);
                indexSignalNames[comparators[i0][j0]]=declare("IndexSignal"+to_string(comparatorCounter),this->indexWordSize);
                this->outPortMap(comparators[i0][j0],"IMax","IndexSignal"+to_string(comparatorCounter),false);


                //inportmap
                if(comparators.size()>i0+1)
                {
                    //connect 1 or 2 comparators of a higher stage
                    if(comparators[i0+1].size()>=(2*(j0+1)))
                    {
                        //connect 2 comparators of a higher stage
                        this->inPortMap(comparators[i0][j0],"X0",outputSignalNames[comparators[i0+1][2*j0]]);
                        this->inPortMap(comparators[i0][j0],"I0",indexSignalNames[comparators[i0+1][2*j0]]);
                        this->inPortMap(comparators[i0][j0],"X1",outputSignalNames[comparators[i0+1][2*j0+1]]);
                        this->inPortMap(comparators[i0][j0],"I1",indexSignalNames[comparators[i0+1][2*j0+1]]);
                    }
                    else if(comparators[i0+1].size()==(2*(j0+1)-1))
                    {
                        //connect 1 comparator of a higher stage and 1 input
                        this->inPortMap(comparators[i0][j0],"X0",outputSignalNames[comparators[i0+1][2*j0]]);
                        this->inPortMap(comparators[i0][j0],"I0",indexSignalNames[comparators[i0+1][2*j0]]);
                        this->inPortMap(comparators[i0][j0],"X1",inputNames[inputCounter]);
                        this->inPortMap(comparators[i0][j0],"I1",indexNames[inputCounter]);
                        inputCounter++;
                    }
                    else
                    {
                        //connect 2 inputs
                        this->inPortMap(comparators[i0][j0],"X0",inputNames[inputCounter]);
                        this->inPortMap(comparators[i0][j0],"I0",indexNames[inputCounter]);
                        inputCounter++;
                        this->inPortMap(comparators[i0][j0],"X1",inputNames[inputCounter]);
                        this->inPortMap(comparators[i0][j0],"I1",indexNames[inputCounter]);
                        inputCounter++;
                    }
                }
                else
                {
                    //connect 2 inputs
                    this->inPortMap(comparators[i0][j0],"X0",inputNames[inputCounter]);
                    this->inPortMap(comparators[i0][j0],"I0",indexNames[inputCounter]);
                    inputCounter++;
                    this->inPortMap(comparators[i0][j0],"X1",inputNames[inputCounter]);
                    this->inPortMap(comparators[i0][j0],"I1",indexNames[inputCounter]);
                    inputCounter++;
                }

                this->vhdl << this->instance(comparators[i0][j0], "comparator_instance"+to_string(comparatorCounter));
                comparatorCounter++;
            }
        }

        nextCycle();
        vhdl << "R <= " << indexSignalNames[comparators[0][0]] << ";" << endl;

        this->setCycle(comparators.size()+1);
    }





    OperatorPtr FindIndexOfMaximum::parseArguments(Target *target, vector<string> &args) {
        int wIn;
        int numberOfInputs;
        UserInterface::parsePositiveInt(args, "wIn", &wIn, false);
        UserInterface::parsePositiveInt(args, "numberOfInputs", &numberOfInputs, false);
        return new FindIndexOfMaximum(target,wIn,numberOfInputs);
    }

    void FindIndexOfMaximum::registerFactory(){
        UserInterface::add(
                "FindIndexOfMaximum", // name
                "Find the index of the input with the maximum value (useful for classification of CNN-outputs)", // description, string
                "NeuralNetworks", // category, from the list defined in UserInterface.cpp
                "", //seeAlso
                // Now comes the parameter description string.
                // Respect its syntax because it will be used to generate the parser and the docs
                // Syntax is: a semicolon-separated list of parameterDescription;
                // where parameterDescription is parameterName (parameterType)[=defaultValue]: parameterDescriptionString
                "wIn(int): word size of the inputs; numberOfInputs(int): the number of inputs that should be compared against each other",
                "",
                FindIndexOfMaximum::parseArguments
        );
    }

}