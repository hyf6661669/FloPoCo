// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

// include the header of the Operator
#include "LayerArguments.hpp"


using namespace std;
namespace flopoco
{
    LayerArguments::LayerArguments()
    {
        //set all values on default values!
        this->layerType = "No Type";
        this->number = -1;
        this->coreSize = -1;
        this->inputHeight = -1;
        this->inputWidth = -1;
        this->inputDepth = -1;
        this->wordSize = -1;
        this->fraction = -1;
        this->weightWordSize = -1;
        this->weightFraction = -1;
        this->numberOfOutputFeatures = -1;
        //weights-vector is still empty!
        this->id = "0";
    }
    LayerArguments::LayerArguments(string layerType_,int number_, int coreSize_, int inputHeight_, int inputWidth_, int inputDepth_, int wordSize_, int fraction_, int weightWordSize_, int weightFraction_, int numberOfOutputFeatures_, vector<double> weights_, string id_)
    {
        this->layerType = layerType_;
        this->number = number_;
        this->coreSize = coreSize_;
        this->inputHeight = inputHeight_;
        this->inputWidth = inputWidth_;
        this->inputDepth = inputDepth_;
        this->wordSize = wordSize_;
        this->fraction = fraction_;
        this->weightWordSize = weightWordSize_;
        this->weightFraction = weightFraction_;
        this->numberOfOutputFeatures = numberOfOutputFeatures_;
        this->weights = weights_;
        this->id = id_;
    }
    string LayerArguments::getLayerType()
    {
        if(this->layerType=="No Type")
        {
            cout << "LayerArguments.getLayerType, Warning: Requested layerType isn't set" << endl;
        }
        return this->layerType;
    }
    int LayerArguments::getNumber()
    {
        if(this->number==-1)
        {
            cout << "LayerArguments.getNumber, Warning: Requested number isn't set" << endl;
        }
        return this->number;
    }
    int LayerArguments::getCoreSize()
    {
        if(this->layerType!="Convolutional" && this->layerType!="Pooling")
        {
            cout << "LayerArguments.getCoreSize, Warning: Requesting coreSize on layerType '" << this->layerType << "'" << endl;
        }
        if(this->coreSize==-1)
        {
            cout << "LayerArguments.getCoreSize, Warning: Requested coreSize isn't set" << endl;
        }
        return this->coreSize;
    }
    int LayerArguments::getInputHeight()
    {
        if(this->inputHeight==-1)
        {
            cout << "LayerArguments.getInputHeight, Warning: Requested inputHeight isn't set" << endl;
        }
        return this->inputHeight;
    }
    int LayerArguments::getInputWidth()
    {
        if(this->inputWidth==-1)
        {
            cout << "LayerArguments.getInputWidth, Warning: Requested inputWidth isn't set" << endl;
        }
        return this->inputWidth;
    }
    int LayerArguments::getInputDepth()
    {
        if(this->inputDepth==-1)
        {
            cout << "LayerArguments.getInputDepth, Warning: Requested inputDepth isn't set" << endl;
        }
        return this->inputDepth;
    }
    int LayerArguments::getWordSize()
    {
        if(this->wordSize==-1)
        {
            cout << "LayerArguments.getWordSize, Warning: Requested wordSize isn't set" << endl;
        }
        return this->wordSize;
    }
    int LayerArguments::getFraction()
    {
        if(this->fraction==-1)
        {
            cout << "LayerArguments.getFraction, Warning: Requested fraction isn't set" << endl;
        }
        return this->fraction;
    }
    int LayerArguments::getWeightWordSize()
    {
        if(this->weightWordSize==-1)
        {
            cout << "LayerArguments.getWeightWordSize, Warning: Requested weightWordSize isn't set" << endl;
        }
        return this->weightWordSize;
    }
    int LayerArguments::getWeightFraction()
    {
        if(this->weightFraction==-1)
        {
            cout << "LayerArguments.getWeightFraction, Warning: Requested weightFraction isn't set" << endl;
        }
        return this->weightFraction;
    }
    int LayerArguments::getNumberOfOutputFeatures()
    {
        if(this->numberOfOutputFeatures==-1)
        {
            cout << "LayerArguments.getNumberOfOutputFeatures, Warning: Requested getNumberOfOutputFeatures isn't set" << endl;
        }
        return this->numberOfOutputFeatures;
    }
    vector<double> LayerArguments::LayerArguments::getWeights()
    {
        if(this->weights.size()==0)
        {
            cout << "LayerArguments.getWeights, Warning: Requested weights-vector is empty" << endl;
        }
        return this->weights;
    }
    double LayerArguments::getSpecificWeight(unsigned int index)
    {
        if(this->weights.size()<=index)
        {
            cout << "LayerArguments.getSpecificWeight, Warning: Requested weights-vector is smaller than given index" << endl;
            return 0;
        }
        return this->weights[index];
    }
    double LayerArguments::getConvWeight(unsigned int inputFeature, unsigned int outputFeature, unsigned int index)
    {
        if(this->layerType!="Convolutional")
        {
            cout << "Requested a Convolutional Weight, but requested LayerType is '" << this->getLayerType() << "' instead of 'Convolutional'";
            return 0;
        }
        if(this->weights.size()<=((inputFeature*this->getCoreSize()*this->getCoreSize()*this->getNumberOfOutputFeatures())+(outputFeature*this->getCoreSize()*this->getCoreSize())+(index)))
        {
            cout << "LayerArguments.getConvWeight, Warning: Requested weights-vector is smaller than given variables" << endl;
            return 0;
        }
        return weights[((inputFeature*this->getCoreSize()*this->getCoreSize()*this->getNumberOfOutputFeatures())+(outputFeature*this->getCoreSize()*this->getCoreSize())+(index))];
    }
    string LayerArguments::LayerArguments::getId()
    {
        return this->id;
    }

    void LayerArguments::setLayerType(string lt)
    {
        this->layerType=lt;
    }
    void LayerArguments::setNumber(int n)
    {
        this->number=n;
    }
    void LayerArguments::setCoreSize(int cs)
    {
        this->coreSize=cs;
    }
    void LayerArguments::setInputHeight(int ih)
    {
        this->inputHeight=ih;
    }
    void LayerArguments::setInputWidth(int iw)
    {
        this->inputWidth=iw;
    }
    void LayerArguments::setInputDepth(int id)
    {
        this->inputDepth=id;
    }
    void LayerArguments::setWordSize(int ws)
    {
        this->wordSize=ws;
    }
    void LayerArguments::setFraction(int f)
    {
        this->fraction=f;
    }
    void LayerArguments::setNumberOfOutputFeatures(int o)
    {
        this->numberOfOutputFeatures=o;
    }
    void LayerArguments::setWeightWordSize(int wws)
    {
        this->weightWordSize=wws;
    }
    void LayerArguments::setWeightFraction(int wf)
    {
        this->weightFraction=wf;
    }
    void LayerArguments::setWeights(vector<double> w)
    {
        this->weights.resize(w.size());
        this->weights=w;
    }
    void LayerArguments::setConvWeights(vector<vector<vector<double>>> w)
    {
        vector <double> tmpVec;
        for(auto inputFIt : w)
        {
            for(auto outputFIt : inputFIt)
            {
                for(auto indexIt : outputFIt)
                {
                    tmpVec.push_back(indexIt);
                }
            }
        }
        this->weights=tmpVec;
    }
    void LayerArguments::addWeight(double w)
    {
        this->weights.push_back(w);
    }
    void LayerArguments::setId(string i)
    {
        this->id=i;
    }

}//namespace flopoco
