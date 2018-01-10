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
        this->paddingTop = -1;
        this->paddingBot = -1;
        this->paddingLeft = -1;
        this->paddingRight = -1;
        this->paddingType = "";
        this->inputFeaturesParallel = false;
        this->outputFeaturesParallel = false;
        this->activationFunction = "No Activation Function";
        this->id = "0";
    }
    LayerArguments::LayerArguments(string layerType_, int coreSize_, int inputHeight_, int inputWidth_, int inputDepth_, int wordSize_, int fraction_, int weightWordSize_, int weightFraction_, int numberOfOutputFeatures_, vector<double> weights_, int padding_, string paddingType_, bool inputFeaturesParallel_, bool outputFeaturesParallel_, string activationFunction_, int stride_, string id_)
    {
        this->layerType = layerType_;
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
        this->paddingTop = padding_;
        this->paddingBot = padding_;
        this->paddingLeft = padding_;
        this->paddingRight = padding_;
        this->paddingType = paddingType_;
        this->inputFeaturesParallel = inputFeaturesParallel_;
        this->outputFeaturesParallel = outputFeaturesParallel_;
        this->activationFunction = activationFunction_;
        this->stride=stride_;
        this->id = id_;
    }

    void LayerArguments::printLayerArguments()
    {
        cout << endl << endl;
        cout << "### Layer Arguments" << endl;
        cout << "###    layerType: " << this->layerType << endl;
        cout << "###    coreSize: " << this->coreSize << endl;
        cout << "###    inputHeight: " << this->inputHeight << endl;
        cout << "###    inputWidth: " << this->inputWidth << endl;
        cout << "###    inputDepth: " << this->inputDepth << endl;
        cout << "###    wordSize: " << this->wordSize << endl;
        cout << "###    fraction: " << this->fraction << endl;
        cout << "###    weightWordSize: " << this->weightWordSize << endl;
        cout << "###    weightFraction: " << this->weightFraction << endl;
        cout << "###    numberOfOutputFeatures: " << this->numberOfOutputFeatures << endl;
        cout << "###    Want me to print weights? Nope" << endl;
        cout << "###    paddingTop: " << this->paddingTop << endl;
        cout << "###    paddingBot: " << this->paddingBot<< endl;
        cout << "###    paddingLeft: " << this->paddingLeft << endl;
        cout << "###    paddingRight: " << this->paddingRight << endl;
        cout << "###    paddingType: " << this->paddingType << endl;
        cout << "###    inputFeaturesParallel: " << this->inputFeaturesParallel << endl;
        cout << "###    outputFeaturesParallel: " << this->outputFeaturesParallel << endl;
        cout << "###    activationFunction: " << this->activationFunction << endl;
        cout << "###    stride: " << this->stride << endl;
        cout << "###    id: " << this->id << endl;
        cout << endl << endl;
    }
    string LayerArguments::getLayerType()
    {
        cout << this->layerType << endl;
        if(this->layerType=="No Type")
        {
            cout << "LayerArguments.getLayerType, Warning: Requested layerType isn't set" << endl;
        }
        return this->layerType;
    }
    int LayerArguments::getCoreSize()
    {
        if(this->layerType!="Convolutional" && this->layerType!="Pooling")
        {
            cout << "LayerArguments.getCoreSize, Warning: Requesting coreSize on layerType \"" << this->layerType << "\"" << endl;
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
            cout << "LayerArguments.getNumberOfOutputFeatures, Warning: Requested numberOfOutputFeatures isn't set" << endl;
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
            cout << "Requested a Convolutional Weight, but requested LayerType is \"" << this->getLayerType() << "\" instead of \"Convolutional\"" << endl;
            return 0;
        }
        if(this->weights.size()<=((inputFeature*this->getCoreSize()*this->getCoreSize()*this->getNumberOfOutputFeatures())+(outputFeature*this->getCoreSize()*this->getCoreSize())+(index)))
        {
            cout << "LayerArguments.getConvWeight, Warning: Requested weights-vector is smaller than given variables" << endl;
            return 0;
        }
        return weights[((inputFeature*this->getCoreSize()*this->getCoreSize()*this->getNumberOfOutputFeatures())+(outputFeature*this->getCoreSize()*this->getCoreSize())+(index))];
    }

    vector<vector<vector<double> > > LayerArguments::getConvWeights()
    {
        if(this->layerType!="Convolutional")
        {
            cout << "Requested  Convolutional Weights, but requested LayerType is \"" << this->getLayerType() << "\" instead of \"Convolutional\"" << endl;
        }

        vector<vector<vector<double>>> tmp;
        unsigned int inputFeatures = this->getInputDepth();
        unsigned int outputFeatures = this->getNumberOfOutputFeatures();
        unsigned int windowSize = this->getCoreSize();

        for(unsigned int inputC=0; inputC<inputFeatures; inputC++)
        {
            vector <vector <double>> tmp2;
            tmp.push_back(tmp2);
            for(unsigned int outputC=0; outputC<outputFeatures; outputC++)
            {
                vector <double> tmp3;
                tmp[inputC].push_back(tmp3);
                for(unsigned int indexC=0; indexC<windowSize*windowSize; indexC++)
                {
                    tmp[inputC][outputC].push_back(this->getConvWeight(inputC,outputC,indexC));
                }
            }
        }

        return tmp;
    }

    int LayerArguments::getPaddingTop()
    {
        return this->paddingTop;
    }

    int LayerArguments::getPaddingBot()
    {
        return this->paddingBot;
    }

    int LayerArguments::getPaddingLeft()
    {
        return this->paddingLeft;
    }

    int LayerArguments::getPaddingRight()
    {
        return this->paddingRight;
    }

    int LayerArguments::getPadding()
    {
        if(this->paddingTop==-1)
        {
            if(this->paddingBot==-1)
            {
                if(this->paddingLeft==-1)
                {
                    if(this->paddingRight==-1)
                    {
                        cout << "LayerArguments.getPadding, Warning: Requested padding isn't set" << endl;
                        return -1;
                    }
                    return this->paddingRight;
                }
                return this->paddingLeft;
            }
            return this->paddingBot;
        }
        return this->paddingTop;
    }

    string LayerArguments::getPaddingType()
    {
        if(this->paddingType=="")
        {
            cout << "LayerArguments.getPaddingType, Warning: Requested padding type isn't set" << endl;
        }
        return this->paddingType;
    }

    bool LayerArguments::getInputFeaturesParallel()
    {
        return this->inputFeaturesParallel;
    }

    bool LayerArguments::getOutputFeaturesParallel()
    {
        return this->outputFeaturesParallel;
    }

    string LayerArguments::getActivationFunction()
    {
        if(this->activationFunction=="")
        {
            cout << "LayerArguments.getActivationFunction, Warning: Requested activation function isn't set" << endl;
        }
        return this->activationFunction;
    }

    int LayerArguments::getStride()
    {
        if(this->stride==-1)
        {
            cout << "LayerArguments.getStride, Warning: Requested stride isn't set" << endl;
        }
        return this->stride;
    }

    string LayerArguments::getId()
    {
        return this->id;
    }

    void LayerArguments::setLayerType(string lt)
    {
        cout << "###    setting layerType to " << lt << endl;
        this->layerType=lt;
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

    void LayerArguments::setPaddingType(string p)
    {
        this->paddingType=p;
    }

    void LayerArguments::setInputFeaturesParallel(bool i)
    {
        this->inputFeaturesParallel = i;
    }

    void LayerArguments::setOutputFeaturesParallel(bool o)
    {
        this->outputFeaturesParallel = o;
    }

    void LayerArguments::setActivationFunction(string a)
    {
        this->activationFunction = a;
    }

    void LayerArguments::setStride(int s)
    {
        this->stride=s;
    }

    void LayerArguments::setPadding(int p)
    {
        this->paddingTop=p;
        this->paddingBot=(this->coreSize&1==1?p:p-1);
        this->paddingLeft=p;
        this->paddingRight=(this->coreSize&1==1?p:p-1);
    }

    void LayerArguments::setPaddingTop(int p)
    {
        this->paddingTop=p;
    }

    void LayerArguments::setPaddingBot(int p)
    {
        this->paddingBot=p;
    }

    void LayerArguments::setPaddingLeft(int p)
    {
        this->paddingLeft=p;
    }

    void LayerArguments::setPaddingRight(int p)
    {
        this->paddingRight=p;
    }
    void LayerArguments::setId(string i)
    {
        this->id=i;
    }

}//namespace flopoco
