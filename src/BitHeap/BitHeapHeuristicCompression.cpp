

#include "BitHeapHeuristicCompression.hpp"
#include <algorithm>

namespace flopoco
{


    BitHeapHeuristicCompression::BitHeapHeuristicCompression(BitHeap* bh, std::string mode, bool useVariableCompressors)
       : bh_(bh), mode(mode)
#ifdef HAVE_SCIP
       , bitHeapILPCompression(bh)
#endif //HAVE_SCIP
    {
        //standard values
        lowerBound = 1.5;
        useHoles = true;
        generateSolution = false;
        passHeuristicSolution = false;
        passMultipleSolutions = false;
        usePreReduction = false;
        reduceILPStageCount = false;
        useMaxHeuristicStageCount = false;
        optimalAndPreReduction = false;
        differentStages = false;
        minSliceCount = 10000000;
        minFixedStage = 3;
        maxFixedStage = 6;
        paSorting = false;
        buildSingleStages = false;
        numberOfBuildStages = 0;
        useSmootherCompressorPlacing = false;
        dontUsePA = false;
//        useVariableCompressors = false;
        this->useVariableCompressors = useVariableCompressors;
        useCompleteHeuristic = false;
        getLowerBoundsFromBitHeap = false;

        //modified values

        //lowerBound = 0.0;
        //useHoles = false;
        generateSolution = true;
        //passMultipleSolutions = true;

		passHeuristicSolution = true;
        reduceILPStageCount = true;
        //useMaxHeuristicStageCount = true;
        //optimalAndPreReduction = true;
        differentStages = true;
        dontUsePA = false;
        minFixedStage = 1;
        maxFixedStage = 8;
        //paSorting = true;
        //buildSingleStages = true;
        //numberOfBuildStages = 8;
        //useSmootherCompressorPlacing = true;
//        useVariableCompressors = bh->useVariableColumnCompressors; //!!
        //useCompleteHeuristic = true;
        //getLowerBoundsFromBitHeap = true;
        usePreReduction = true;


        //now some dependencies

        if(reduceILPStageCount && !passHeuristicSolution){
            reduceILPStageCount = false;
        }
        cout << "mode is " << mode << endl;
        //get the compressiontype from mode

		buildVariableCompressors();
#ifdef HAVE_SCIP

		bitHeapILPCompression.variableBCompressors = variableBCompressors;

        if(mode.compare("heuristic_parandeh-afshar_modified") == 0){
            bitHeapILPCompression.compressionType = 6;
        }
        else if(mode.compare("heuristic_pa") == 0){
            bitHeapILPCompression.compressionType = 7;
        }
        else if(mode.compare("optimalMinStages") == 0){
            bitHeapILPCompression.compressionType = 5;
        }

#endif //HAVE_SCIP

		vector<BasicCompressor *>* possibleCompressors_ = bh->getPossibleCompressors();
		compOutputWordSizeMax = 0;
		for(unsigned e=0; e < possibleCompressors_->size(); e++){
			if(compOutputWordSizeMax < (*possibleCompressors_)[e]->getOutputSize()){
				compOutputWordSizeMax = (*possibleCompressors_)[e]->getOutputSize();
			}
		}
		compOutputWordSizeMax -= 1;//every compressor has at least one column of inputs
		//REPORT(DEBUG,"compOutputWordSizeMax=" << compOutputWordSizeMax);
    }



    BitHeapHeuristicCompression::~BitHeapHeuristicCompression(){
        //cleanUp();
    }

    int BitHeapHeuristicCompression::getMaxStageCount(int maxHeight){

        /*      old, only for one starting cycle
        int height=2;
        int stageCount=0;

        //the following loop computes the Dadda sequence [Dadda 1965] which is [used] as max. no of stages
        //(typically less are selected by the optimization using higher order compressors)
        while(height < maxHeight)
        {
            //REPORT(DEBUG, "stageCount=" << stageCount << ", height=" << height);
            height = floor(height*3.0/2.0);
            stageCount++;
        }
        return stageCount;

        */

        //new one: "reverse" the dadda algorithm: compute for all bits of the current stage: ceil(bitamount * 2.0/3.0) for the next cycle.
        //in the next cycle add the bits of the current cycle to the remaining


        unsigned int stages = 0;
        vector<int> tempVector(bh_->bits.size());
        for(unsigned w = 0; w < newBits[0].size(); w++){
            tempVector[w] = ceil((((float)newBits[0][w]) *  (2.0/3.0)) - 0.00001);
        }

        //cout << "stage " << stages << endl;
        for(unsigned j = 0; j < tempVector.size(); j++){
            //cout << tempVector[j] << " ";
        }
        //cout << endl << endl;


        while(!maximalTwoBitsInCurrentStage(tempVector) || (stages + 1) < newBits.size()){
            stages++;

            if(stages < newBits.size()){
                for(unsigned w = 0; w < newBits[0].size(); w++){
                    tempVector[w] += newBits[stages][w];
                }
            }
            for(unsigned w = 0; w < newBits[0].size(); w++){
                tempVector[w] = ceil((((float)tempVector[w]) *  (2.0/3.0)) - 0.00001);
            }
            /*
            cout << "stage " << stages << endl;
            for(unsigned j = 0; j < tempVector.size(); j++){
                cout << tempVector[j] << " ";
            }
            cout << endl << endl;
            */
        }


        return stages;

    }

    bool BitHeapHeuristicCompression::maximalTwoBitsInCurrentStage(vector<int> currentBits){

        for(unsigned w = 0; w < currentBits.size(); w++){
            if(currentBits[w] > 2){
                return false;
            }
        }
        return true;
    }

    int BitHeapHeuristicCompression::generateProblem(){

        if(passHeuristicSolution){
            lowerBound = 0.0;
        }
#ifdef HAVE_SCIP
        bitHeapILPCompression.useVariableCompressors = useVariableCompressors;
#endif //HAVE_SCIP



        if(!getLowerBoundsFromBitHeap){
            for(unsigned i = 0; i < (sizeof(lowerBounds) / sizeof(lowerBounds[0])); i++){
                lowerBounds[i] = 10.0;
                if(useCompleteHeuristic || passMultipleSolutions || (mode.compare("heuristic_parandeh-afshar_modified") == 0)){
                    lowerBounds[i] = 0;     //use PA in every stage
                }
            }
        }





        //now fill lowerBounds with real values
        //the preset value is infinity (a.k.a. 10.0)
		if(mode.compare("heuristic_parandeh-afshar_modified") != 0){
            //lowerBounds[0] = 0;
            //lowerBounds[1] = 0;
			//lowerBounds[2] = 1.75;
			//lowerBounds[3] = 1.75;
			//lowerBounds[4] = 0;
			//lowerBounds[5] = 0;
			//lowerBounds[6] = 0;
			//lowerBounds[7] = 0;
		}
        printLowerBounds();

        //set the minFixedStage to at least the number of stages which are solved by the heuristic
        unsigned lowerBoundZeros = 0;
        for(unsigned i = 0; i < (sizeof(lowerBounds) / sizeof(lowerBounds[0])); i++){
            if(fabs(lowerBounds[i]) < 0.0001){
                lowerBoundZeros++;
            }
            else{
                break;
            }
        }
        if(lowerBoundZeros > minFixedStage){
            minFixedStage = lowerBoundZeros;
            if(maxFixedStage < minFixedStage){
                maxFixedStage = minFixedStage;
                cout << "maxFixedStage set to minFixedStage" << endl;
            }
        }
        cout << "minFixedStage = " << minFixedStage << " maxFixedStage = " << maxFixedStage << endl;
		cout << "compOutputWordSizeMax = " << compOutputWordSizeMax << endl;



        //generate newBits

        //first get maximal cycle
        unsigned int maxCycle = 0;
        for(unsigned w = 0; w < bh_->bits.size(); w++){
            unsigned i = 0;
            for(list<WeightedBit*>::iterator it = bh_->bits[w].begin(); it != bh_->bits[w].end(); it++){
                if((*it)->getCycle() > maxCycle){
                    maxCycle = (*it)->getCycle();
                }
            }
        }
        cout << "maxCycle is " << maxCycle << endl;


        for(unsigned int cy = 0; cy <= maxCycle; cy++){
            vector<int> cycleVector(bh_->bits.size() + compOutputWordSizeMax * cy);
            for(unsigned w = 0; w < bh_->bits.size() + compOutputWordSizeMax * cy; w++){
                int amount = 0;
                if(w < bh_->bits.size()){   //every bitposition > bh_->bits.size() == 0
                    for(list<WeightedBit*>::iterator it = bh_->bits[w].begin(); it != bh_->bits[w].end(); it++){
                        if((*it)->getCycle() == cy){
                            amount++;
                        }
                    }
                }
                cycleVector[w] = amount;
            }
            newBits.push_back(cycleVector);
            originalNewBits.push_back(cycleVector);
        }

        //delete stages, where in that stage and all stages before no bits are being added

        unsigned int delStages = 0;
        bool foundInput = false;
        for(unsigned int s = 0; s <= newBits.size(); s++){
            for(unsigned c = 0; c < newBits[s].size(); c++){
                if(newBits[s][c] > 0){
                    foundInput = true;
                    break;
                }
            }
            if(foundInput){
                break;
            }
            delStages++;    //else: new zerostage found
        }

        //if(delStages > 0){  //deleting the zervectors
        if(0){
            newBits.erase(newBits.begin(), newBits.begin() + delStages);
            originalNewBits.erase(originalNewBits.begin(), originalNewBits.begin() + delStages);
            //resize the other vectors
            for(unsigned int s = 0; s < newBits.size(); s++){
                newBits[s].resize(bh_->bits.size() + compOutputWordSizeMax * s);
                originalNewBits[s].resize(bh_->bits.size() + compOutputWordSizeMax * s);
            }
        }



/*
        //debug output:: newBits
        unsigned int sumOfBits = 0;
        for(unsigned i = 0; i < newBits.size(); i++){
            for(unsigned j = 0; j < newBits[i].size(); j++){
                cout << newBits[i][j] << " ";
                sumOfBits += newBits[i][j];
            }
            cout << endl;
        }
        cout << "The total amount of initial Bits in the Bitheap is " << sumOfBits << endl;
*/
        noOfStages_ = getMaxStageCount(bh_->getMaxHeight());

//        noOfStages_ += 3;   //to be sure. Can be deleted after testing.
        noOfStages_ += 4;   //to be even more sure. Can be deleted after testing.

        noOfStages_++;	//we need s+1 vectors for s stages. DONT DELETE THIS

        cout << "noOfStages: " << noOfStages_ << endl;
        //fill solution with empty lists
        solution.resize(noOfStages_);


        /*
        vector<int> firstCycle(bh_->bits.size());
        for(unsigned w = 0; w < bh_->bits.size(); w++){
            firstCycle[w] = bh_->bits[w].size();
        }
        newBits.push_back(firstCycle);
        originalNewBits.push_back(firstCycle);
        */

        for(unsigned w = 0; w < bh_->bits.size(); w++){
            unsigned i = 0;
            for(list<WeightedBit*>::iterator it = bh_->bits[w].begin(); it != bh_->bits[w].end(); it++){
                //cout << "bit " << i << " in column " << w << " has the cycle " << (*it)->getCycle() << endl;
                i++;
            }

        }

        //now fill the rest with zeros

        for(unsigned v = newBits.size(); v < (unsigned)noOfStages_ + 1; v++){
            vector<int> tempZeroVector(compOutputWordSizeMax * v + bh_->bits.size());
            for(unsigned w = 0; w < tempZeroVector.size(); w++){
                tempZeroVector[w] = 0;
            }
            newBits.push_back(tempZeroVector);
            originalNewBits.push_back(tempZeroVector);
        }

        for(unsigned i = 0; i < newBits.size(); i++){
            for(unsigned j = 0; j < newBits[i].size(); j++){
                cout << newBits[i][j] << " ";
            }
            cout << endl;
        }

        //cout << mode << endl;





        if(mode.compare("heuristic_parandeh-afshar_modified") == 0 || mode.compare("optimalMinStages") == 0){
            //this is an improved algorithm of the heuristic descriped by parandeh-afshar

            if(mode.compare("optimalMinStages") == 0){
                dontUsePA = true;
            }
            vector<BasicCompressor *> * compUnsorted = bh_->getPossibleCompressors();

/*
            //now add flipflop so that a lowerBound of 0 creates a valid solution
            BasicCompressor *flipflop;
            vector<int> col(1);
            col[0]=1;
            flipflop = new BasicCompressor(bh_->getOp()->getTarget(),col);
            flipflop->areaCost = 0.5;
            compUnsorted->push_back(flipflop);
*/
            if(paSorting){
                cout << "we sort the compressors according to the pa-efficiency: inputBits / outputBits" << endl;
            }
            else{
                cout << "we sort the compressors according to the efficiency : (inputBits - outputBits)/area" << endl;
            }
            //now sort them
            bool used[compUnsorted->size()];
            for(unsigned i = 0; i < compUnsorted->size(); i++){
                used[i] = false;
            }
            unsigned count = 0;
            unsigned coveringDots = 0;
            double efficiency = 0;
            while(compUnsorted->size() > count){
                efficiency = 0;
                coveringDots = 0;
                BasicCompressor * bc;
                unsigned pos;
                for(unsigned i = 0; i < compUnsorted->size(); i++){
                    if(!paSorting){ //our sorting - efficiency
                        if(computeEfficiency(compUnsorted->at(i)) > efficiency && !used[i]){
                            bc = compUnsorted->at(i);
                            efficiency = computeEfficiency(compUnsorted->at(i));
                            pos = i;
                        }
                    }
                    else{
                        unsigned tempCoveringDots = 0;
                        for(unsigned j = 0; j < compUnsorted->at(i)->height.size(); j++){
                            tempCoveringDots += compUnsorted->at(i)->height[j];
                        }

                        if(computeCompressionRatioPA(compUnsorted->at(i)) - efficiency > 0.001 && !used[i]){
                            //compression ratio is bigger
                            bc = compUnsorted->at(i);
                            efficiency = computeCompressionRatioPA(compUnsorted->at(i));
                            coveringDots = tempCoveringDots;
                            pos = i;
                        }
                        else if(fabs(computeCompressionRatioPA(compUnsorted->at(i)) - efficiency) < 0.001 && !used[i]){
                            //has compression ratio is the same, has it more covering dots?
                            if(tempCoveringDots > coveringDots){
                                bc = compUnsorted->at(i);
                                efficiency = computeCompressionRatioPA(compUnsorted->at(i));
                                coveringDots = tempCoveringDots;
                                pos = i;
                            }
                        }

                    }
                }

                bhCompressor tempCompressor;
                tempCompressor.pointer = bc;
                tempCompressor.originalPosition = pos;
                tempCompressor.maxEfficiency = efficiency;

                cout << "compressor " << tempCompressor.originalPosition << " " << *bc << " has the efficiency of ";
                cout << tempCompressor.maxEfficiency << "  and outputsize " << tempCompressor.pointer->getOutputSize() << endl;

                //now put it at the back of the compressors-vector
                compressors.push_back(tempCompressor);
                used[pos] = true;
                count++;
            }

            //now add flipflop so that a lowerBound of 0 creates a valid solution
            BasicCompressor *flipflop;
            vector<int> col(1);
            col[0]=1;
            flipflop = new BasicCompressor(bh_->getOp()->getTarget(),col);
            if(bh_->getOp()->getTarget()->isPipelined()){
                std::string targetID = bh_->getOp()->getTarget()->getID();
                if((targetID == "Virtex6") || (targetID == "Virtex7") || (targetID == "Spartan6")){
                    flipflop->areaCost = 0.5; //there are two flip-flops per LUT for Virtex 6/7 and Spartan6
                }
                else{
                    flipflop->areaCost = 1; //assume 1 for unknown device //!!!!
                }
            }
            else{
                flipflop->areaCost = 0.01; //nearly 0 for unpipelined designs
            }
            compUnsorted->push_back(flipflop);

            bhCompressor tempCompressor;
            tempCompressor.pointer = flipflop;
            tempCompressor.originalPosition = compressors.size();
            tempCompressor.maxEfficiency = 0.0;
            compressors.push_back(tempCompressor);


                  //debug compressors
            cout << "ordering of sorted compressor list:" << endl;
            for(unsigned i = 0; i < compressors.size(); i++){
                cout << "compressor " << i << " (original position " << compressors[i].originalPosition << "): " << *(compressors[i].pointer) << " has efficiency " << compressors[i].maxEfficiency << endl;
            }


            if(usePreReduction && (mode.compare("heuristic_parandeh-afshar_modified") != 0)){
                preReduction(-100, 1.75);
            }

            if((!passHeuristicSolution && !dontUsePA) || !dontUsePA){
                cout << "using algorithm" << endl;

                if(passMultipleSolutions){
#ifdef HAVE_SCIP
                    cout << "solutions in heuristicSolutions before first round are " << bitHeapILPCompression.heuristicSolutions.size() << endl;

                    generateAllHeu(false, 1.5);

                    bhCompressor tempCompressor = compressors[0];
                    compressors.erase(compressors.begin());
                    tempCompressor.maxEfficiency = -1.0;
                    compressors.push_back(tempCompressor);
                    cout << "solutions in heuristicSolutions are " << bitHeapILPCompression.heuristicSolutions.size() << endl;
                    generateAllHeu(false, 1.5);
#endif //HAVE_SCIP
                }
                else{
                    cout << "dontUsePA should be 0, equals " << dontUsePA << endl;
                    algorithm();            //dontUsePA == false
                }




            }
            else if(!dontUsePA || passMultipleSolutions){
#ifdef HAVE_SCIP
                cout << "!dontUsePA or passMultipleSolutions" << endl << endl;
                if(passMultipleSolutions){
                    cout << "solutions in heuristicSolutions before first round are " << bitHeapILPCompression.heuristicSolutions.size() << endl;

                    generateAllHeu(false, 1.5);

                    bhCompressor tempCompressor = compressors[0];
                    compressors.erase(compressors.begin());
                    tempCompressor.maxEfficiency = -1.0;
                    compressors.push_back(tempCompressor);
                    cout << "solutions in heuristicSolutions are " << bitHeapILPCompression.heuristicSolutions.size() << endl;
                    generateAllHeu(false, 1.5);

                }
                else if(!optimalAndPreReduction){
                    algorithm();
                    addHeuristicSolutionToILPSolver();
                }
                else{
                    //optimal, preReduction is done
                    cout << "preReduction done. We now use the optimal method" << endl;
					cout << "________________________________________________" << endl << endl << endl << endl;
                    printBitHeap();
                }
#endif //HAVE_SCIP

            }

            if(preSolution.size() > 0){
                //no merging just replacing because previous solution should be empty?
                solution = preSolution;
            }



        }
        else if(mode.compare("heuristic_pa") == 0){
            //real parandeh-afshar

            //setting the flags so that the ilp-solver does not "work"
            passHeuristicSolution = false;
            passMultipleSolutions = false;
            reduceILPStageCount = false;


            cout << "in heuristic_pa" << endl;

            if(paSorting){
                cout << "we sort the compressors according to the pa-efficiency: inputBits / outputBits" << endl;
            }
            else{
                cout << "we sort the compressors according to the efficiency : (inputBits - outputBits)/area" << endl;
            }

            vector<BasicCompressor *> * compUnsorted = bh_->getPossibleCompressors();
            //now sort them
            bool used[compUnsorted->size()];
            for(unsigned i = 0; i < compUnsorted->size(); i++){
                used[i] = false;
            }
            unsigned count = 0;
            double compressionRatio = 0;
            while(compUnsorted->size() > count){
                compressionRatio = 0;
                BasicCompressor * bc;
                unsigned pos = 0;
                for(unsigned i = 0; i < compUnsorted->size(); i++){
                    if(paSorting){
                        if(computeCompressionRatioPA(compUnsorted->at(i)) > (compressionRatio + 0.0001) && !used[i]){
                            bc = compUnsorted->at(i);
                            compressionRatio = computeCompressionRatioPA(compUnsorted->at(i));
                            pos = i;
                        }
                        else if((fabs(computeCompressionRatioPA(compUnsorted->at(i)) - compressionRatio) < 0.0001) && !used[i]){
                            //now check which one has more covering bits:

                            int oldCoveringDots = 0;
                            int newCoveringDots = 0;
                            for(unsigned l = 0; l < compUnsorted->at(pos)->height.size(); l++){
                                oldCoveringDots += compUnsorted->at(pos)->height[l];
                            }
                            for(unsigned l = 0; l < compUnsorted->at(i)->height.size(); l++){
                                newCoveringDots += compUnsorted->at(i)->height[l];
                            }

                            if(newCoveringDots > oldCoveringDots){
                                bc = compUnsorted->at(i);
                                compressionRatio = computeCompressionRatioPA(compUnsorted->at(i));
                                pos = i;
                            }

                        }
                    }
                    else{
                        if(computeEfficiency(compUnsorted->at(i)) > (compressionRatio + 0.0001) && !used[i]){
                            bc = compUnsorted->at(i);
                            compressionRatio = computeEfficiency(compUnsorted->at(i));
                            pos = i;
                        }
                    }

                }

                bhCompressor tempCompressor;
                tempCompressor.pointer = bc;
                tempCompressor.originalPosition = pos;
                tempCompressor.maxEfficiency = compressionRatio;

                cout << "compressor " << tempCompressor.originalPosition << " " << *bc << " has the efficiency of ";
                cout << tempCompressor.maxEfficiency << "  and outputsize " << tempCompressor.pointer->getOutputSize() << endl;

                //now put it at the back of the compressors-vector
                compressors.push_back(tempCompressor);
                used[pos] = true;
                count++;


            }

            //now add flipflop so that a lowerBound of 0 creates a valid solution
            BasicCompressor *flipflop;
            vector<int> col(1);
            col[0]=1;
            flipflop = new BasicCompressor(bh_->getOp()->getTarget(),col);
            if(bh_->getOp()->getTarget()->isPipelined()){
                std::string targetID = bh_->getOp()->getTarget()->getID();
                if((targetID == "Virtex6") || (targetID == "Virtex7") || (targetID == "Spartan6")){
                    flipflop->areaCost = 0.5; //there are two flip-flops per LUT for Virtex 6/7 and Spartan6
                }
                else{
                    flipflop->areaCost = 1; //assume 1 for unknown device //!!!!
                }
            }
            else{
                flipflop->areaCost = 0.01; //nearly 0 for unpipelined designs
            }
            compUnsorted->push_back(flipflop);

            bhCompressor tempCompressor;
            tempCompressor.pointer = flipflop;
            tempCompressor.originalPosition = compressors.size();
            tempCompressor.maxEfficiency = 0.0;
            compressors.push_back(tempCompressor);
#ifdef HAVE_SCIP
            bitHeapILPCompression.dontAddFlipFlop = true;
#endif //HAVE_SCIP
            //sorting is the same for the moment TODO: change sorting to parandeh-afshar (sort for bitreduction)

            //debug compressors
            //debug compressors
            for(unsigned i = 0; i < compressors.size(); i++){
                cout << "(";
                for(unsigned j = 0; j < compressors[i].pointer->getNumberOfColumns(); j++){
                    cout << compressors[i].pointer->getColumnSize(j) << ",";
                }

                cout << ";";
                cout << compressors[i].pointer->getOutputSize() << ") with efficiency " << compressors[i].maxEfficiency << endl;
            }



            parandehAfshar();


        }

        //debug output:: newBits
        cout << "newbits:" << endl;
        for(unsigned i = 0; i < newBits.size(); i++){
            for(unsigned j = 0; j < newBits[i].size(); j++){

                cout << newBits[i][j] << " ";
            }
            cout << endl;
        }

        if((mode.compare("heuristic_parandeh-afshar_modified") == 0) || (mode.compare("heuristic_pa") == 0)){
            cout << "size of solution is " << solution.size() << endl;
            return 0;
        }
/*
        if((mode.compare("heuristic_pa") == 0)){
            cout << "size of solution is " << solution.size() << endl;
            return 0;
        }
*/
        //debug output:: compressors:
        /*
        for(unsigned s=0; s < solution.size(); s++)
        {
            list<pair<int,int> >::iterator it;
            for(it = solution[s].begin(); it != solution[s].end(); ++it)
            {
                cout << "applying compressor " << (*it).first << " to column " << (*it).second << " in stage " << s << endl;

            }
        }
        */

#ifdef HAVE_SCIP

        //now initialize BitHeapILPCompression
        bitHeapILPCompression.solution = solution;
        bitHeapILPCompression.newBits = newBits;
        bitHeapILPCompression.useHeuristic = true;
        bitHeapILPCompression.zeroStages = 0;               //generally assume that there are none presolved stages
        if(reduceILPStageCount){
            setUpNewStageCountOfILP();
        }
        if(solution.size() > 0){
            //compute the area size of the already chosen compressors for better comparability
            double area = computeAreaofSolution();
            bitHeapILPCompression.preReductionAreaCost = area;
        }

        if(differentStages){
            preSolution.clear();
            preSolution = solution;         //saving the solution, maybe we did a reduction (e.g. buildSingleStages)
            solution.clear();
            cout << "noOfStages_: " << noOfStages_ << endl;
            solution.resize(noOfStages_);

            if(!useCompleteHeuristic){
                //check if the heuristic generated a solution already
                //we assume a finaladderheight of 2
                bool bitFound = false;
                bool finished = true;
                bool zerosSet = false;
  //              unsigned zeros;     //not needed because it doesn't make it faster
                for(unsigned s = 0; s < newBits.size() && finished; s++){
                    //if bitFound == true : a previous stage has unfinished bits
                    //if we find now bits as well, we are not ready

                    //if we found bits in the previous stage as well in the current Stage, we are not finished
                    bool previousStageBitsFound = bitFound;

                    for(unsigned c = 0; c < newBits[s].size(); c++){
                        if(newBits[s][c] > 2){
                            //break
                            minFixedStage = s;
                            if(minFixedStage == 0){
                                minFixedStage++;
                            }
                            finished = false;
                            break;
                        }
                        else if(newBits[s][c] > 0 && newBits[s][c] < 3){
                            if(previousStageBitsFound){
                                minFixedStage = s;
                                finished = false;
                                break;
                            }
                            bitFound = true;
                        }
                    }

                    if(!zerosSet && bitFound){  //set zeros at previous stage and lock that value
                        if(s > 0){
//                            zeros = s - 1;
                            zerosSet = true;
                        }
                    }
                    if(!zerosSet && !bitFound && !finished){
                        //we found a stage with more than 2 bits, but no bits in previous stage
                        if(s > 0){
//                            zeros = s - 1;
                            zerosSet = true;
                        }
                    }
                }
                cout << "minFixedStage is now " << minFixedStage << endl;

                if(finished){
                    //heuristic found complete solution therefore set flag useCompleteHeuristic
                    cout << "useCompleteHeuristic set to TRUE" << endl;
                    useCompleteHeuristic = true;
                }

            }

            if(useCompleteHeuristic){
                int stagesUsedByHeuristic = 0;
                //check the size of the preSolution:
                for(unsigned i = 0; i < preSolution.size(); i++){
                    if(preSolution[i].size() > 0){
                        stagesUsedByHeuristic++;
                    }
                    else{
                        break;
                    }
                }

                minFixedStage = stagesUsedByHeuristic;
                maxFixedStage = stagesUsedByHeuristic;
                cout << "set minFixedStage and maxFixedStage equal" << endl;
                bitHeapILPCompression.zeroStages = 0; //trivial - therefore problemreduction is not necessary (1 second)
                setUpILPForMoreStages(minFixedStage, true);
            }
            setUpILPForMoreStages(minFixedStage, true);


            /*
            if(buildSingleStages){
                bitHeapILPCompression.zeroStages = numberOfBuildStages;
            }
            else{
                bitHeapILPCompression.zeroStages = 0;
            }
            bitHeapILPCompression.getExternalStageCount = true;
            bitHeapILPCompression.noOfStages_ = minFixedStage;
            bitHeapILPCompression.useFixedStageCount = true;
            bitHeapILPCompression.solution.clear();
            bitHeapILPCompression.solution.resize(noOfStages_);

            */
        }


        cout << "initialisation of bitHeapILPCompression finished" << endl;
        bitHeapILPCompression.generateProblem();

        if(passHeuristicSolution){
            int count = bitHeapILPCompression.passHeuristicSolutions();

            cout << "there were " << count << " heuristic solutions passed" << endl;
        }
#endif //HAVE_SCIP
        return 0;
    }


    void BitHeapHeuristicCompression::algorithm(){
        cout << "in algorithm__________________" << endl;
        cout << "compresssors.size() = " << compressors.size() << endl;


		printLowerBounds();
		cout << "___________" << endl;

        //if debugLoop == true, there are only debugMax steps
        bool debugLoop = false;
        int debugMax = 50;
        int debugCount = 0;

        unsigned numberOfStages = newBits.size() - 1;

        if(buildSingleStages && (numberOfBuildStages > 0) ){
            numberOfStages = numberOfBuildStages;
        }
		bool compressionNeeded = false;
		for(unsigned b = 0; b < newBits[0].size(); b++){
			if(newBits[0][b] > 2){
				compressionNeeded = true;
				break;
			}
			else{
				for(unsigned a = 1; a < newBits.size(); a++){
					if(newBits[a][b] > 0){
						compressionNeeded = true;
					}
				}
			}
		}
		if(compressionNeeded == false){
			cout << "no compression necessary" << endl;
			return;
		}

        //find a compressor which fits the best and use it.
        for(unsigned s = 0; s < numberOfStages && !(debugLoop && debugCount >= debugMax); s++){
            debugCount++;
            bool found = true;
			bool foundVariableCompressor = false;

            while(found == true && !(debugLoop && debugCount >= debugMax)){
                debugCount++;

                found = false;
				foundVariableCompressor = false;

                double achievedEfficiencyBest = -1.0;
                unsigned compressor = 0;
				unsigned variableCompressor = 0;
                unsigned column = 0;
				unsigned variableCompressorMidLength = 0;


                if(useVariableCompressors == true){
					for(unsigned tempVariableLowCompressor = 0; tempVariableLowCompressor < variableBCompressors.size(); tempVariableLowCompressor += 3){
						pair<double, pair<unsigned, unsigned> > variableResult = variableCompEffBitHeap(s, tempVariableLowCompressor);
						double variableAchievedEfficiency = variableResult.first;
						unsigned tempColumn = variableResult.second.first;
						unsigned tempVariableCompressorMidLength = variableResult.second.second;

						//cout << "for s = " << s << " c = " << tempColumn <<  ", lowCompressor = " << tempVariableLowCompressor << " and midLength = " << tempVariableCompressorMidLength << endl;
						//cout << "variableAchievedEfficiency is " << variableAchievedEfficiency << endl;

						if(variableAchievedEfficiency > achievedEfficiencyBest + 0.0001){
							bool necessary = variableCompressorNecessary(s, tempColumn, tempVariableCompressorMidLength, tempVariableLowCompressor);

							if(fabs(lowerBounds[s]) < 0.001){
								necessary = true;
							}
							if(necessary && variableAchievedEfficiency > (lowerBounds[s] - 0.0001)){
								//cout << "    variableAchievedEfficiency = " << variableAchievedEfficiency << "and lowerBounds[s] = " << lowerBounds[s] << endl;

								achievedEfficiencyBest = variableAchievedEfficiency;
								column = tempColumn;
								variableCompressorMidLength = tempVariableCompressorMidLength;
								variableCompressor = tempVariableLowCompressor;

								foundVariableCompressor = true;
							}


						}
					}
				}


                for(unsigned i = 0; i < compressors.size()  && !(debugLoop && debugCount >= debugMax); i++){
                    debugCount++;
                    //first get maximal efficiency of this compressor.
                    //if its lower than the efficiency already found, then stop
                    double maxEfficiency = compressors[i].maxEfficiency;
                    if((maxEfficiency - 0.0001) < achievedEfficiencyBest || (maxEfficiency + 0.0001) < lowerBounds[s]){    //accuracy
                        break;
                    }

                    //unsigned maxWidth = newBits[s].size() - (compressors[i].pointer->getOutputSize() - 1);
					unsigned maxWidth = newBits[s].size();
                    bool used[maxWidth];        //used to check whether we already tried the compressor in this column
                    for(unsigned k = 0; k < maxWidth; k++){
                        used[k] = false;
                    }

                    unsigned columnsAlreadyChecked = 0;
                    bool foundCompressor = false;
                    while((columnsAlreadyChecked < maxWidth) && !(foundCompressor && (fabs(maxEfficiency - achievedEfficiencyBest) < 0.0001)) ){ //maxEfficiency accuracy
                        //get max stage
                        unsigned currentMaxColumn = 0;
                        int currentSize = 0;
                        for(unsigned c = 0; c < maxWidth; c++){
                            //cout << "s = " << s << " c = " << c << " and newBits = " << newBits[s][c] << "  ---old- maxWidth = " << currentSize << " currentMaxColumn = " << currentMaxColumn<< endl;
                            if(!used[c] && newBits[s][c] > currentSize){
                                currentMaxColumn = c;
                                currentSize = newBits[s][c];
                            }
                        }
                        used[currentMaxColumn] = true;
                        //cout << "final maxColumn = " << currentMaxColumn <<" with width " << currentSize << " for compressor " << i << endl;

                        double achievedEfficiencyCurrent = compEffBitHeap(s, currentMaxColumn, i);
                        //cout << "achieved efficiency is " << achievedEfficiencyCurrent << endl;

                        //check if it is necessary
                        bool necessary = true;
                        if(!fabs(lowerBounds[s]) < 0.001){
                            necessary = false;
                            for(unsigned k = 0; k < currentMaxColumn + compressors[i].pointer->getNumberOfColumns(); k++){
                                unsigned sum = 0;
                                for(unsigned l = s; l < newBits.size(); l++){
                                    sum += newBits[l][k];
                                }
                                if(sum > 2){
                                    necessary = true;
                                    break;
                                }
                            }
                        }
                        if(achievedEfficiencyCurrent > (achievedEfficiencyBest + 0.0001) && achievedEfficiencyCurrent > (lowerBounds[s] - 0.0001) && necessary){       //accuracy !!
                            achievedEfficiencyBest = achievedEfficiencyCurrent;
                            compressor = i;
                            column = currentMaxColumn;
                            foundCompressor = true;
                            found = true;
							foundVariableCompressor = false;
                        }
                        columnsAlreadyChecked++;
                    }
                }

                if(found){
                    //we found a compressor and use it now.

                    cout << "-- using compressor " <<  compressor  << " in stage " << s << " and in column ";
					cout << column << " with efficiency " << achievedEfficiencyBest << endl;

                    useCompressor(s, column, compressor);
                    //printBitHeap();

                }
				if(foundVariableCompressor){
					//printBitHeap();
                    cout << "-- using variable compressor " << variableCompressor / 3 << " in stage " << s << " and in column ";	//the variable "variableCompressor" is the position of the variable basic compressor
                    cout << column << " and midLength " << variableCompressorMidLength << " with efficiency " << achievedEfficiencyBest << endl;
					useVariableCompressor(s, column, variableCompressorMidLength, variableCompressor);

					//printBitHeap();
					found = true;
				}

            }


            //now we are done with this stage. if we need to generate a solution, then we
            //have to put the remaining bits into the next stage.
            //if we reached the stage were we can use the two-bit-adder, we are done
            if(fabs(lowerBounds[s]) < 0.0001 ){
/*              //not needed because in the compressorlist is already a flipflop at the back.
                cout << "generate solution - therefore put remaining bits into the next stage" << endl;
                for(unsigned j = 0; j < newBits[s].size(); j++){
                    if(newBits[s][j] > 0){
                        int temp = newBits[s][j];
                        newBits[s + 1][j] += newBits[s][j];
                        newBits[s][j] = 0;

                        for(int k = 0; k < temp; k++){
                            solution[s].push_back(pair<int,int>((compressors.size() - 1),j));

                        }

                    }
                }
*/

                //check exit-condition

                bool twoBitAdderReached = true;
                for(unsigned a = 0; a < newBits.size(); a++){
                    if(a != s + 1){
                        //check if those stages are empty
                        for(unsigned b = 0; b < newBits[a].size(); b++){
                            if(newBits[a][b] > 0){
                                twoBitAdderReached = false;
                                break;
                            }
                        }
                    }
                    else{
                        //check if s + 1 stage is the final one (only has a maximum bitcount of 2
                        for(unsigned b = 0; b < newBits[a].size(); b++){
                            if(newBits[a][b] > 2){
                                twoBitAdderReached = false;
                                break;
                            }
                        }


                    }
                }

                if(twoBitAdderReached){
                    //cout << endl << "twoBitAdderStage reached" << endl;
                    //printBitHeap();
                    break;
                }
            }
        }

        /*
        //workaround
        if(passHeuristicSolution){
            bitHeapILPCompression.heuristicSolutions.push_back(solution);
            solution.clear();
            solution.resize(noOfStages_);
            cout << endl << "bitHeap number " << (bitHeapILPCompression.heuristicSolutions.size() - 1 ) << endl;
            printBitHeap();

            newBits = originalNewBits;
        }

        */



        if(useHoles){
            //now delete the holes in the first stages, because they can't be used there
            for(unsigned i = 0; i < newBits[0].size(); i++){
                if(newBits[0][i] < 0){
                    newBits[0][i] = 0;
                }
            }

            //cout << "newBits after deleting holes in stage 0" << endl;
            //printBitHeap();
        }





    }

    pair<int,int> BitHeapHeuristicCompression::parandehAfsharSearch(int stage, int column){

        int compressor = -1;
        int resultColumn = -1;
        double achievedEfficiencyBest = -1;
        bool found = false;



        for(unsigned i = 0; i < compressors.size(); i++){
            double achievedEfficiencyCurrentLeft = compEffBitHeap(stage, column, i);

            double achievedEfficiencyCurrentRight = -1;
            int rightStartPoint = column - (compressors[i].pointer->getNumberOfColumns() - 1);
            //cout << "rightStartPoint is " << rightStartPoint << " and length of compressorColumns -1 is " << (compressors[i].pointer->getNumberOfColumns() - 1) << endl;
            if(rightStartPoint >= 0){
                achievedEfficiencyCurrentRight = compEffBitHeap(stage, (unsigned) rightStartPoint, i);
            }
            //cout << "column is " << currentMaxColumn << " and compressor is " << i << endl;
            //cout << " right Eff = " << achievedEfficiencyCurrentRight << " and left Eff = " << achievedEfficiencyCurrentLeft << endl;

            if(achievedEfficiencyCurrentLeft > 0.0001 && achievedEfficiencyCurrentRight <= achievedEfficiencyCurrentLeft + 0.0001){
                //normal (left) search is successful - prefer it if left and right search has equal efficiency

                if(achievedEfficiencyBest + 0.0001 < achievedEfficiencyCurrentLeft){
                    achievedEfficiencyBest = achievedEfficiencyCurrentLeft;
                    compressor = i;
                    resultColumn = column;
                }
                found = true;
            }
            else if(achievedEfficiencyCurrentRight > 0.0001){
                //right search is successful

                if(achievedEfficiencyBest + 0.0001 < achievedEfficiencyCurrentRight){
                    achievedEfficiencyBest = achievedEfficiencyCurrentRight;
                    compressor = i;
                    resultColumn = rightStartPoint;
                }
                found = true;
            }
            //cout << "current best compressor = " << compressor << " and column  = " << column <<  " with efficiency " << achievedEfficiencyBest << endl;
        }
        //cout << "after pa search achievedEfficiencyBest is " << achievedEfficiencyBest << " and column " << resultColumn << " and compressor is " << compressor << endl;

        pair<int, int> result;
        if(found == true){
            result.first = resultColumn;
            result.second = compressor;
        }
        else{
            result.first = -1;
            result.second = -1;
        }
        return result;
    }

    void BitHeapHeuristicCompression::parandehAfshar(){

        cout << "in parandeh-Afshar" << endl;
		double threshold = 0.0;
		if(paSorting){
			cout << "paSorting is true" << endl;
			threshold += 1.0;		//paSorting: to reduce bits, the metric of a compressor has to be bigger than one  (metric is here: inputbits / outputbits)
		}
		else{
			cout << "paSorting is false" << endl;
		}

        bool exit = false;

        for(unsigned s = 0; s < (newBits.size() - 1); s++){
            //cout << "stage = " << s << endl;
            bool found = true;
            while(found){
                found = false;
                unsigned compressor = 0;
                unsigned column = 0;
                pair<int,int> result;
                //cout << "now find compressor" << endl;
                unsigned currentMaxColumn = 0;
                int maxSize = 0;
#if 0
                for(unsigned c = 0; c < newBits[s].size(); c++){
                    //cout << "s = " << s << " c = " << c << " and newBits = " << newBits[s][c] << "  ---old- maxWidth = " << maxSize << " currentMaxColumn = " << currentMaxColumn<< endl;
                    if(newBits[s][c] > maxSize){
                        currentMaxColumn = c;
                        maxSize = newBits[s][c];
                    }
                }
                //cout << "maxColumn is " << currentMaxColumn << " and maxSize is " << maxSize << endl;
                if(maxSize > 0){ //otherwise there are no bits left in this stage
                    result = parandehAfsharSearch(s, currentMaxColumn);
                    if(result.first >= 0){
                        found = true;
                    }
                }
#endif

                bool used[newBits[s].size()];
                for(unsigned k = 0; k < newBits[s].size(); k++){
                    used[k] = 0;
                }
                for(unsigned a = 0; a < newBits[s].size(); a++){
                    //find maxColumn
                    for(unsigned c = 0; c < newBits[s].size(); c++){
                        if(used[c] == false && newBits[s][c] > maxSize){
                            currentMaxColumn = c;
                            maxSize = newBits[s][c];
                        }
                    }
                    used[currentMaxColumn] = true;

                    //check found column
                    if(maxSize > 0){
                        result = parandehAfsharSearch(s, currentMaxColumn);
                        if(result.first >= 0){
                            found = true;
                        }
                    }
                    if(found){
                        break;
                    }
                }


                //now use the compressor
                if(found){
                    column = result.first;
                    compressor = result.second;
                    useCompressor(s, column, compressor);
                    cout << "used compressor "<< compressor << " in column " << column << " and stage " << s << endl;
                    //printBitHeap();
                    //cout << "new iteration" << endl;
                    //cout << endl;
                }
            }

            //now we are done with this stage. if we need to generate a solution, then we
            //have to put the remaining bits into the next stage.
            //if we reached the stage were we can use the two-bit-adder, we are done
            if(generateSolution){
                cout << "generate solution - therefore put remaining bits into the next stage" << endl;
                for(unsigned j = 0; j < newBits[s].size(); j++){
                    if(newBits[s][j] > 0){
                        int temp = newBits[s][j];
                        newBits[s + 1][j] += newBits[s][j];
                        newBits[s][j] = 0;

                        //assuming: flip-flop is not in compressors, but after that at the back of bh_->possibleCompressors
                        //which is compressors.size
                        for(int k = 0; k < temp; k++){
                            solution[s].push_back(pair<int,int>((compressors.size() - 1),j));

                        }

                    }
                }

                //check exit-condition

                bool twoBitAdderReached = true;
                for(unsigned a = 0; a < newBits.size(); a++){
                    if(a != s + 1){
                        //check if those stages are empty
                        for(unsigned b = 0; b < newBits[a].size(); b++){
                            if(newBits[a][b] > 0){
                                twoBitAdderReached = false;
                                break;
                            }
                        }
                    }
                    else{
                        //check if s + 1 stage is the final one (only has a maximum bitcount of 2
                        for(unsigned b = 0; b < newBits[a].size(); b++){
                            if(newBits[a][b] > 2){
                                twoBitAdderReached = false;
                                break;
                            }
                        }


                    }
                }

                if(twoBitAdderReached){
                    //cout << endl << "twoBitAdderStage reached" << endl;
                    //printBitHeap();
                    break;
                }

                if(!twoBitAdderReached && !found && s == newBits.size() - 1){
                    //add new stage;
                    vector<int> zeroVector;
                    zeroVector.resize(newBits[s].size() + compOutputWordSizeMax);
                    newBits.push_back(zeroVector);

                    cout << "added a new stage because compression tree generation isn't finished yet" << endl;
                }


            }

            cout << "after stage " << s << " the bitheap looks like" << endl;
            printBitHeap();

            //if we reached the last stage, and there are still more than 3 bits in any column of the last stage, add an additional stage
            if(s == newBits.size() - 2){
                cout << "last stage" << endl;
                bool addNewStage = false;
                for(unsigned k = 0; k < newBits[s].size(); k++){
                    if(newBits[s + 1][k] > 2){
                        addNewStage = true;
                        break;
                    }
                }

                if(addNewStage){
                    cout << "we weren't able to compute a solution within " << s << " stages. Therefore we need an additional stage." << endl;
                    vector<int> tempList;
                    tempList.resize(newBits[s].size() + compOutputWordSizeMax);
                    for(unsigned k = 0; k < tempList.size(); k++){
                        tempList[k] = 0;
                    }
                    cout << "old newBitsSize  = " << newBits.size();
                    cout << "; noOfStages_ is << " << noOfStages_ << " old" << endl;
                    newBits.push_back(tempList);
                    solution.resize(newBits.size());
                    noOfStages_++;
#ifdef HAVE_SCIP
                    bitHeapILPCompression.solution.resize(newBits.size());
                    bitHeapILPCompression.noOfStages_ = noOfStages_;
                    bitHeapILPCompression.getExternalStageCount = true;
#endif //HAVE_SCIP
                    cout << "NEW newBitsSize is " << newBits.size() << " and noOfStages_ is " << noOfStages_ << endl;
                }
            }




        }



    }


    //this reduces the problem so that a smaller problem is passed to the ilp-solver
    void BitHeapHeuristicCompression::preReduction(unsigned minHeight, double minEff){
        //minimumHeightReduction = minHeight;
        //minimumEffReduction = minEff;
        //lowerBound = minEff;
        cout << "in prereduction " << endl << endl << endl;

        //disabling generate solution by not putting bits into the next stage with ffs
        bool rememberGenerateSolution = false;
        if(generateSolution){
            generateSolution = false;
            rememberGenerateSolution = true;

        }

        //the real checking if the conditions are met is done in compEffBitHeap()
        algorithm();

        cout << "preReduction done. The problem is now reduced to the following bitheap:" << endl;
        printBitHeap();

        //enabling generateSolution
        if(rememberGenerateSolution){
            generateSolution = true;
        }

        originalNewBits = newBits;
        preSolution = solution;
        solution.clear();
        solution.resize(noOfStages_);
        usePreReduction = false;



    }


    //this function places and (6,0,6;5) compressor in every column starting in column 0.
    //this leads to a more smoother outputBits distribution of 5 in every column (except the start & end). On the other hand the inputbits are 6,12,12,12,12,.... 12,6
    //if it reaches the maximum bits of the bitheap and there is at least one compressor placed, it will start at column 0 again.
    //it will be only placed, if the efficiency is reached.
    bool BitHeapHeuristicCompression::smootherCompressorPlacing(unsigned s, double efficiency){
        //assumption: (6,0,6:5) is first compressor
        bool atLeastOnceUsed = false;

        bool found = true;

        while(found){
            found = false;
            for(unsigned c = 0; c < newBits[s].size() - 2; c++){
                if((compEffBitHeap(s, c, 0) + 0.00001) >= efficiency){
                    found = true;
                    atLeastOnceUsed = true;
                    useCompressor(s, c, 0);
                }
            }
        }

        return atLeastOnceUsed;
    }


    //this function generates all Heuristic solutions by switching the priority of compressors with the same efficiency
    //if allCombinations == true, then it will compute all possible rankings and generate solutions
    //otherwise it will only switch the compressors with the given parameter efficiency
    int BitHeapHeuristicCompression::generateAllHeu(bool allCombinations, double efficiency){

        if(!allCombinations){
            unsigned count = 0;
            unsigned firstOccurance = 0;
            bool foundCompressor = false;
            for(unsigned i = 0; i < compressors.size(); i++){
                if(fabs(compressors[i].maxEfficiency - efficiency) < 0.001){
                    count++;
                    if(!foundCompressor){
                        firstOccurance = i;
                    }
                    foundCompressor = true;

                }



            }
            cout << count << " found for efficiency " << efficiency << endl;
            cout << "first occurance is " << firstOccurance << endl;

            //storing compressors and originalCompressorPosition;
            vector<bhCompressor> unmodifiedCompressors;
            //vector<unsigned> unmodifiedPositions;


            for(unsigned i = 0; i < compressors.size(); i++){
                unmodifiedCompressors.push_back(compressors[i]);
            }


            int combinations[count];
            for(unsigned i = 0; i < count; i++){
                combinations[i] = firstOccurance + i;
            }
            int times = 0;

            do{
                algorithm();
                addHeuristicSolutionToILPSolver();
                times++;

                for(unsigned i = 0; i < count; i++){
                    compressors[firstOccurance + i] = unmodifiedCompressors[combinations[i]];
                    //originalCompressorPosition[firstOccurance + i] = unmodifiedPositions[combinations[i]];
                }

            } while (std::next_permutation(combinations, combinations + count));


            //cout << "we were " << times <<" times in the do-while loop" << endl;
        }
        else{
            algorithm();
            addHeuristicSolutionToILPSolver();
        }




        return 0;
    }

    //this adds the solution to the heuristicSolutions in the bitHeapILPCompression
    //and resets the problem
    void BitHeapHeuristicCompression::addHeuristicSolutionToILPSolver(){
#ifdef HAVE_SCIP
        bitHeapILPCompression.heuristicSolutions.push_back(solution);
        solution.clear();
        solution.resize(noOfStages_);
        cout << endl << "bitHeap number " << (bitHeapILPCompression.heuristicSolutions.size() - 1 ) << endl;
        printBitHeap();

        newBits = originalNewBits;
#endif //HAVE_SCIP
    }



    //this is the standard efficiency E = delta / k; delta is Bits which are reduced and k is the size of the compressor in Luts
    double BitHeapHeuristicCompression::computeEfficiency(BasicCompressor* comp){
        int inputBits = 0;
        int outputBits = 0;

        for(unsigned i = 0; i < comp->height.size(); i++){
            inputBits += comp->height[i];
        }
        for(unsigned i = 0; i < comp->outputs.size(); i++){
            outputBits += comp->outputs[i];
        }

        int reducedBits = inputBits - outputBits;

        return  (double) reducedBits / comp->areaCost;
    }



    //ranking of Parandeh-Afshar descriped in his paper "Efficient Synthesis of Compressor Trees 2008"
    double BitHeapHeuristicCompression::computeCompressionRatioPA(BasicCompressor* comp){
        int inputBits = 0;
        int outputBits = 0;

        for(unsigned i = 0; i < comp->height.size(); i++){
            inputBits += comp->height[i];
        }
        for(unsigned i = 0; i < comp->outputs.size(); i++){
            outputBits += comp->outputs[i];
        }


        return (double) inputBits / outputBits;
    }


	//checks if a variable compressor is necessary. that is the case if in one of its inputcolumns are more than two bits to compress
	bool BitHeapHeuristicCompression::variableCompressorNecessary(unsigned s, unsigned column, unsigned midLength, unsigned compType){

		//low part
		unsigned sum = 0;
		for(unsigned l = s; l < newBits.size(); l++){
			if(newBits[l].size() > column){
				sum += newBits[l][column];
			}
		}
		if(sum > 2){
			return true;
		}

		//mid part
		for(unsigned m = 0; m < midLength; m++){
			sum = 0;
			for(unsigned l = s; l < newBits.size(); l++){
				if(newBits[l].size() > column + m){
					sum += newBits[l][column + m];
				}
			}
			if(sum > 2){
				return true;
			}
		}

		//high part
		//no need to consider that heights of high compressor are reversed
		for(unsigned highLength = 0; highLength < variableBCompressors[compType + 2].height.size(); highLength++){
			sum = 0;
			for(unsigned l = s; l < newBits.size(); l++){
				if(newBits[l].size() > column + midLength + highLength){
					sum += newBits[l][column + midLength + highLength];
				}
			}
			if(sum > 2){
				return true;
			}
		}

		//compressor not necessery
		return false;
	}



	//returns the efficiency (as double), as well as starting column and length of middle Part of the variable Compressor compPos with the highest efficiency. if efficiency = -10, then we didn't find a suitable compressor
	pair<double, pair<unsigned, unsigned> > BitHeapHeuristicCompression::variableCompEffBitHeap(unsigned s, unsigned compType){
		//because we are before the normal compressors, we also have to go through all the starting columns
		//cout << "in variableCompEffBitHeap" << endl;
		unsigned column = 0;
		unsigned middleLength = 0;
		double achievedEfficiencyBest = -10.0;

		unsigned maxWidth = newBits[s].size();
		bool used[maxWidth];
		for(unsigned k = 0; k < maxWidth; k++){
			used[k] = false;
		}

		unsigned columnsAlreadyChecked = 0;
		//cout << " setup done" << endl;
		while(columnsAlreadyChecked < maxWidth){

			unsigned currentMaxColumn = 0;
			int currentSize = 0;
			for(unsigned c = 0; c < maxWidth; c++){		//searching for max Column
				if(!used[c] && newBits[s][c] > currentSize){
					currentMaxColumn = c;
					currentSize = newBits[s][c];
				}
			}
			used[currentMaxColumn] = true;

			//cout << "variableBitEfficiency: currentMaxColumn is " << currentMaxColumn << endl;

			for(unsigned tempMiddleLength = 0; currentMaxColumn + tempMiddleLength < newBits[s].size(); tempMiddleLength++){
				double achievedEfficiencyCurrent = variableCompEffBitHeapBasic(s, currentMaxColumn, tempMiddleLength, compType);
				//cout << "for middleLength of " << tempMiddleLength << " is the efficiency " << achievedEfficiencyCurrent << endl;
				if(achievedEfficiencyCurrent > achievedEfficiencyBest + 0.00001){
					//cout << "updated achievedEfficiency from " << achievedEfficiencyBest << " to " << achievedEfficiencyCurrent << endl;
					column = currentMaxColumn;
					middleLength = tempMiddleLength;
					achievedEfficiencyBest = achievedEfficiencyCurrent;
				}
			}
			columnsAlreadyChecked++;
		}


		pair<unsigned, unsigned> tempResult (column, middleLength);

		pair<double, pair<unsigned, unsigned> > result (achievedEfficiencyBest,tempResult);

		return result;
	}


	double BitHeapHeuristicCompression::variableCompEffBitHeapBasic(unsigned s, unsigned c, unsigned midLength, unsigned compType){
		int sum = 0;
		//cout << "compType is " << compType << endl;
		double cost = variableBCompressors[compType].areaCost + midLength * variableBCompressors[compType + 1].areaCost + variableBCompressors[compType + 2].areaCost;
		int outputSize = variableBCompressors[compType].outputs[0] + midLength * variableBCompressors[compType + 1].outputs[0];
		for(unsigned i = 0; i < variableBCompressors[compType + 2].outputs.size(); i++){
			outputSize += variableBCompressors[compType + 2].outputs[i];
		}

		//cout << "the sum of the outputs is " << outputSize << endl;

		unsigned maxSize = newBits[s].size();
		//computing sum
		if(maxSize > c){
			if(newBits[s][c] > 0){

				if((unsigned) newBits[s][c] >= (unsigned) variableBCompressors[compType].height[0]){
					sum += variableBCompressors[compType].height[0];
				}
				else{
					sum += newBits[s][c];
				}
			}
		}

		//middle part
		for(unsigned m = 0; m < midLength; m++){
			if(maxSize > c + m + 1){
				if(newBits[s][c + m + 1] > 0){
					if((unsigned) newBits[s][c + m + 1] >= (unsigned) variableBCompressors[compType + 1].height[0]){
						sum += variableBCompressors[compType + 1].height[0];
					}
					else{
						sum += newBits[s][c + m + 1];
					}
				}
			}
		}
		//high part
		unsigned bitPosition = c + midLength + 1;	//position in newBits
		//multicolumns: inputs and outputs reversed
		for(int compPosition = variableBCompressors[compType + 2].height.size() - 1; compPosition >= 0; compPosition--){
			//cout << "compPosition is " << compPosition << endl;
			if(maxSize > bitPosition){
				if(newBits[s][bitPosition] > 0){

					if((unsigned) newBits[s][bitPosition] >= (unsigned) variableBCompressors[compType + 2].height[compPosition]){
						sum += variableBCompressors[compType + 2].height[compPosition];
					}
					else{
						sum += newBits[s][bitPosition];
					}
				}
			}


			bitPosition++;
		}

		//cout << "the sum of used inputs is " << sum << endl;
		int bitsReduced = ((int) sum) - ((int) outputSize);
		//cout << "bits Reduced is " << bitsReduced << " and the cost is " << cost << endl;
		double eff = ((double) bitsReduced) / cost;
		return eff;
	}





    //conditions met disabled!!

    //this is the standardEfficiency E you would achieve if the compressor is used in stage s column c
    //the conditions for compressors from parandeh-afshar are checked as well
    double BitHeapHeuristicCompression::compEffBitHeap(unsigned s, unsigned c, unsigned compPos){
        int sum = 0;
        bool conditionsMet = true;
//        bool preReductionConditionsMet = true;

        //make sure the width of the current stage is large enough
        if(newBits[s].size() - 1 < c + (compressors[compPos].pointer->getNumberOfColumns() - 1)){
            //cout << "we got over the edge of the current stage " << s << " at column " << c << " and original compPos " << compressors[compPos].originalPosition << endl;
        }

        for(unsigned i = 0; i < compressors[compPos].pointer->getNumberOfColumns(); i++){

            if(i == 0 && newBits[s][c] == 1){
                //the input of the last bit is only connected to the lowest output
                conditionsMet = true;           //TODO: check whether we need conditionsMet
            }

            if(i + c < newBits[s].size()){  //we only add something to "sum" if there is a column

                if(newBits[s][c + i] > 0){
                    //we do not count "holes"
                    if((unsigned) newBits[s][c + i] >= compressors[compPos].pointer->getColumnSize(i)){
                        sum += compressors[compPos].pointer->getColumnSize(i);
                    }
                    else{
                        sum += newBits[s][c + i];
                    }
                }

/*
                if(((unsigned) newBits[s][c + i]) - compressors[compPos].pointer->getColumnSize(i) < (int) minimumHeightReduction){
                    preReductionConditionsMet = false;      //deactivated
                }
*/
            }


        }

		if(paSorting){
			double ratio = ((double) sum) / ((double) compressors[compPos].pointer->getOutputSize());
			return ratio;
		}

        //now in sum are the bits which are covered with this compressor
        sum -= compressors[compPos].pointer->getOutputSize();

        if(compressors[compPos].pointer->areaCost < 0.00001){
            cout << " compressor " << compPos << " has areaCost of " << compressors[compPos].pointer->areaCost << endl;
        }

        double eff = ((double) sum) / compressors[compPos].pointer->areaCost;


        if(eff <= 0.0001 && !generateSolution){   //not shure about accuracy - normally <= 0.0
            //conditionsMet = false;                //TODO: clean up
        }
/*
        if(usePreReduction && (!preReductionConditionsMet || eff < minimumEffReduction - 0.0001)){
            return -1.0;
        }
*/
        if(!conditionsMet){
            return -1.0;
        }
        else{
            return eff;
        }
    }

	//uses variableCompressor at stage s. Note that compType is the low compressor at variableBCompressors. (0, 3, 6,...)
	void BitHeapHeuristicCompression::useVariableCompressor(unsigned s, unsigned column, unsigned midLength, unsigned compType){


		//compType is first compressor in variableBCompressors
		if(newBits[s].size() > column){
			newBits[s][column] -= variableBCompressors[compType].height[0];
            if(!useHoles){
                if(newBits[s][column] < 0){
                    newBits[s][column] = 0;
                }
            }
		}
		//printBitHeap();
		for(unsigned m = 0; m < midLength; m++){
			if(newBits[s].size() > column + m + 1){
				newBits[s][column + m + 1] -= variableBCompressors[compType + 1].height[0];
				//cout << "subtracted " << variableBCompressors[compType + 1].height[0] << " in column " << column + m + 1 << endl;
				//printBitHeap();
				//cout << "__" << endl;
				if(!useHoles){
					if(newBits[s][column + m + 1] < 0){
						newBits[s][column + m + 1] = 0;
					}
				}
			}
		}

		unsigned bitPosition = column + midLength + 1;
		//heights of high compressor reversed
		for(int compPosition = variableBCompressors[compType + 2].height.size() - 1; compPosition >= 0; compPosition--){
			if(newBits[s].size() > bitPosition){
				newBits[s][bitPosition] -= variableBCompressors[compType + 2].height[compPosition];
				if(!useHoles){
					if(newBits[s][bitPosition] < 0){
						newBits[s][bitPosition] = 0;
					}
				}
			}
			bitPosition++;
		}


		//adding outputs:

		if(newBits[s + 1].size() > column){
			newBits[s + 1][column] += variableBCompressors[compType].outputs[0];
		}
		for(unsigned m = 0; m < midLength; m++){
			if(newBits[s + 1].size() > column + m + 1){
				newBits[s + 1][column + m + 1] += variableBCompressors[compType + 1].outputs[0];
			}
		}

		bitPosition = column + midLength + 1;
		//outputs of high compressor reversed
		for(int compPosition = variableBCompressors[compType + 2].outputs.size() - 1; compPosition >= 0; compPosition--){
			if(newBits[s + 1].size() > bitPosition){
				newBits[s + 1][bitPosition] += variableBCompressors[compType + 2].outputs[compPosition];
			}
			bitPosition++;
		}


		//add compressors to solution

		unsigned offset = compressors.size();
		solution[s].push_back(pair<int,int>(offset + compType + 0, column));	//low
		for(unsigned m = 0; m < midLength; m++){
			solution[s].push_back(pair<int,int>(offset + compType + 1, column + m + 1));	//mid
		}
		solution[s].push_back(pair<int,int>(offset + compType + 2, column + midLength + 1));	//high
	}

    void BitHeapHeuristicCompression::useCompressor(unsigned s, unsigned column, unsigned newCompPos){

        //make changes to newBits
        for(unsigned i = 0; i < compressors[newCompPos].pointer->getNumberOfColumns(); i++){
            newBits[s][column + i] -= compressors[newCompPos].pointer->getColumnSize(i);
            if(!useHoles){
                if(newBits[s][column + i] < 0){
                    newBits[s][column + i] = 0;
                }
            }

        }
        for(unsigned i = 0; i < (unsigned) compressors[newCompPos].pointer->getOutputSize(); i++){
            newBits[s + 1][column + i]++;
        }

        //add compressor to solution
        solution[s].push_back(pair<int,int>(compressors[newCompPos].originalPosition,column));
    }

    int BitHeapHeuristicCompression::printBitHeap(){
        for(unsigned i = 0; i < newBits.size(); i++){
            for(unsigned j = 0; j < newBits[i].size(); j++){
                cout << newBits[i][j] << " ";
            }
            cout << endl;
        }

        return 0;
    }

    //this function reduces the stages s of the solution. this is only possible if at least one presolution is passed
    //if useMaxHeuristicStageCount is false, the smallest s of all heuristic solutions is the s for the ILP-problem and
    //all presolutions with more stages are being deleted
    //otherwise the s in the ILP-problem is the s of the presolution, which needs the most stages
    void BitHeapHeuristicCompression::setUpNewStageCountOfILP(){
#ifdef HAVE_SCIP
        unsigned maxS = 0;

        if(!useMaxHeuristicStageCount){
            maxS = 1000;
        }

        unsigned sOfSolution[bitHeapILPCompression.heuristicSolutions.size()];

        for(unsigned i = 0; i < bitHeapILPCompression.heuristicSolutions.size(); i++){
            unsigned tempS = 0;

            for(unsigned s = 0; s < bitHeapILPCompression.heuristicSolutions[i].size(); s++){
                if(!bitHeapILPCompression.heuristicSolutions[i][s].empty()){
                    tempS = s;
                }
            }




            sOfSolution[i] = tempS;

            if(useMaxHeuristicStageCount){      //find maximal stage
                if(tempS > maxS){
                    maxS = tempS;
                }
            }
            else{                               //find minimal stage
                if(tempS < maxS){
                    maxS = tempS;
                }
            }
        }

        //now delete all solutions, which are bigger than maxS and resize all others
        // but only if !useMaxHeuristicStageCount
        if(!useMaxHeuristicStageCount){
            for(int i = (int) (bitHeapILPCompression.heuristicSolutions.size() - 1); i >= 0; i--){
                if(sOfSolution[i] != maxS){
                    bitHeapILPCompression.heuristicSolutions.erase(bitHeapILPCompression.heuristicSolutions.begin() + i);
                }
            }
        }

        //setting up noOfStages in ILP

        bitHeapILPCompression.getExternalStageCount = true;
        bitHeapILPCompression.noOfStages_ = (maxS + 2); //two - 1 because in the next stage(the last stage) is no compressor and one because we start at counting with 0

        cout << "new numberOfStages_ is " << (maxS + 2) << endl;
        if(useMaxHeuristicStageCount){
            cout << "numberOfStages_ was reduced so that all heuristic solutions are still valid" << endl;
        }
        else{
            cout << "numberOfstages_ was reduced as much as possible. heuristic solutions with more stages could have been deleted" << endl;
        }
#endif //HAVE_SCP
    }


    double BitHeapHeuristicCompression::computeAreaofSolution(){

        vector<BasicCompressor *> * comps = bh_->getPossibleCompressors();
        double areaSize = 0.0;
        double flipFlopCost = 1.0;
        if(bh_->getOp()->getTarget()->isPipelined())        {
            std::string targetID = bh_->getOp()->getTarget()->getID();
            if((targetID == "Virtex6") || (targetID == "Virtex7") || (targetID == "Spartan6")){
                flipFlopCost = 0.5; //there are two flip-flops per LUT for Virtex 6/7 and Spartan6
            }
            else{
                flipFlopCost = 1.0; //assume 1 for unknown device //!!!!
            }
        }
        else{
            flipFlopCost = 0.01; //nearly 0 for unpipelined designs
        }
		unsigned offset = compressors.size();
        for(unsigned i = 0; i < solution.size(); i++){
            list<pair<int, int> >:: iterator it;

            for(it = solution[i].begin(); it != solution[i].end(); it++){
                double areaOfCurrentCompressor = 0.0;
                if((*it).first == (int) compressors.size()){
                    areaOfCurrentCompressor = flipFlopCost; //cost of flipflop
                }
				else if((unsigned) (*it).first > compressors.size()){
					//variable basic Compressor
					if(((*it).first) - offset < variableBCompressors.size()){
						areaOfCurrentCompressor = variableBCompressors[(*it).first - offset].areaCost;
					}
					else{
						cout << "warning: non specified compressor number " << (*it).first << endl;
						cout << "there are " << variableBCompressors.size() << " variable basic compressors. ";
						cout << "but requesting variable basic compressor " << (*it).first - offset << endl;
					}
				}
                else{
                    areaOfCurrentCompressor = comps->at((*it).first)->areaCost;
                }
                //cout << "found compressor " << (*it).first << " with cost " << areaOfCurrentCompressor << " in stage " << i << endl;
                areaSize += areaOfCurrentCompressor;

            }
        }

		if(varCompSolution.size() > 0){

			for(unsigned s = 0; s < varCompSolution.size(); s++){
				list<variableCompressor>::iterator it;
				for(it = varCompSolution[s].begin(); it != varCompSolution[s].end(); it++){
					unsigned type = (*it).type * 3;
                    double areaOfCurrentCompressor = 0;
					areaOfCurrentCompressor += variableBCompressors[type].areaCost;
					areaOfCurrentCompressor += (*it).middleCompressorWidth * variableBCompressors[type + 1].areaCost;
					areaOfCurrentCompressor += variableBCompressors[type + 2].areaCost;
                    //cout << " found variable compressor " << (*it).type << " with cost " << areaOfCurrentCompressor << " in stage " << s << endl;
                    areaSize += areaOfCurrentCompressor;
				}

			}
		}

        return areaSize;
    }

    //now just functions for the ilp module

    int BitHeapHeuristicCompression::writeProblem(std::string filename){
#ifdef HAVE_SCIP
        bitHeapILPCompression.writeProblem(filename);
#endif //HAVE_SCIP
        return 0;
    }

    void BitHeapHeuristicCompression::setUpILPForMoreStages(unsigned stages, bool firstOne){
#ifdef HAVE_SCIP
        cout  << "in setupILPForMoreStages with stages = " << stages << endl;
        bitHeapILPCompression.solution.clear();
        bitHeapILPCompression.solution.resize(stages);
        if(!firstOne){
            bitHeapILPCompression.cleanUp();
        }
        unsigned zeros = 0;
        for(unsigned i = 0; i < sizeof(lowerBounds) / sizeof(lowerBounds[0]); i++){
           if(fabs(lowerBounds[i]) < 0.001){
               zeros++;
           }
        }

        //zeros does not result in a much faster ilp-solving process
        if(!useCompleteHeuristic ){
/*            if(!getExternalZero){
                bitHeapILPCompression.zeroStages = zeros;
                cout << "zeros added" << endl;
            }
            else{
                bitHeapILPCompression.zeroStages = exZero;
                cout << "external zeros added" << endl;
            }
*/
        }


        bitHeapILPCompression.getExternalStageCount = true;
        bitHeapILPCompression.useFixedStageCount = true;
        bitHeapILPCompression.useHeuristic = true;
        bitHeapILPCompression.dontAddFlipFlop = true; //!firstOne;  //first one: add flipflop otherwise not
        bitHeapILPCompression.noOfStages_ = stages;

        cout << "setting up done" << endl;
#endif //HAVE_SCIP
    }

    int BitHeapHeuristicCompression::solve(){
        cout << "in solve()" << endl;
        if((mode.compare("heuristic_parandeh-afshar_modified") == 0) || (mode.compare("heuristic_pa") == 0)) {
			if(useVariableCompressors){
				cout << "extracting the variable compressors" << endl;
				buildVariableCompressorSolution();
				cout << "extracting done" << endl;
			}
			double finalCost = computeAreaofSolution();
			cout << "the compressors have a cost of " << finalCost << endl;
            removeEmptyStages();
            return 0;
        }
#ifdef HAVE_SCIP
        if(differentStages){
            cout << "differentStages == true" << endl;
        }
        else{
            cout << "differentStages == false" << endl;
        }
        if(differentStages){
            printBitHeap();
            cout << "minFixedStage = " << minFixedStage << " and maxFixedStage = " << maxFixedStage << endl;
            bool takeFirstSolution = false; //we also check the next stages for a better solution

            for(unsigned s = minFixedStage; s <= maxFixedStage && !takeFirstSolution; s++){
                cout << "number of stages for ilp problem is " << s << endl;
                if(s != minFixedStage){
                    //setting up new bitHeapILPCompressioin
                    //bitHeapILPCompression = BitHeapILPCompression(bh_);
  //                  vector<list<pair<int,int> > > emptySolution;
    //                bitHeapILPCompression.solution = emptySolution;
                    cout << "calling cleanUp" << endl;
/*
                    bitHeapILPCompression.solution.clear();
                    bitHeapILPCompression.solution.resize(noOfStages_);
                    bitHeapILPCompression.cleanUp();
                    if(buildSingleStages){
                        bitHeapILPCompression.zeroStages = numberOfBuildStages;
                    }
                    else{
                        bitHeapILPCompression.zeroStages = 0;
                    }
                    bitHeapILPCompression.newBits = newBits;
                    bitHeapILPCompression.useHeuristic = true;
                    bitHeapILPCompression.dontAddFlipFlop = true;
                    bitHeapILPCompression.noOfStages_ = s;
                    bitHeapILPCompression.getExternalStageCount = true;
                    bitHeapILPCompression.useFixedStageCount = true;
                    bitHeapILPCompression.generateProblem();

*/
                    cout << "in loop" << endl;
                    setUpILPForMoreStages(s, false);
                    bitHeapILPCompression.generateProblem();
                }
                bitHeapILPCompression.writeProblem();
                bitHeapILPCompression.solve();
                cout << "solving done" << endl;
                if(!bitHeapILPCompression.infeasible){
                    cout << "solution is not infeasible" << endl;
                    if(bitHeapILPCompression.costOfCurrentSolution < minSliceCount){
                        minSliceCount = bitHeapILPCompression.costOfCurrentSolution;
                        solution.clear();
                        solution = bitHeapILPCompression.solution;
                        takeFirstSolution = true;      //uncomment, if you want to take the first solution, ilp finds (no restarts with a higher stagecount)
                    }
                }

                //bitHeapILPCompression.cleanUp();
            }

            if(bitHeapILPCompression.infeasible){       //if we have no solution after we went all the stages, exit
                exit(-1);
            }
            /*
            cout << endl << "before merging, solution" << endl;
            for(int j = 0; j < solution.size(); j++){
                list<pair<int,int> >:: iterator it;
                for(it = solution[j].begin(); it != solution[j].end(); it++){
                    cout << "applying compressor " << (*it).first << " to column " << (*it).second << " in stage " << j << endl;
                }

            }

            cout << endl << endl;

            cout << endl << "before merging, presolution" << endl;
            for(int j = 0; j < preSolution.size(); j++){
                list<pair<int,int> >:: iterator it;
                for(it = preSolution[j].begin(); it != preSolution[j].end(); it++){
                    cout << "applying compressor " << (*it).first << " to column " << (*it).second << " in stage " << j << endl;
                }

            }

            cout << endl << endl;
            */

            if(preSolution.size() > solution.size()){
                solution.resize(preSolution.size());
            }
            for(unsigned s = 0; s < preSolution.size(); s++){   //merge solution and presolution at the end

                solution[s].splice(solution[s].end(), preSolution[s]);
            }


            /*
            cout << "after SCIP" << endl;
            for(int j = 0; j < solution.size(); j++){
                list<pair<int,int> >:: iterator it;
                for(it = solution[j].begin(); it != solution[j].end(); it++){
                    cout << "applying compressor " << (*it).first << " to column " << (*it).second << " in stage " << j << endl;
                }

            }
            */
        }
        else{
            bitHeapILPCompression.solve();
            //copying solution from ILPCompression
            solution = bitHeapILPCompression.solution;
        }

        buildVariableCompressorSolution();

        //to fix the stagecount: delete zero stages

        unsigned tempMaxStage = solution.size();
        for(unsigned i = 0; i < solution.size(); i++){
            if(!solution[i].empty()){
                tempMaxStage = i;
            }
        }
        if(tempMaxStage + 1 < solution.size()){
            //first check if stages are really empty:
            cout << "we are deleting stages. initial stagesize is " << solution.size() << endl;
            unsigned leastEmptyStage = solution.size();
            for(unsigned i = solution.size() - 1; i >= tempMaxStage + 1; i--){
                if(solution[i].size() == 0){
                    leastEmptyStage = i;
                }
                else{
                    break;
                }
            }
            cout << "leastEmptyStage " << leastEmptyStage << endl;
            solution.resize(leastEmptyStage);          //cut empty stages off
            cout << "cut off of unused stages" << endl;
        }



        cout << "done" << endl;
#endif //HAVE_SCIP

		double finalCost = computeAreaofSolution();
        cout << "the compressors have a cost of " << finalCost << endl;

        removeEmptyStages();
        return 0;
    }

    //deletes empty stages if there are no later stages with compressors
    void BitHeapHeuristicCompression::removeEmptyStages(){
        bool foundEmpty = true;

        while(foundEmpty){
            foundEmpty = false;

            if(solution.size() >= 1 && solution[solution.size() - 1].size() == 0){
                solution.pop_back();
                foundEmpty = true;
            }
        }

        foundEmpty = true;
        while(foundEmpty){
            foundEmpty = false;

            if(varCompSolution.size() >= 1 && varCompSolution[varCompSolution.size() - 1].size() == 0){
                varCompSolution.pop_back();
                foundEmpty = true;
            }
        }

    }

    //this function goes through the solution, finds variableBasicCompressors, deletes them and
    //fills the variable compressor solution
    void BitHeapHeuristicCompression::buildVariableCompressorSolution(){



        int offset = compressors.size();

		for(unsigned variableOffset = 0; variableOffset < variableBCompressors.size(); variableOffset += 3){

			//right now only RCA
			int rcaLow = 0;			//change these three if necessary
			int rcaMid = 1;			//value is the position in ilpCompressions variableBCompressors
			int rcaHigh = 2;
			rcaLow += offset + variableOffset;
			rcaMid += offset + variableOffset;
			rcaHigh += offset + variableOffset;

			varCompSolution.resize(solution.size());

			for(unsigned s = 0; s < solution.size(); s++){

				bool found = true;
				while(found){
					found = false;
					variableCompressor tVarComp;
					tVarComp.type = (unsigned) (variableOffset / 3);	//three variablebasicCompressors (high mid low) build one variableCompressor
					tVarComp.startCompressorWidth = 1;
					tVarComp.middleCompressorWidth = 0;
					tVarComp.endCompressorWidth = 1;

					//now search for rca low basicVarCompressor
					bool lowFound = false;
					int currentColumn = 0;
					std::list<pair<int, int> >::iterator it;
					for(it = solution[s].begin(); it != solution[s].end(); it++){
						if((*it).first == rcaLow){
							//we found the start of a RCA
							tVarComp.column = (*it).second;
							currentColumn = (*it).second + 1;		//search for middle or high varBasicComp starting with next column
							solution[s].erase(it);
							lowFound = true;
							break;
						}
					}

					if(lowFound){
						bool highFound = false;
						bool foundSomething = true;
						while(!highFound && foundSomething){
							foundSomething = false;
							for(it = solution[s].begin(); it != solution[s].end(); it++){
								if((*it).second == currentColumn){
									if((*it).first == rcaMid){		//we found a suitable middle Compressor
										tVarComp.middleCompressorWidth++;
										currentColumn++;
										solution[s].erase(it);
										foundSomething = true;
										break;				//start again to search
									}
									else if((*it).first == rcaHigh){	//we found the high compressor
										highFound = true;
										solution[s].erase(it);
										foundSomething = true;
										break;
									}
								}
							}
						}

						//if we found a RCA compressor, add it to the varCompSolution
						if(highFound){
							varCompSolution[s].push_back(tVarComp);
							found = true;

						}
						else{
							cout << "=============" << endl;
							cout << "something went horrible wrong. a RCA compressor is not closed with a high compressor" << endl;
							cout << "=============" << endl;
						}
					}

				}

				//place here the search for other compressors with variable width

			}

		}


        //debug
        cout << "there were the following variable compressors:" << endl;
        for(unsigned s = 0; s < varCompSolution.size(); s++){
            for(list<variableCompressor>::iterator it = varCompSolution[s].begin(); it != varCompSolution[s].end(); it++){
                cout << "variable Compressor no. " << (*it).type << " starting at column " << (*it).column  << " and stage " << s  << " with middleSize of " << (*it).middleCompressorWidth << endl;
            }

        }


        //now print the other variables
        cout << "other compressors are:" << endl;
        for(unsigned s = 0; s < solution.size(); s++){
            std::list<pair<int, int> >::iterator it;
            for(it = solution[s].begin(); it != solution[s].end(); it++){
                unsigned type = (*it).first;
                unsigned column = (*it).second;
                cout << "normal Compressor " << type << " in column " << column  << " and stage " << s << endl;
            }
        }
    }

    void BitHeapHeuristicCompression::printLowerBounds(){
        cout << "PA-vector is (";
        for(unsigned i = 0; i < (sizeof(lowerBounds) / sizeof(lowerBounds[0])); i++){
            if(lowerBounds[i] > 4.1){
                cout << "inf";
                break;
            }
            else{
                cout << lowerBounds[i];
            }
            if(i + 1 != sizeof(lowerBounds) / sizeof(lowerBounds[0])){
                cout << "; ";
            }
        }
        cout << ")" << endl;
    }

    int BitHeapHeuristicCompression::plotSolution(){
#ifdef HAVE_SCIP
        bitHeapILPCompression.plotSolution();
#endif //HAVE_SCIP
        return 0;
    }

	void BitHeapHeuristicCompression::buildVariableCompressors(){

        if(useVariableCompressors)
        {
            std::string targetID = bh_->getOp()->getTarget()->getID();
            if((targetID == "Virtex5") || (targetID == "Virtex6") || (targetID == "Virtex7") || (targetID == "Spartan6"))
            {
                variableBasicCompressor c4_1;
                c4_1.areaCost = 1.0;
                c4_1.height = vector<int> (1);
                c4_1.outputs = vector<int> (1);
                c4_1.height[0] = 4;
                c4_1.outputs[0] = 1;
                variableBCompressors.push_back(c4_1);

                variableBasicCompressor c4_2;
                c4_2.areaCost = 1.0;
                c4_2.height = vector<int> (1);
                c4_2.outputs = vector<int> (1);
                c4_2.height[0] = 4;
                c4_2.outputs[0] = 2;
                variableBCompressors.push_back(c4_2);

                variableBasicCompressor c0_2;
                c0_2.areaCost = 1.0;
                c0_2.height = vector<int> (2);
                c0_2.outputs = vector<int> (2);
                c0_2.height[0] = 0;
                c0_2.height[1] = 2;         //inputs are reversed
                c0_2.outputs[0] = 1;
                c0_2.outputs[1] = 2;        //outputs are reversed
                variableBCompressors.push_back(c0_2);



                //second variable compressor. high is a (0;2) with cost of 0.5
                variableBasicCompressor c4_1b;
                c4_1b.areaCost = 1.0;
                c4_1b.height = vector<int> (1);
                c4_1b.outputs = vector<int> (1);
                c4_1b.height[0] = 4;
                c4_1b.outputs[0] = 1;
                variableBCompressors.push_back(c4_1b);

                variableBasicCompressor c4_2b;
                c4_2b.areaCost = 1.0;
                c4_2b.height = vector<int> (1);
                c4_2b.outputs = vector<int> (1);
                c4_2b.height[0] = 4;
                c4_2b.outputs[0] = 2;
                variableBCompressors.push_back(c4_2b);

                variableBasicCompressor c0_2b;
                c0_2b.areaCost = 0.5;
                c0_2b.height = vector<int> (1);
                c0_2b.outputs = vector<int> (1);
                c0_2b.height[0] = 0;
                c0_2b.outputs[0] = 2;
                variableBCompressors.push_back(c0_2b);

                cout << " there are " << variableBCompressors.size() << " variable basic compressors in variableBCompressors " << endl;
            }
        }
	}

    void BitHeapHeuristicCompression::setLowerBounds(string s){
        //parses the string efficiencyPerStage and writes the values to lowerBounds
        unsigned int sizeOfLowerBounds = (sizeof(lowerBounds) / sizeof(lowerBounds[0]));
        if(getLowerBoundsFromBitHeap){

            if(s.compare("") == 0){
                //efficiencyPerStage wasn't specified. therefore use complete heuristic ( = (0,0,0,...))
                for(unsigned j = 0; j < sizeOfLowerBounds; j++){
                    lowerBounds[j] = 0.0;
                }
                return;
            }


            unsigned int pos = 0;
            bool fraction = false;
            bool valueRead = false;
            double currentValue = 0.0;
            unsigned int fractionDigitCount = 0;
            bool foundInf = false;
            for(unsigned i = 0; i < s.size(); i++){
                if(valueRead && (s.at(i) == ' ' || s.at(i) == ',' || s.at(i) == ';')){    //new value
                    if(pos < sizeOfLowerBounds){
                        lowerBounds[pos] = currentValue;
                    }
                    else{
                        return;
                    }
                    pos++;
                    fraction = false;
                    currentValue = 0;
                    valueRead = false;
                    fractionDigitCount = 0;
                }
                else if(valueRead && s.at(i) == '.'){                   //fraction start
                    fraction = true;
                    valueRead = true;
                }
                else if(s.at(i) >= 48 && s.at(i) <= 57){                //digit
                    valueRead = true;
                    int tempValue = s.at(i) - 48;
                    if(!fraction){
                        currentValue *= 10;
                        currentValue += tempValue;
                    }
                    else{
                        double divider = 1.0;
                        fractionDigitCount++;
                        for(unsigned int j = 0; j < fractionDigitCount; j++){
                            divider *= 10;
                        }
                        double tempFraction = 1.0 / divider;
                        tempFraction *= tempValue;
                        currentValue += tempFraction;
                    }
                }
                                                                        //infinity
                else if((i + 2 < s.size()) && s.at(i) == 'i' && s.at(i + 1) == 'n' && s.at(i + 2) == 'f'){
                    if(pos < sizeOfLowerBounds){
                        lowerBounds[pos] = 10.0;
                    }
                    else{
                        return;
                    }
                    pos++;
                    foundInf = true;
                    break;      //we don't need to look further because there will be no bits in the next stage
                }

            }
            if(!foundInf && valueRead){
                //write last value
                if(pos < sizeOfLowerBounds){
                    lowerBounds[pos] = currentValue;
                }
                else{
                    return;
                }
                pos++;
            }




            //fill rest with infinity
            for(unsigned int j = pos; j < sizeOfLowerBounds; j++){
                lowerBounds[j] = 10.0;
            }
        }
    }

    int BitHeapHeuristicCompression::cleanUp(){
        if((mode.compare("heuristic_parandeh-afshar_modified") == 0) || (mode.compare("heuristic_pa") == 0)){
            return 0;
        }
#ifdef HAVE_SCIP
        return bitHeapILPCompression.cleanUp();
#else
        return 0;
#endif  //HAVE_SCIP
    }
}


//#endif
