#ifdef HAVE_SCIP

#include "BitHeapILPCompression.hpp"

namespace flopoco
{

BitHeapILPCompression::BitHeapILPCompression(BitHeap *bh)
{
        this->bh_ = bh;

        possibleCompressors_ = bh->getPossibleCompressors();
        sol = NULL;
        scip = NULL;
        noOfStagesUsed = -1;
    useVariableCompressors = true;

        srcFileName = bh->getOp()->getSrcFileName() + ":Bitheap:BitHeapILPCompression";

useFixedStageCount = false;

uniqueName_ = "BitHeapILPCompression for " + bh->getName();
}

BitHeapILPCompression::~BitHeapILPCompression()
{
        cleanUp();
}

int BitHeapILPCompression::getMaxStageCount(int maxHeight)
{
        int height=2;
        int stageCount=0;

        //the following loop computes the Dadda sequence [Dadda 1965] which is used as max. no of stages
        //(typically less are selected by the optimization using higher order compressors)
        while(height < maxHeight)
        {
                REPORT(DEBUG, "stageCount=" << stageCount << ", height=" << height);
                height = floor(height*3.0/2.0);
                stageCount++;
        }
        return stageCount;
}

int BitHeapILPCompression::generateProblem(){
/*
    //check for solution:
    if(solution.size() > 3){
        if(solution[2].size() > 0){
            cout << "solution is not cut off" << endl;
        }
    }
    cout << endl << endl << endl;
            */

    //set up first stage of U if there is no heuristic

    if(!useHeuristic){

        vector<int> firstCycle(bh_->bits.size());

        for(unsigned w=0; w < bh_->bits.size(); w++)
        {
            firstCycle[w] = bh_->bits[w].size();
        }

        newBits.push_back(firstCycle);

        useVariableCompressors = false;
        getExternalStageCount = false;
        zeroStages = 0;
    }

    //add flip flop to "possible compressors"
    vector<int> col(1);
    col[0]=1;
    flipflop = new BasicCompressor(bh_->getOp()->getTarget(),col);

    if(bh_->getOp()->getTarget()->isPipelined())
    {
        std::string targetID = bh_->getOp()->getTarget()->getID();
        if((targetID == "Virtex6") || (targetID == "Virtex7") || (targetID == "Spartan6"))
        {
            flipflop->areaCost = 0.5; //there are two flip-flops per LUT for Virtex 6/7 and Spartan6
        }
        else
        {
            flipflop->areaCost = 1; //assume 1 for unknown device //!!!!
        }
    }
    else
    {
        flipflop->areaCost = 0.01; //nearly 0 for unpipelined designs
    }
    REPORT(LIST, "Area cost for flip-flop is set to " << flipflop->areaCost);
    REPORT(DEBUG, "possibleCompressors before adding the flipflop: " << possibleCompressors_->size());
    if(!dontAddFlipFlop || !useHeuristic){
        possibleCompressors_->push_back(flipflop);
    }
    if(!getExternalStageCount){
        REPORT(DEBUG, "we set up the stages as bh_->getMaxHeight()");
        noOfStages_=getMaxStageCount(bh_->getMaxHeight());
    }
            unsigned noOfCompressors=possibleCompressors_->size();
    if(useVariableCompressors){
        noOfCompressors += variableBCompressors.size();
    }
    REPORT(DEBUG, "no of stages=" << noOfStages_);
    REPORT(DEBUG, "no of compressors=" << noOfCompressors);

    unsigned compOutputWordSizeMax = 0;
    for(unsigned e=0; e < possibleCompressors_->size(); e++)
    {
            if(compOutputWordSizeMax < (*possibleCompressors_)[e]->getOutputSize())
            {
                    compOutputWordSizeMax = (*possibleCompressors_)[e]->getOutputSize();
            }
    }
	compOutputWordSizeMax -= 1;
    REPORT(DEBUG,"compOutputWordSizeMax=" << compOutputWordSizeMax);

    noOfColumnsMax = newBits[0].size()+noOfStages_*(compOutputWordSizeMax-1);

    //if heuristic isn't used, fill up the U with zero - vectors
    if(!useHeuristic){
        for(int s = 1; s <= noOfStages_; s++){ //first one is filled with initial bits

            vector<int> zero;
            for(unsigned t = 0; t < bh_->bits.size(); t++){
                zero.push_back(0);
            }

            newBits.push_back(zero);
        }
    }

    if(scip != 0) cleanUp(); //clean up if used before

    /* initialize SCIP */
    SCIP_CALL( SCIPcreate(&scip) );
    SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
    SCIP_CALL( SCIPcreateProbBasic(scip, "example") );



    //create variables:
    compCountVars.clear();
    compCountVars.resize(noOfStages_+1);
    for(unsigned s=0; s < compCountVars.size(); s++)
    {
        compCountVars[s].resize(noOfCompressors);
    }

    SCIP_VAR* tmpvar;

    stageVars.resize(noOfStages_+1);
    for(unsigned s=0; s < stageVars.size(); s++)
    {
        stringstream varName;
        varName << "D_" << s;
        SCIP_CALL( SCIPcreateVarBasic(scip, &tmpvar, varName.str().c_str(), 0, 1, 0, SCIP_VARTYPE_BINARY ) );
        stageVars[s] = tmpvar;
        SCIP_CALL( SCIPaddVar(scip, tmpvar) );
    }


    for(unsigned s=0; s < compCountVars.size(); s++)
    {
        for(unsigned e=0; e < compCountVars[s].size(); e++)
        {
            for(unsigned c=0; c < newBits[0].size()+s*(compOutputWordSizeMax-1); c++)
            {
                //cout << "compCountVars[s][e][c] is being processed " << s << " " << e << " " << c << " areaCost: " << (*possibleCompressors_)[e]->areaCost << endl;
                stringstream varName;

                varName << "k_" << s << "_" << e << "_" << c;

                if(!useVariableCompressors){
                    SCIP_CALL( SCIPcreateVarBasic(scip, &tmpvar, varName.str().c_str(), 0.0, SCIPinfinity(scip), (*possibleCompressors_)[e]->areaCost, SCIP_VARTYPE_INTEGER) );
                }
                else{
                    if(e < possibleCompressors_->size()){
                        SCIP_CALL( SCIPcreateVarBasic(scip, &tmpvar, varName.str().c_str(), 0.0, SCIPinfinity(scip), (*possibleCompressors_)[e]->areaCost, SCIP_VARTYPE_INTEGER) );
                    }
                    else{           //we are now adding variable basic compressors
                        SCIP_CALL( SCIPcreateVarBasic(scip, &tmpvar, varName.str().c_str(), 0.0, SCIPinfinity(scip), variableBCompressors[e - possibleCompressors_->size()].areaCost, SCIP_VARTYPE_INTEGER) );
                    }
                }
                compCountVars[s][e].push_back(tmpvar);
                SCIP_CALL( SCIPaddVar(scip, tmpvar) );
                if(useVariableCompressors && e >= possibleCompressors_->size()){
                    unsigned offset = possibleCompressors_->size();
                    //assume that every complete variable compressor exists of three parts: low middle high. Every high-compressor is at the last position of those three.
                    //therefore high is at offset + 2, + 5, + 8 ... (=> % 3 == 2)
                    if(   ((e - offset) % 3 == 2)   && c == (newBits[0].size()+s*(compOutputWordSizeMax-1)) - 1){ //problem with the high-compressor

                        //add for the high - compressor a k-variable which exceeds the boundary.
                        stringstream varName;
                        varName << "k_" << s << "_" << e << "_" << c+1;
                        SCIP_CALL( SCIPcreateVarBasic(scip, &tmpvar, varName.str().c_str(), 0.0, SCIPinfinity(scip), variableBCompressors[e - possibleCompressors_->size()].areaCost, SCIP_VARTYPE_INTEGER) );
                        compCountVars[s][e].push_back(tmpvar);
                        SCIP_CALL( SCIPaddVar(scip, tmpvar) );
                    }
                }


            }
        }
    }
    //cout << "declaration of k done" << endl;


    //start at second stage. bits of first stage are in U
    columnBitCountVars.clear();
    columnBitCountVars.resize(noOfStages_);
    for(unsigned s=1; s < columnBitCountVars.size() + 1; s++){
        for(unsigned c=0; c < newBits[0].size()+s*(compOutputWordSizeMax-1); c++){
            stringstream varName;
            varName << "N_" << s << "_" << c;
            SCIP_CALL( SCIPcreateVarBasic(scip, &tmpvar, varName.str().c_str(), 0.0, SCIPinfinity(scip), 0, SCIP_VARTYPE_INTEGER) );
            columnBitCountVars[s - 1].push_back(tmpvar);
            SCIP_CALL( SCIPaddVar(scip, tmpvar) );
        }
    }
    //cout << "declaration of N done" << endl;
    //now create U and store them in newBitsCountVars
    newBitsCountVars.clear();
    newBitsCountVars.resize(noOfStages_+1);
    for(unsigned s=0; s < newBitsCountVars.size(); s++){
        vector<SCIP_VAR*> tempVec;
        for(unsigned c=0; c < (((unsigned)bh_->bits.size())+s*(compOutputWordSizeMax-1)); c++){
            stringstream varName;
            varName << "U_" << s << "_" << c;
            SCIP_CALL( SCIPcreateVarBasic(scip, &tmpvar, varName.str().c_str(), -SCIPinfinity(scip), SCIPinfinity(scip), 0, SCIP_VARTYPE_INTEGER) );
            tempVec.push_back(tmpvar);
            SCIP_CALL( SCIPaddVar(scip, tmpvar) );
        }
        newBitsCountVars[s] = tempVec;
    }

    //cout << "declaration of U done" << endl;
    //add constraints:
    SCIP_CONS* tmpcons;
    //add constraint C0: fill all the U's with values


    //cout << newBitsCountVars.size() << " " << newBits.size() << endl;
    for(unsigned i = 0; i < newBitsCountVars.size(); i++){
        vector<SCIP_VAR*> tempVectorCons = newBitsCountVars.at(i);
        vector<int> tempVectorValue;
        if(newBits.size() <= i){
            //zero vector
            tempVectorValue.resize(tempVectorCons.size());
        }
        else{
            tempVectorValue = newBits.at(i);
        }

        //cout << i << endl;
        //cout << tempVectorCons.size() << " " << tempVectorValue.size() << endl;
        for(unsigned j = 0; j < tempVectorCons.size(); j++){
            stringstream consName;
            consName << "C0_" << i << "_" << j;

            //set the value to zero if it exceeds the newBits
            int value;
            if(j >= tempVectorValue.size()){
                value = 0;
            }
            else{
                value = tempVectorValue.at(j);
            }
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &tmpcons, consName.str().c_str(), 0, NULL, NULL, value, value) );
            SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons,tempVectorCons.at(j), 1.0) );
            SCIP_CALL( SCIPaddCons(scip, tmpcons) );
        }

    }
    //cout << compCountVars.size() << " is compCountVars.size()" << endl;
    if(compCountVars[0][0][0] != NULL){
        //cout << "0-0-0 is full" << endl;
    }
    else{
        //cout << "0-0-0 is not full" << endl;
    }
    //cout << "declaration of C0 done" << endl;
    const int LARGE_NUMBER = 10000;
//  int LARGE_NUMBER = bh_->getMaxHeight()+1; //!!!









    //add constraints C1 and C2:


    for(unsigned s=0; s < compCountVars.size()-1; s++){



        //add constraint C1:
        for(int c=0; c < ((int) (newBits[0].size()+s*(compOutputWordSizeMax-1))); c++){
            stringstream consName;
            consName << "C1_" << s << "_" << c;
            //cout << "s = " << s << " and c = " << c << endl;
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &tmpcons, consName.str().c_str(), 0, NULL, NULL, 0, SCIPinfinity(scip)) );
            for(unsigned e=0; e < possibleCompressors_->size(); e++){
                for(int ce=0; ce < (int) (*possibleCompressors_)[e]->height.size(); ce++){
                    if(c-ce >= 0){
                        //cout << "call s#e#c-ce: " << s << e << c-ce << endl;
                        SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons, compCountVars[s][e][c-ce] , (*possibleCompressors_)[e]->height[(*possibleCompressors_)[e]->height.size()-ce-1]) );
                        //cout << "after call" << endl;
                    }
                }
            }
            if(useVariableCompressors){
                for(unsigned e=possibleCompressors_->size(); e < possibleCompressors_->size() + variableBCompressors.size(); e++){
                    for(int ce=0; ce < (int) variableBCompressors[e - possibleCompressors_->size()].height.size(); ce++){
                        if(c-ce >= 0){
                            //cout << "call s#e#c-ce: " << s << e << c-ce << endl;
                            SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons, compCountVars[s][e][c-ce] , variableBCompressors[e - possibleCompressors_->size()].height[variableBCompressors[e - possibleCompressors_->size()].height.size()-ce-1]) );
                            //cout << "after call" << endl;
                        }
                    }
                }
            }
            //cout << "after loop" << endl;
            if(s != 0){
                //first stage is in U therefore N starts with s-1
                SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons, columnBitCountVars[s - 1][c] , -1) );
            }
            //subtracting U as well
            SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons, newBitsCountVars[s][c] , -1) );
            SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons, stageVars[s] , LARGE_NUMBER) );
            SCIP_CALL( SCIPaddCons(scip, tmpcons) );
        }

        //cout << "filling of C1 at stage " << s << " done" << endl;

        //-----------------------------------------------------------------------

        //add constraint C2:
        for(int c=0; c < ((int) (newBits[0].size()+(s+1)*(compOutputWordSizeMax-1))); c++){
            stringstream consName;
            consName << "C2_" << s+1 << "_" << c;

            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &tmpcons, consName.str().c_str(), 0, NULL, NULL, 0, 0) );

            for(unsigned e=0; e < possibleCompressors_->size(); e++){
                for(int ce=(int) (*possibleCompressors_)[e]->outputs.size()-1; ce >= 0; ce--){
                    if((c-ce >= 0) && (c-ce < ((int) (newBits[0].size()+s*(compOutputWordSizeMax-1))))){
                        SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons, compCountVars[s][e][c-ce] , (*possibleCompressors_)[e]->outputs[(*possibleCompressors_)[e]->outputs.size()-ce-1]) );
                    }
                }
            }
            if(useVariableCompressors){
                unsigned offset = possibleCompressors_->size();
                for(unsigned e=offset; e < offset + variableBCompressors.size(); e++){
                    for(int ce=(int) variableBCompressors[e - offset].outputs.size()-1; ce >= 0; ce--){
                        if((c-ce >= 0) && (c-ce < ((int) (newBits[0].size()+s*(compOutputWordSizeMax-1))))){
                            SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons, compCountVars[s][e][c-ce] , variableBCompressors[e - offset].outputs[variableBCompressors[e - offset].outputs.size()-ce-1]) );
                            /*
                            if(e == 14){
                                cout << "added " << e << " " << c-ce << " normal" << endl;
                            }
                            */
                        }
                    }
                    
                    
                    //if(e == noOfCompressors - 1 && variableBCompressors[e - offset].height[0] == 0 && c == ((int) (newBits[0].size()+s*(compOutputWordSizeMax-1))) ){
                    
                    //assume that every complete variable compressor exists of three parts: low middle high. Every high-compressor is at the last position of those three.
                    //therefore high is at offset + 2, + 5, + 8 ... (=> % 3 == 2)
                    if(   ((e - offset) % 3 == 2)     && c == ((int) (newBits[0].size()+s*(compOutputWordSizeMax-1))) ){
                        for(unsigned tempOutput = 0; tempOutput < variableBCompressors[e - offset].outputs.size(); tempOutput++){
                            //cout << "s = " << s << " and tempOutput = " << tempOutput << endl;
                            //cout << "c = " << c << " and compCountVars[s][e].size() = " << compCountVars[s][e].size() << endl;
                            SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons, compCountVars[s][e][c] , variableBCompressors[e - offset].outputs[tempOutput]) );
                        }
                        
                    }
                }
            }
            SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons, columnBitCountVars[s][c] , -1) );
            SCIP_CALL( SCIPaddCons(scip, tmpcons) );


        }
        //cout << "filling of C2 at stage " << s << " done" << endl;
    }
    
    //cout << "C1 and C2 totally done" << endl;


    //------------------------------------------------------------

    //add constraint C3:
    for(unsigned s=0; s < compCountVars.size(); s++){
        for(unsigned c=0; c < newBitsCountVars[s].size(); c++){
            stringstream consName;
            consName << "C3_" << s;
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &tmpcons, consName.str().c_str(), 0, NULL, NULL, -SCIPinfinity(scip), LARGE_NUMBER+2) );

            if(s != 0){//in first stage there is no N
                SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons, columnBitCountVars[s - 1][c] , 1) );
            }
            SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons, newBitsCountVars[s][c] , 1) );

            //now make sure that the U's in later stages are empty
            for(unsigned z = s + 1 ; z < compCountVars.size(); z++){
                SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons, newBitsCountVars[z][c] , 4) );
            }
            SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons, stageVars[s] , LARGE_NUMBER) );
            SCIP_CALL( SCIPaddCons(scip, tmpcons) );
        }

    }

    //cout << "C3 filling done" << endl;

    //-----------------------------------------------------------

    //add constraint C4:
    SCIP_CALL( SCIPcreateConsBasicLinear(scip, &tmpcons, "C4", 0, NULL, NULL, 1, 1) );
    for(unsigned s=0; s < stageVars.size(); s++){
        SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons, stageVars[s] , 1) );
    }
    SCIP_CALL( SCIPaddCons(scip, tmpcons) );


    //-----------------------------------------------------------

    //add constraint C4_1: setting up the 2-bit-adder stage

    if(useFixedStageCount){
        SCIP_CALL( SCIPcreateConsBasicLinear(scip, &tmpcons, "C4_1", 0, NULL, NULL, 1, 1) );
        SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons, stageVars[noOfStages_], 1) );
        SCIP_CALL( SCIPaddCons(scip, tmpcons) );
    }

    //cout << "C4 filling done" << endl;


    //----------------------------------------------------------

    if(useVariableCompressors){

        cout << "in useVariableCompressors " << endl << endl << endl;
        //assume that the order in variableBCompressors is low - middle -high ( (3;1) - (2;1) - (0;1) )

        unsigned offset = possibleCompressors_->size();

        for(unsigned s = 0; s < compCountVars.size() - 1; s++){

            for(int c=0; c < ((int) (newBits[0].size()+s*(compOutputWordSizeMax-1))); c++){
                //C5:

                //format : k_high_c=i+1 + k_mid_c=i+1 = k_mid_c=i + k_low_c=i
                //setting k_high_0 and k_mid_0 to zero is done after this loop

                stringstream consName;
                consName << "C5_" << s << "_" << c;
                SCIP_CALL( SCIPcreateConsBasicLinear(scip, &tmpcons, consName.str().c_str(), 0, NULL, NULL, 0, 0) );

                SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons,compCountVars[s][offset + 2][c + 1], 1.0) );
                if(c != ((int) (newBits[0].size()+s*(compOutputWordSizeMax-1))) - 1){   //if it is not the highest column, in the next higher column can be a middle compressor
                    SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons,compCountVars[s][offset + 1][c + 1], 1.0) );
                }
                SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons,compCountVars[s][offset + 0][c], -1.0) );
                SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons,compCountVars[s][offset + 1][c], -1.0) );            //new: subtract the middle compressors as well

                SCIP_CALL( SCIPaddCons(scip, tmpcons) );

                /*  not longer needed. C5 and the restrictions regarding high and mid part in lowest column is enough
                //C6:
                stringstream consName2;
                consName2 << "C6_" << s << "_" << c;
                SCIP_CALL( SCIPcreateConsBasicLinear(scip, &tmpcons, consName2.str().c_str(), 0, NULL, NULL, 0, 0) );

                SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons,compCountVars[s][offset + 2][c + 1], 1.0) );
                if(c != ((int) (newBits[0].size()+s*(compOutputWordSizeMax-1))) - 1){   //if it is not the highest column, in the next higher column can be a middle compressor
                    SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons,compCountVars[s][offset + 1][c + 1], 1.0) );
                }
                SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons,compCountVars[s][offset + 1][c], -1.0) );    //(offset + 1 ) only difference to C5
                SCIP_CALL( SCIPaddCons(scip, tmpcons) );

                */
            }

            //make sure that the high compressor and middle compressor are not used in the lowest column.

            stringstream consNameB;
            consNameB << "C5b_" << s;
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &tmpcons, consNameB.str().c_str(), 0, NULL, NULL, 0, 0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons,compCountVars[s][offset + 2][0], 1.0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons,compCountVars[s][offset + 1][0], 1.0) );
            SCIP_CALL( SCIPaddCons(scip, tmpcons) );

        }
    }




    //add constraint C7: set maximum of slicecount. the maximum is received by setting lower D_s to 1.

/*
    double tempMax = 28;

    SCIP_CALL( SCIPcreateConsBasicLinear(scip, &tmpcons, "C7", 0, NULL, NULL, 0, tempMax - 0.5) );

    for(unsigned s = 0; s < lastStage; s++){
        for(unsigned e = 0; e < compCountVars[s].size(); e++){
            for(unsigned c = 0; c < compCountVars[s][e].size(); c++){
                SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons, compCountVars[s][e][c], (*possibleCompressors_)[e]->areaCost));
            }
        }
    }
    SCIP_CALL( SCIPaddCons(scip, tmpcons));
*/

    // now set k = 0 if we solved some of the stages with a heuristic
    //C8

    //cout << "zerostages = " << zeroStages << endl;
    for(unsigned i = 0; i < zeroStages; i++){
        vector<vector<SCIP_VAR*> > tempCompCountVars = compCountVars.at(i);

        for(unsigned j = 0; j < tempCompCountVars.size(); j++){

            vector<SCIP_VAR*> tempCompCountVars2 = tempCompCountVars.at(j);

            for(unsigned k = 0; k < tempCompCountVars2.size(); k++){
                stringstream consName;
                consName << "C8_" << i << "_" << j << "_" << k;
                SCIP_CALL( SCIPcreateConsBasicLinear(scip, &tmpcons, consName.str().c_str(), 0, NULL, NULL, 0, 0) );
                SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons,tempCompCountVars2.at(k), 1.0) );
                SCIP_CALL( SCIPaddCons(scip, tmpcons) );
            }


        }
    }
    //cout << "added constraint C8" << endl;


    //C9: settings N's to zero
    for(unsigned i = 0; i < zeroStages; i++){       //i = 0 are the N's from stage 1
        vector<SCIP_VAR*> tempColumnBitCountVars = columnBitCountVars.at(i);

        for(unsigned j = 0; j < tempColumnBitCountVars.size(); j++){
            stringstream consName;
            consName << "C9_" << i << "_" << j;
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &tmpcons, consName.str().c_str(), 0, NULL, NULL, 0, 0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons,tempColumnBitCountVars.at(j), 1.0) );
            SCIP_CALL( SCIPaddCons(scip, tmpcons) );
        }
    }

    //cout << "end of SCIP problem description" << endl;
    return 0;
}





int BitHeapILPCompression::writeProblem(std::string filename){
    const char *fname=NULL;;
    if(!filename.empty()){
        fname = filename.c_str();
    }
    //print model in LP format:
    SCIP_RESULT result;
    SCIP_CALL(SCIPwriteLp(scip, NULL, fname, FALSE, SCIP_OBJSENSE_MINIMIZE, 1.0, 0.0, SCIPgetVars(scip), SCIPgetNVars(scip), SCIPgetNBinVars(scip), SCIPgetNIntVars(scip), SCIPgetNImplVars(scip), SCIPgetNContVars(scip), SCIPgetConss(scip), SCIPgetNOrigConss(scip), &result));
    return 0;
}

bool BitHeapILPCompression::solve(){
    //set a time limit:
    float timeout;
    timeout = bh_->getOp()->getTarget()->ilpTimeout();

    SCIP_CALL( SCIPsetRealParam(scip,"limits/time",timeout) );

    REPORT(INFO, "ILP solver timeout is set to " << timeout);

//  SCIP_CALL ( SCIPsetRealParam(scip, "limits/memory", 10240) );

//  REPORT(INFO, "ILP solver memory limit is set to 10240 MB");

    //prioritize heurisitics:
/*
    SCIP_CALL ( SCIPsetIntParam(scip, "heuristics/crossover/priority", 50000000) );
    SCIP_CALL ( SCIPsetIntParam(scip, "heuristics/pscostdiving/priority", 40000000) );
    SCIP_CALL ( SCIPsetIntParam(scip, "heuristics/linesearchdiving/priority", 30000000) );
    SCIP_CALL ( SCIPsetIntParam(scip, "heuristics/veclendiving/priority", 10000000) );
*/

    //solve the problem:
    SCIP_CALL( SCIPsolve(scip) );

    SCIP_STATUS status;
    status = SCIPgetStatus(scip);
    if(status == SCIP_STATUS_INFORUNBD && !useFixedStageCount){
        THROWERROR("No optimal solution found (problem infeasible or unbounded!)");
    }
    else if(status == SCIP_STATUS_INFORUNBD && useFixedStageCount){
        REPORT(INFO, "No optimal solution for current stage count found (problem infeasible or unbound!).\n We will try it again with more stages");
        infeasible = true;
    }

    sol = SCIPgetBestSol(scip);

    if(!sol && !useFixedStageCount){     //we need to run ilp-solving several times when trying to determine the minimal amount of stages
        THROWERROR("No feasible solution found within the ILP timeout of " << timeout << " seconds");
    }
    else{
        SCIP_Real area = SCIPgetPrimalbound(scip);
        preReductionAreaCost += (double)area;
        REPORT(DEBUG, "The size of the complete Solution is " << preReductionAreaCost << " LUTs");
        costOfCurrentSolution = (double)area;
        infeasible = false;
    }

    /*
    SCIP_Bool feasible;
    SCIP_CALL( SCIPcheckSol ( scip,sol,true,true,true,true,&feasible) );

    if(!feasible)
    {
            cerr << "No feasible solution found within the ILP timeout of " << timeout << " seconds" << endl;
            exit(-1);
    }
    */

    SCIP_VAR** allVariables = SCIPgetVars(scip);
    int noOfallVariables = SCIPgetNVars(scip);
    SCIP_Real allValues[noOfallVariables];

    //display solution:
    SCIP_CALL( SCIPgetSolVals   (scip,sol,noOfallVariables,allVariables,allValues));

    REPORT(DEBUG, "solution:");
    for(int i=0; i < noOfallVariables; i++){
        if(fabs((allValues[i])) > 1E-4){
            REPORT(DEBUG, SCIPvarGetName(allVariables[i]) << "=" << allValues[i]);
        }
    }

    noOfStagesUsed=-1; //the actual number of used compressor stages
    for(unsigned s=0; s < stageVars.size(); s++){
        if(round(SCIPgetSolVal(scip,sol,stageVars[s])) == 1){
            noOfStagesUsed = s;
        }
    }

    if(!infeasible){
        if(!useHeuristic){
            solution.resize(compCountVars.size());
        }
        //convert solution to a more portable data structure:
        REPORT(DEBUG, "adding compressors to solution with noOfStages = " << noOfStages_);
        SCIP_Real val;
        for(int s=0; s <= noOfStagesUsed; s++)
        {
            for(unsigned e=0; e < compCountVars[s].size(); e++)
            {
                for(int c=0; c < noOfColumnsMax; c++)
                {
                    if(c < ((int) compCountVars[s][e].size()))
                    {
                        val = SCIPgetSolVal(scip,sol,compCountVars[s][e][c]);
                        int k_max = ((int) round(val));
                        if(k_max > 0)
                        {
                            for(int k=0; k < k_max; k++)
                            {
                                solution[s].push_back(pair<int,int>(e,c));
                            }
                        }
                    }
                }
            }
        }
    }
    REPORT(DEBUG, "SCIP done in BitHeapILPCompression");
    //cout << solution.size() << endl;
    for(unsigned j = 0; j < solution.size(); j++){
        list<pair<int,int> >:: iterator it;
        for(it = solution[j].begin(); it != solution[j].end(); it++){
            REPORT(DEBUG, "adding compressor to solution: " << (*it).first << " to column " << (*it).second << " in stage " << j);
        }

    }
    //quick fix: delete empty stages
    unsigned int realStagesUsed = 0;
    for(unsigned j = 0; j < solution.size(); j++){
        if(solution[j].size() > 0){
            realStagesUsed = j;
        }
    }
    //cout << "before resizing solution: in ilpCompression " << solution.size() << endl;
    solution.resize(realStagesUsed + 1);
    //cout << "after resizing solution: in ilpCompression " << solution.size() << endl;

    return infeasible;
}


//this function checks the vector heuristicSolutions for solutions and passes them to
//the scip solver
int BitHeapILPCompression::passHeuristicSolutions(){
    int successFullPasses = 0;

    SCIP_SOL* heuSol;

    //cout << "in passHeuristicSolution" << endl;

    for(unsigned i = 0; i < heuristicSolutions.size(); i++){
        REPORT(DEBUG, "passing heuristic solution number " << i);
        SCIP_CALL( SCIPcreateSol(scip, &heuSol, NULL) );
        computeCompressorCount(i);      //computes also heuristicN

        //U: is the same for every solution
        for(unsigned s = 0; s < newBits.size(); s++){
            for(unsigned c = 0; c < newBits[s].size(); c++){

                //check whether we have the SCIP_VAR in newBitsCountVars
                if(newBitsCountVars.size() > s){
                    if(newBitsCountVars[s].size() > c){
                        //set the Var in heuSol
                        SCIP_CALL( SCIPsetSolVal(scip, heuSol, newBitsCountVars[s][c], newBits[s][c]));
                    }
                }

            }
        }
        //cout << "filling U values done" << endl;
        //N: only has s - 1 stages; s = 0 is the second stage
        for(unsigned s = 0; s < heuristicN.size(); s++){
            for(unsigned c = 0; c < heuristicN[s].size(); c++){
                if(columnBitCountVars.size() > s){
                    if(columnBitCountVars[s].size() > c){
                        SCIP_CALL( SCIPsetSolVal(scip, heuSol, columnBitCountVars[s][c], heuristicN[s][c]));
                    }
                }
            }
        }
        //cout << "filling N values done" << endl;
        //k:
        for(unsigned s = 0; s < compressorCount.size(); s++){
            for(unsigned e = 0; e < compressorCount[s].size(); e++){
                for(unsigned c = 0; c < compressorCount[s][e].size(); c++){
                    if(compCountVars.size() > s){
                        if(compCountVars[s].size() > e){
                            if(compCountVars[s][e].size() > c){
                                SCIP_CALL( SCIPsetSolVal(scip, heuSol, compCountVars[s][e][c], compressorCount[s][e][c]));
                            }
                        }
                    }
                }
            }
        }
        //cout << "filling k values done" << endl;
        //D: find D by going through the solution. D_i is true, if stage i-1 is the last stage with compressors in the current heuristicSolution
        unsigned outputStage = heuristicSolutions[i].size() - 1;
        for(int s = (int) heuristicSolutions[i].size() - 1; s >= 0; s--){
            if(heuristicSolutions[i][s].size() > 0){
                break;
            }
            else{
                outputStage = s;
            }
        }

        if(stageVars.size() > outputStage){
            SCIP_CALL( SCIPsetSolVal(scip, heuSol, stageVars[outputStage], 1));
        }
        cout << "filling D values done" << endl;

        SCIP_Bool stored;

        SCIP_CALL (SCIPaddSolFree(scip, &heuSol, &stored));

        if(stored == 1){
            successFullPasses++;
            cout << "pass was successfull" << endl;
        }
        else{
            cout << "pass was not successfull" << endl;
        }



    }



    return successFullPasses;
}


//this counts the compressors used in the following format:
//s - e - c : s = stage, e = compressing element, c = column
void BitHeapILPCompression::computeCompressorCount(int pos){
    compressorCount.clear();

    //generate size
    compressorCount.resize(heuristicSolutions[pos].size());
    for(unsigned s = 0; s < heuristicSolutions[pos].size(); s++){

        compressorCount[s].resize(possibleCompressors_->size());
        for(unsigned e = 0; e < possibleCompressors_->size(); e++){

            compressorCount[s][e].resize(newBits[s].size());
        }

    }
    //now we have the right size and start to fill it up with values
    //cout << "sizing in computeCompressorCount done" << endl;
    for(unsigned s = 0; s < heuristicSolutions[pos].size(); s++){
        list<pair<int,int> >:: iterator it;
        for(it = heuristicSolutions[pos][s].begin(); it != heuristicSolutions[pos][s].end(); it++){
            int tempE = (*it).first;
            int tempC = (*it).second;
            compressorCount[s][tempE][tempC]++;
        }
    }
    //cout << "filling in computeCompressorCount done" << endl;
    computeHeuristicN();
}

//computes the N of the heuristic solution for passing it to the scip-solver
//you have to fill compressorCount first
//NOTE: works only with compressors, whose output is one in every column
void BitHeapILPCompression::computeHeuristicN(){
    heuristicN.clear();
    //printNewBits();

    //resizing:
    heuristicN.resize(compressorCount.size());
    for(unsigned s = 0; s < heuristicN.size(); s++){
        heuristicN[s].resize(newBits[s].size());
    }


    for(unsigned e = 0; e < possibleCompressors_->size(); e++){
        unsigned outputSize = (*possibleCompressors_)[e]->getOutputSize();
        for(unsigned s = 0; s < heuristicN.size(); s++){
            for(unsigned c = 0; c < heuristicN[s].size(); c++){
                if(compressorCount[s][e][c] > 0){
                    for(unsigned i = 0; i < outputSize; i++){
                        //cout << heuristicN[s + 1].size() << endl;
                        heuristicN[s + 1][c + i] += compressorCount[s][e][c];
                    }
                }
            }
        }
    }

    //now delete the first stage because there are only zeros
    heuristicN.erase(heuristicN.begin());
}




    void BitHeapILPCompression::plotSolution()
    {
        if(sol == NULL) return;

        SCIP_Real obj;
        obj = SCIPgetPrimalbound(scip);

        //display bit heap stages:
        cerr << "bit col: ";
        for(int c=noOfColumnsMax-1; c >= 0; c--)
        {
          cerr << c << '\t';
        }
        cerr << endl;

        SCIP_Real val;
        int overheadTotal=0;
        for(int s=0; s <= noOfStagesUsed; s++)
        {
          cerr << "-----------------------------------------------------------------------------------------------------------------" << endl;
          cerr << "stage " << s << ": ";
          int colmnHeight[noOfColumnsMax];

          for(int c=noOfColumnsMax-1; c >= 0; c--)
          {
            if(c < ((int) columnBitCountVars[s].size()))
            {
              val = round(SCIPgetSolVal(scip,sol,columnBitCountVars[s][c]));
            }
            else
            {
              val = 0;
            }
            colmnHeight[c] = (int) val;
            cerr << val << "\t";
          }
          cerr << endl;

          for(unsigned e=0; e < compCountVars[s].size(); e++)
          {
            for(int c=noOfColumnsMax-1; c >= 0; c--)
            {
              if(c < ((int) compCountVars[s][e].size()))
              {
                val = SCIPgetSolVal(scip,sol,compCountVars[s][e][c]);
                int k_max = ((int) round(val));
                if(k_max > 0)
                {
                  for(int k=0; k < k_max; k++)
                  {
                    cerr << "       - ";
                    for(int i=0; i < noOfColumnsMax-((int) (*possibleCompressors_)[e]->height.size())-c; i++)
                    {
                      cerr << " \t";
                    }
                    for(int ce=0; ce < (int) (*possibleCompressors_)[e]->height.size(); ce++)
                    {
                      colmnHeight[c+(*possibleCompressors_)[e]->height.size()-ce-1] -= (*possibleCompressors_)[e]->height[ce];
                      cerr << (*possibleCompressors_)[e]->height[ce] << '\t';
                    }
                    cerr << endl;
                  }
                }
              }
            }
          }

          cerr << "       = ";
          for(int c=noOfColumnsMax-1; c >= 0; c--)
          {
            cerr << colmnHeight[c] << '\t';
          }
          cerr << endl;

          int overheadStage=0;
          for(int c=0; c < noOfColumnsMax; c++)
          {
            if(s == noOfStagesUsed)
            {
              if(colmnHeight[c] > 0)
                overheadStage += 2-colmnHeight[c];
            }
            else
            {
              overheadStage += -colmnHeight[c];
            }
          }
          cerr << "overhead (stage " << s << ")=" << overheadStage << endl;
          overheadTotal += overheadStage;
        }


        cerr << "overhead (total)=" << overheadTotal << endl;
        cerr << "no of compressor stages: " << noOfStagesUsed << endl;
        cerr << "objective value: " << obj << endl;
    }

    int BitHeapILPCompression::cleanUp()
    {
        REPORT(DEBUG, "cleaning up SCIP... ");
        //clean up:


        if(compressionType == 6 || compressionType == 7){
            //heuristic only. therefore no variables for ilp has been produces
            //and we don't need to clean up
            BMScheckEmptyMemory();

            scip = NULL;
            sol = NULL;
            return 0;
        }

        if(!compCountVars.empty()){
            //cout << "not empty" << endl;
            for(unsigned s=0; s < compCountVars.size(); s++)
            {
              for(unsigned e=0; e < compCountVars[s].size(); e++)
              {
                for(vector<SCIP_VAR*>::iterator it = compCountVars[s][e].begin(); it != compCountVars[s][e].end(); ++it)
                {
                    //cout << "SCIPfreeBufferArray" << endl;
                    //SCIPfreeBufferArray(scip, &(*it));
                    //cout << "SCIPreleaseVar" << endl;
                    SCIP_CALL( SCIPreleaseVar(scip, &(*it)) );

                }
              }
            }
        }


        for(unsigned s = 0; s < columnBitCountVars.size(); s++){
            for(vector<SCIP_VAR*>::iterator it = columnBitCountVars[s].begin(); it != columnBitCountVars[s].end(); it++){
                SCIP_CALL( SCIPreleaseVar(scip, &(*it)) );
            }
        }

        for(vector<SCIP_VAR*>::iterator it = stageVars.begin(); it != stageVars.end(); it++){
            SCIP_CALL( SCIPreleaseVar(scip, &(*it)) );
        }

        for(unsigned s = 0; s < newBitsCountVars.size(); s++){
            for(vector<SCIP_VAR*>::iterator it = newBitsCountVars[s].begin(); it != newBitsCountVars[s].end(); it++){
                SCIP_CALL( SCIPreleaseVar(scip, &(*it)) );
            }
        }

        SCIP_CALL( SCIPfree(&scip) );

        BMScheckEmptyMemory();

        scip = NULL;
        sol = NULL;

        //cout << "cleaning up done " << endl;
        return 0;
    }

void BitHeapILPCompression::printNewBits(){
    for(unsigned i = 0; i < newBits.size(); i++){
        for(unsigned j = 0; j < newBits[i].size(); j++){
            cerr << newBits[i][j] << " ";
        }
        cerr << endl;
    }
}






}   //end namespace flopoco

#endif // HAVE_SCIP
