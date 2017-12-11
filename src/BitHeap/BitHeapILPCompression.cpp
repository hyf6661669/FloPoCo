#ifdef HAVE_SCIP

#include "BitHeapILPCompression.hpp"
/*
#ifdef HAVE_SCALP
#include <ScaLP/Solver.h>
#include <ScaLP/Exception.h>    // ScaLP::Exception
#include <ScaLP/SolverDynamic.h> // ScaLP::newSolverDynamic
#endif //HAVE_SCALP
*/
namespace flopoco
{

BitHeapILPCompression::BitHeapILPCompression(BitHeap *bh)
{
        this->bh_ = bh;

        possibleCompressors_ = bh->getPossibleCompressors();
#ifdef HAVE_SCALP
        sol = ScaLP::Result();
#else
        sol = NULL;
#endif //HAVE_SCALP
        scip = NULL;
        noOfStagesUsed = -1;
        useVariableCompressors = bh->useVariableColumnCompressors;

        srcFileName = bh->getOp()->getSrcFileName() + ":Bitheap:BitHeapILPCompression";

        useFixedStageCount = false;

        uniqueName_ = "BitHeapILPCompression for " + bh->getName();

#ifdef HAVE_SCALP
        string usedILPSolver = UserInterface::ilpSolver;
        if(usedILPSolver.compare("Gurobi") == 0 && usedILPSolver.compare("CPLEX") == 0 && usedILPSolver.compare("SCIP") == 0 && usedILPSolver.compare("LPSolve") == 0){
            THROWERROR("ilpSolver " << usedILPSolver << " unknown!");
        }
        ScaLP::Solver s = ScaLP::Solver(ScaLP::newSolverDynamic({usedILPSolver}));
        //ScaLP::Solver s = ScaLP::Solver(ScaLP::newSolverDynamic({"Gurobi","CPLEX","SCIP","LPSolve"}));
		s.quiet=true; // disable solver output

		// declare the Variables
        ScaLP::Variable x = ScaLP::newIntegerVariable("x"); // x is free
        // ScaLP::Variable x = ScaLP::newIntegerVariable("x",-ScaLP::INF(),ScaLP::INF()); // alternate
        ScaLP::Variable y = ScaLP::newRealVariable("y",12.5,26);

		// Set the Objective
        ScaLP::Term t = x;
        ScaLP::Objective o = ScaLP::maximize(t);
		s.setObjective(o); // alternate: s<<o;

		// print objective
		std::cout << "Objective: " << o << std::endl;

		// add the Constraints
        ScaLP::Constraint c1 = x+y<=30;
        ScaLP::Constraint c2 = 5<=x<=30;
		s<<c1<<c2;
		//s.addConstraint(c1); // alternate

		// write or print a LP-Format-representation
		//std::cout << s.showLP() << std::endl;
		s.writeLP("simple_gurobi.lp");

		// Try to solve
        ScaLP::status stat = s.solve();

		// print results
		std::cout << "The result is " << stat << std::endl;
        if(stat==ScaLP::status::OPTIMAL || stat==ScaLP::status::FEASIBLE)
		{
          ScaLP::Result r = s.getResult();
		  std::cout << r << std::endl;
		}
/*
        s.reset();

        ScaLP::Variable u = ScaLP::newIntegerVariable("u"); // u is free
        // ScaLP::Variable u = ScaLP::newIntegerVariable("u",-ScaLP::INF(),ScaLP::INF()); // alternate
        ScaLP::Variable v = ScaLP::newRealVariable("v",12.5,26);

		// Set the Objective
        ScaLP::Term w = u;
        ScaLP::Objective oo = ScaLP::maximize(w);
		s.setObjective(oo); // alternate: s<<o;

		// print objective
		std::cout << "Objective: " << oo << std::endl;

		// add the Constraints
        ScaLP::Constraint c11 = u+v<=30;
        ScaLP::Constraint c22 = 5<=u<=30;
		s<<c11<<c22;
		//s.addConstraint(c1); // alternate

		// write or print a LP-Format-representation
		//std::cout << s.showLP() << std::endl;
		s.writeLP("simple_gurobi.lp");

		// Try to solve
        ScaLP::status statt = s.solve();

		// print results
		std::cout << "The result is " << statt << std::endl;
        if(stat==ScaLP::status::OPTIMAL || stat==ScaLP::status::FEASIBLE)
		{
          ScaLP::Result r = s.getResult();
		  std::cout << r << std::endl;
		}
*/
#endif //HAVE_SCALP

}

BitHeapILPCompression::~BitHeapILPCompression()
{
        cleanUp();
}

int BitHeapILPCompression::getMaxStageCount(int maxHeight)
{
        int height=2;
        int stageCount=0;
        cout << "maxHeight is " << maxHeight << endl;
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
    cout << "getExternalStageCount is " << getExternalStageCount << endl;
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

#ifdef HAVE_SCALP
    cout << "before setting pointer " << endl;
    string usedILPSolver = UserInterface::ilpSolver;
    if(usedILPSolver.compare("Gurobi") == 0 && usedILPSolver.compare("CPLEX") == 0 && usedILPSolver.compare("SCIP") == 0 && usedILPSolver.compare("LPSolve") == 0){
        THROWERROR("ilpSolver " << usedILPSolver << " unknown!");
    }
    problemSolver = new ScaLP::Solver(ScaLP::newSolverDynamic({usedILPSolver}));
    //problemSolver = new ScaLP::Solver(ScaLP::newSolverDynamic({"Gurobi","CPLEX","SCIP","LPSolve"}));
    solvers.push_back(problemSolver);
    //problemSolver = solvers[solvers.size() - 1];
    //problemSolver = ScaLP::Solver(ScaLP::newSolverDynamic({"CPLEX","SCIP","LPSolve"}));
    //problemSolver.reset();

#else
    /* initialize SCIP */
    SCIP_CALL( SCIPcreate(&scip) );
    SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
    SCIP_CALL( SCIPcreateProbBasic(scip, "example") );
    SCIP_VAR* tmpvar;
#endif //HAVE_SCALP


    //create variables:
    compCountVars.clear();
    compCountVars.resize(noOfStages_+1);
    cout << "noOfStages_ is " << noOfStages_ << endl << endl << endl;
    for(unsigned s=0; s < compCountVars.size(); s++)
    {
        compCountVars[s].resize(noOfCompressors);
    }



    stageVars.resize(noOfStages_+1);
    for(unsigned s=0; s < stageVars.size(); s++)
    {
        stringstream varName;
        varName << "D_" << s;
#ifdef HAVE_SCALP
        ScaLP::Variable tempD = ScaLP::newBinaryVariable(varName.str());
        stageVars[s] = tempD;
#else
        SCIP_CALL( SCIPcreateVarBasic(scip, &tmpvar, varName.str().c_str(), 0, 1, 0, SCIP_VARTYPE_BINARY ) );
        stageVars[s] = tmpvar;
        SCIP_CALL( SCIPaddVar(scip, tmpvar) );
#endif // HAVE_SCALP
    }


    for(unsigned s=0; s < compCountVars.size() - 1; s++)
    {
        for(unsigned e=0; e < compCountVars[s].size(); e++)
        {
            for(unsigned c=0; c < newBits[0].size()+s*(compOutputWordSizeMax-1); c++)
            {
                //cout << "compCountVars[s][e][c] is being processed " << s << " " << e << " " << c << " areaCost: " << (*possibleCompressors_)[e]->areaCost << endl;
                stringstream varName;

                varName << "k_" << s << "_" << e << "_" << c;
#ifdef HAVE_SCALP
                //in scalp, we put the area (the factor) when setting the objective. therefore the conditions are not necessary.
                ScaLP::Variable tempK = ScaLP::newIntegerVariable(varName.str(), 0, ScaLP::INF());
                compCountVars[s][e].push_back(tempK);
                if(useVariableCompressors){
                    if(e >= possibleCompressors_->size()){
                        unsigned offset = possibleCompressors_->size();
                        if(   ((e - offset) % 3 == 2)   && c == (newBits[0].size()+s*(compOutputWordSizeMax-1)) - 1){
                            stringstream varName;
                            varName << "k_" << s << "_" << e << "_" << c+1;
                            ScaLP::Variable tempK = ScaLP::newIntegerVariable(varName.str(), 0, ScaLP::INF());
                            compCountVars[s][e].push_back(tempK);
                        }
                    }
                }
#else
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
                if(useVariableCompressors){
                    if(e >= possibleCompressors_->size())
                    {
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
#endif //HAVE_SCALP

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
#ifdef HAVE_SCALP
            ScaLP::Variable tempN = ScaLP::newIntegerVariable(varName.str(), 0, ScaLP::INF());
            columnBitCountVars[s - 1].push_back(tempN);
#else
            SCIP_CALL( SCIPcreateVarBasic(scip, &tmpvar, varName.str().c_str(), 0.0, SCIPinfinity(scip), 0, SCIP_VARTYPE_INTEGER) );
            columnBitCountVars[s - 1].push_back(tmpvar);
            SCIP_CALL( SCIPaddVar(scip, tmpvar) );
#endif //HAVE_SCALP
        }
    }
    //cout << "declaration of N done" << endl;
    //now create U and store them in newBitsCountVars
    newBitsCountVars.clear();
    newBitsCountVars.resize(noOfStages_+1);
    for(unsigned s=0; s < newBitsCountVars.size(); s++){
#ifdef HAVE_SCALP
        vector<ScaLP::Variable> tempVec;
#else
        vector<SCIP_VAR*> tempVec;
#endif  //HAVE_SCALP
        for(unsigned c=0; c < (((unsigned)bh_->bits.size())+s*(compOutputWordSizeMax-1)); c++){
            stringstream varName;
            varName << "U_" << s << "_" << c;
#ifdef HAVE_SCALP
            ScaLP::Variable tempU = ScaLP::newIntegerVariable(varName.str(), -ScaLP::INF(), ScaLP::INF());
            tempVec.push_back(tempU);
#else
            SCIP_CALL( SCIPcreateVarBasic(scip, &tmpvar, varName.str().c_str(), -SCIPinfinity(scip), SCIPinfinity(scip), 0, SCIP_VARTYPE_INTEGER) );
            tempVec.push_back(tmpvar);
            SCIP_CALL( SCIPaddVar(scip, tmpvar) );
#endif //HAVE_SCALP
        }
        newBitsCountVars[s] = tempVec;
    }

    //cout << "declaration of U done" << endl;
    //add constraints:
#ifdef HAVE_SCALP
    //defining the objective
    ScaLP::Term objectiveSum;
    cout << compCountVars.size() << endl;
    for(unsigned int s = 0; s < compCountVars.size() - 1; s++){
        for(unsigned int e = 0; e < compCountVars[s].size(); e++){
            for(unsigned int c = 0; c < compCountVars[s][e].size(); c++){
                if(!useVariableCompressors){
                    objectiveSum = objectiveSum + (*possibleCompressors_)[e]->areaCost * compCountVars[s][e][c];
                }
                else{
                    if(e < possibleCompressors_->size()){
                        objectiveSum = objectiveSum + (*possibleCompressors_)[e]->areaCost * compCountVars[s][e][c];
                    }
                    else{
                        //cout << "e is " << e << " and variableBCompressors[e - possibleCompressors_->size()].areaCost is " << variableBCompressors[e - possibleCompressors_->size()].areaCost << endl;
                        objectiveSum = objectiveSum + variableBCompressors[e - possibleCompressors_->size()].areaCost * compCountVars[s][e][c];
                    }
                }
            }
        }
    }
    cout << "after setting objective" << endl;
    ScaLP::Objective obj = ScaLP::minimize(objectiveSum);
    problemSolver->setObjective(obj);
    //problemSolver->writeLP("compressorTree.lp");
#else
    SCIP_CONS* tmpcons; //needed for scip constraints
#endif //HAVE_SCALP

    //add constraint C0: fill all the U's with values


    //cout << newBitsCountVars.size() << " " << newBits.size() << endl;
    for(unsigned i = 0; i < newBitsCountVars.size(); i++){
#ifdef HAVE_SCALP
        vector<ScaLP::Variable> tempVectorCons = newBitsCountVars.at(i);
#else
        vector<SCIP_VAR*> tempVectorCons = newBitsCountVars.at(i);
#endif //HAVE_SCALP
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
#ifdef HAVE_SCALP
            ScaLP::Constraint tempConstraint = tempVectorCons.at(j) - value == 0;
            //tempConstraint.setName(consName.str());
            tempConstraint.name = consName.str();
            problemSolver->addConstraint(tempConstraint);

#else
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &tmpcons, consName.str().c_str(), 0, NULL, NULL, value, value) );
            SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons,tempVectorCons.at(j), 1.0) );
            SCIP_CALL( SCIPaddCons(scip, tmpcons) );
#endif //HAVE_SCALP
        }

    }
    //cout << compCountVars.size() << " is compCountVars.size()" << endl;
    if(compCountVars[0][0][0] != NULL){
        //cout << "0-0-0 is full" << endl;
    }
    else{
        //cout << "0-0-0 is not full" << endl;
    }
    cout << "declaration of C0 done" << endl;
    const int LARGE_NUMBER = 10000;
//  int LARGE_NUMBER = bh_->getMaxHeight()+1; //!!!









    //add constraints C1 and C2:


    for(unsigned s=0; s < compCountVars.size()-1; s++){



        //add constraint C1:
        for(int c=0; c < ((int) (newBits[0].size()+s*(compOutputWordSizeMax-1))); c++){
            stringstream consName;
            consName << "C1_" << s << "_" << c;
            //cout << "s = " << s << " and c = " << c << endl;
#ifdef HAVE_SCALP
            ScaLP::Term c1Term;
#else
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &tmpcons, consName.str().c_str(), 0, NULL, NULL, 0, SCIPinfinity(scip)) );
#endif // HAVE_SCALP
            for(unsigned e=0; e < possibleCompressors_->size(); e++){
                for(int ce=0; ce < (int) (*possibleCompressors_)[e]->height.size(); ce++){
                    if(c-ce >= 0){
                        //cout << "call s#e#c-ce: " << s << e << c-ce << endl;
#ifdef HAVE_SCALP
                        c1Term = c1Term + (*possibleCompressors_)[e]->height[(*possibleCompressors_)[e]->height.size()-ce-1] * compCountVars[s][e][c-ce];
#else
                        SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons, compCountVars[s][e][c-ce] , (*possibleCompressors_)[e]->height[(*possibleCompressors_)[e]->height.size()-ce-1]) );
#endif // HAVE_SCALP
                        //cout << "after call" << endl;
                    }
                }
            }
            if(useVariableCompressors){
                for(unsigned e=possibleCompressors_->size(); e < possibleCompressors_->size() + variableBCompressors.size(); e++){
                    for(int ce=0; ce < (int) variableBCompressors[e - possibleCompressors_->size()].height.size(); ce++){
                        if(c-ce >= 0){
                            //cout << "call s#e#c-ce: " << s << e << c-ce << endl;
#ifdef HAVE_SCALP
                            c1Term = c1Term + variableBCompressors[e - possibleCompressors_->size()].height[variableBCompressors[e - possibleCompressors_->size()].height.size()-ce-1] * compCountVars[s][e][c-ce];
#else
                            SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons, compCountVars[s][e][c-ce] , variableBCompressors[e - possibleCompressors_->size()].height[variableBCompressors[e - possibleCompressors_->size()].height.size()-ce-1]) );
#endif //HAVE_SCALP
                            //cout << "after call" << endl;
                        }
                    }
                }
            }
            //cout << "after loop" << endl;
#ifdef HAVE_SCALP
            if(s != 0){
                //first stage is in U therefore N starts with s-1
                c1Term = c1Term - columnBitCountVars[s - 1][c];
            }
            c1Term = c1Term - newBitsCountVars[s][c];
            c1Term = c1Term + LARGE_NUMBER * stageVars[s];

            ScaLP::Constraint c1Constraint = c1Term >= 0;
            //c1Constraint.setName(consName.str());
            c1Constraint.name = consName.str();
            problemSolver->addConstraint(c1Constraint);
            //problemSolver->writeLP("compressorTreec1.lp");


#else
            if(s != 0){
                //first stage is in U therefore N starts with s-1
                SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons, columnBitCountVars[s - 1][c] , -1) );
            }
            //subtracting U as well
            SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons, newBitsCountVars[s][c] , -1) );
            SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons, stageVars[s] , LARGE_NUMBER) );
            SCIP_CALL( SCIPaddCons(scip, tmpcons) );
#endif //HAVE_SCALP
        }
        cout << "wrote c1Constraint" << endl;
        //cout << "filling of C1 at stage " << s << " done" << endl;

        //-----------------------------------------------------------------------

        //add constraint C2:
        for(int c=0; c < ((int) (newBits[0].size()+(s+1)*(compOutputWordSizeMax-1))); c++){
            stringstream consName;
            consName << "C2_" << s+1 << "_" << c;
#ifdef HAVE_SCALP
            ScaLP::Term c2Term;
#else
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &tmpcons, consName.str().c_str(), 0, NULL, NULL, 0, 0) );
#endif
            for(unsigned e=0; e < possibleCompressors_->size(); e++){
                for(int ce=(int) (*possibleCompressors_)[e]->outputs.size()-1; ce >= 0; ce--){
                    if((c-ce >= 0) && (c-ce < ((int) (newBits[0].size()+s*(compOutputWordSizeMax-1))))){
#ifdef HAVE_SCALP
                        c2Term = c2Term + (*possibleCompressors_)[e]->outputs[(*possibleCompressors_)[e]->outputs.size()-ce-1] * compCountVars[s][e][c-ce];
#else
                        SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons, compCountVars[s][e][c-ce] , (*possibleCompressors_)[e]->outputs[(*possibleCompressors_)[e]->outputs.size()-ce-1]) );
#endif //HAVE_SCALP
                    }
                }
            }
            if(useVariableCompressors){
                unsigned offset = possibleCompressors_->size();
                for(unsigned e=offset; e < offset + variableBCompressors.size(); e++){
                    for(int ce=(int) variableBCompressors[e - offset].outputs.size()-1; ce >= 0; ce--){
                        if((c-ce >= 0) && (c-ce < ((int) (newBits[0].size()+s*(compOutputWordSizeMax-1))))){
#ifdef HAVE_SCALP
                            c2Term = c2Term + variableBCompressors[e - offset].outputs[variableBCompressors[e - offset].outputs.size()-ce-1] * compCountVars[s][e][c-ce];
#else
                            SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons, compCountVars[s][e][c-ce] , variableBCompressors[e - offset].outputs[variableBCompressors[e - offset].outputs.size()-ce-1]) );
#endif
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
#ifdef HAVE_SCALP
                            c2Term = c2Term + variableBCompressors[e - offset].outputs[tempOutput] * compCountVars[s][e][c];
#else
                            SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons, compCountVars[s][e][c] , variableBCompressors[e - offset].outputs[tempOutput]) );
#endif //HAVE_SCALP
                        }

                    }
                }
            }
#ifdef HAVE_SCALP
            c2Term = c2Term - columnBitCountVars[s][c];
            ScaLP::Constraint c2Constraint = c2Term == 0;
            //c2Constraint.setName(consName.str());
            c2Constraint.name = consName.str();
            problemSolver->addConstraint(c2Constraint);
            //problemSolver->writeLP("compressorTreec2.lp");
#else
            SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons, columnBitCountVars[s][c] , -1) );
            SCIP_CALL( SCIPaddCons(scip, tmpcons) );
#endif //HAVE_SCALP

        }
        //cout << "filling of C2 at stage " << s << " done" << endl;
    }

    cout << "C1 and C2 totally done" << endl;


    //------------------------------------------------------------

    //add constraint C3:
    for(unsigned s=0; s < compCountVars.size(); s++){
        for(unsigned c=0; c < newBitsCountVars[s].size(); c++){
            stringstream consName;
            consName << "C3_" << s << "_" << c;
#ifdef HAVE_SCALP
            ScaLP::Term c3Term;
            if(s != 0){
                c3Term = c3Term + columnBitCountVars[s - 1][c];
            }
            c3Term = c3Term + newBitsCountVars[s][c];
            for(unsigned z = s + 1 ; z < compCountVars.size(); z++){
                //make sure that that in the follwoing stages all the U's are empty.
                //this is done by adding all of the U-variables of later stages with the factor of 4. Therefore e.g. if s = 3 -> D_3 = 1 -> in constraint C3_3_c: 1000 + 4 * U_later <= 1002 -> if at least one U_later is >0, then constraint is not met
                c3Term = c3Term + 4 * newBitsCountVars[z][c];
            }
            c3Term = c3Term + LARGE_NUMBER * stageVars[s];
            ScaLP::Constraint c3Constraint = c3Term <= LARGE_NUMBER + 2;
            //c3Constraint.setName(consName.str());
            c3Constraint.name = consName.str();
            problemSolver->addConstraint(c3Constraint);
            //problemSolver->writeLP("compressorTreec3.lp");

#else
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
#endif //HAVE_SCALP
        }

    }

    //cout << "C3 filling done" << endl;

    //-----------------------------------------------------------

    //add constraint C4:
#ifdef HAVE_SCALP
    ScaLP::Term c4Term;
    for(unsigned s=0; s < stageVars.size(); s++){
        c4Term = c4Term + stageVars[s];
    }
    ScaLP::Constraint c4Constraint = c4Term - 1 == 0;
    //c4Constraint.setName("C4");
    c4Constraint.name = "C4";
    problemSolver->addConstraint(c4Constraint);


#else
    SCIP_CALL( SCIPcreateConsBasicLinear(scip, &tmpcons, "C4", 0, NULL, NULL, 1, 1) );
    for(unsigned s=0; s < stageVars.size(); s++){
        SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons, stageVars[s] , 1) );
    }
    SCIP_CALL( SCIPaddCons(scip, tmpcons) );
#endif

    //-----------------------------------------------------------

    //add constraint C4_1: setting up the 2-bit-adder stage


    if(useFixedStageCount){
#ifdef HAVE_SCALP
        ScaLP::Constraint c4_1Constraint = stageVars[noOfStages_] - 1 == 0;
        //c4_1Constraint.setName("C4_1");
        c4_1Constraint.name = "C4_1";
        problemSolver->addConstraint(c4_1Constraint);
        //problemSolver->writeLP("compressorTreec4_1.lp");
#else
        SCIP_CALL( SCIPcreateConsBasicLinear(scip, &tmpcons, "C4_1", 0, NULL, NULL, 1, 1) );
        SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons, stageVars[noOfStages_], 1) );
        SCIP_CALL( SCIPaddCons(scip, tmpcons) );
#endif
    }

    //cout << "C4 filling done" << endl;


    //----------------------------------------------------------

    if(useVariableCompressors){

        cout << "in useVariableCompressors " << endl << endl << endl;
        //assume that the order in variableBCompressors is low - middle -high ( (3;1) - (2;1) - (0;1) )

        unsigned offset = possibleCompressors_->size();

        for(unsigned s = 0; s < compCountVars.size() - 1; s++){

			for(unsigned e = 0; e < variableBCompressors.size(); e += 3){	//C5 only needed for every variable Compressor, not every variable Basic Compressor

				for(int c=0; c < ((int) (newBits[0].size()+s*(compOutputWordSizeMax-1))); c++){
					//C5:

					//format : k_high_c=i+1 + k_mid_c=i+1 = k_mid_c=i + k_low_c=i
					//setting k_high_0 and k_mid_0 to zero is done after this loop

					stringstream consName;
					consName << "C5_" << s << "_" << e << "_" << c;
#ifdef HAVE_SCALP
                    ScaLP::Term c5Term;
                    c5Term = c5Term + compCountVars[s][offset + e + 2][c + 1];  //high
                    if(c != ((int) (newBits[0].size()+s*(compOutputWordSizeMax-1))) - 1){
                        c5Term = c5Term + compCountVars[s][offset + e + 1][c + 1];  //middle
                    }
                    c5Term = c5Term - compCountVars[s][offset + e + 0][c];  //low
                    c5Term = c5Term - compCountVars[s][offset + e + 1][c];  //middle
                    ScaLP::Constraint c5Constraint = c5Term == 0;
                    //c5Constraint.setName(consName.str());
                    c5Constraint.name = consName.str();
                    problemSolver->addConstraint(c5Constraint);
                    //problemSolver->writeLP("compressorTreec5.lp");
#else
					SCIP_CALL( SCIPcreateConsBasicLinear(scip, &tmpcons, consName.str().c_str(), 0, NULL, NULL, 0, 0) );

					SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons,compCountVars[s][offset + e + 2][c + 1], 1.0) );
					if(c != ((int) (newBits[0].size()+s*(compOutputWordSizeMax-1))) - 1){   //if it is not the highest column, in the next higher column can be a middle compressor
						SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons,compCountVars[s][offset + e + 1][c + 1], 1.0) );
					}
					SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons,compCountVars[s][offset + e + 0][c], -1.0) );
					SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons,compCountVars[s][offset + e + 1][c], -1.0) );            //new: subtract the middle compressors as well

					SCIP_CALL( SCIPaddCons(scip, tmpcons) );
#endif //HAVE_SCALP
				}

				//make sure that the high compressor and middle compressor are not used in the lowest column.

				stringstream consNameB;
				consNameB << "C5b_" << e << "_" << s;
#ifdef HAVE_SCALP
                ScaLP::Constraint c5bConstraint;
                c5bConstraint = compCountVars[s][offset + e + 2][0] + compCountVars[s][offset + e + 1][0] == 0;
                //c5bConstraint.setName(consNameB.str());
                c5bConstraint.name = consNameB.str();
                problemSolver->addConstraint(c5bConstraint);
                //problemSolver->writeLP("compressorTreec5_b.lp");

#else
				SCIP_CALL( SCIPcreateConsBasicLinear(scip, &tmpcons, consNameB.str().c_str(), 0, NULL, NULL, 0, 0) );
				SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons,compCountVars[s][offset + e + 2][0], 1.0) );
				SCIP_CALL( SCIPaddCoefLinear(scip, tmpcons,compCountVars[s][offset + e + 1][0], 1.0) );
				SCIP_CALL( SCIPaddCons(scip, tmpcons) );
#endif //HAVE_SCALP
			}
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

/*  //we dont use zero stages any longer?


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

*/

#ifdef HAVE_SCALP
    //set objective     //TODO remove
    ScaLP::Term objectiveTerm;

    for(unsigned s=0; s < compCountVars.size(); s++)
    {
        for(unsigned e=0; e < compCountVars[s].size(); e++)
        {
            for(unsigned c=0; c < compCountVars[s][e].size(); c++)
            {
                if(!useVariableCompressors){
                    objectiveTerm = objectiveTerm + (*possibleCompressors_)[e]->areaCost * compCountVars[s][e][c];
                }
                else{
                    if(e < possibleCompressors_->size()){
                        objectiveTerm = objectiveTerm + (*possibleCompressors_)[e]->areaCost * compCountVars[s][e][c];
                    }
                    else{
                        objectiveTerm = objectiveTerm + variableBCompressors[e - possibleCompressors_->size()].areaCost * compCountVars[s][e][c];
                    }
                }
            }
        }
    }

    ScaLP::Objective o = ScaLP::minimize(objectiveTerm);
    //problemSolver->setObjective(o);
    cout << "before writeLP" << endl;
    problemSolver->writeLP("compressorTree.lp");
    cout << "after writeLP" << endl;

#endif


    //cout << "end of SCIP problem description" << endl;
    return 0;
}





int BitHeapILPCompression::writeProblem(std::string filename){
    const char *fname=NULL;;
    if(!filename.empty()){
        fname = filename.c_str();
    }
    //print model in LP format:
#ifdef HAVE_SCALP
    cout << "before writeLP" << endl;
    //problemSolver->writeLP(fname);
    //problemSolver->showLP();
    cout << "after writeLP" << endl;
#else
    SCIP_RESULT result;
    SCIP_CALL(SCIPwriteLp(scip, NULL, fname, FALSE, SCIP_OBJSENSE_MINIMIZE, 1.0, 0.0, SCIPgetVars(scip), SCIPgetNVars(scip), SCIPgetNBinVars(scip), SCIPgetNIntVars(scip), SCIPgetNImplVars(scip), SCIPgetNContVars(scip), SCIPgetConss(scip), SCIPgetNOrigConss(scip), &result));
#endif
    return 0;
}

bool BitHeapILPCompression::solve(){
    //set a time limit:
    int timeout;
    timeout = UserInterface::ilpTimeout;
    //timeout = bh_->getOp()->getTarget()->ilpTimeout();
#ifdef HAVE_SCALP
    problemSolver->timeout = (int) timeout;
    cout << " timeout is set to " << timeout;
    problemSolver->threads = 1;  //not needed with scip because scip uses only one thread
    problemSolver->quiet = false;

#else
    SCIP_CALL( SCIPsetRealParam(scip,"limits/time",timeout) );
#endif //HAVE_SCALP
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

#ifdef HAVE_SCALP

    cout << "backend is " << problemSolver->getBackendName() << endl;
    ScaLP::status stat = problemSolver->solve();
    cout << "after solving in scalp" << endl;
    cout << "status is " << stat << endl;

    if(stat == ScaLP::status::INFEASIBLE_OR_UNBOUND || stat == ScaLP::status::INFEASIBLE || stat == ScaLP::status::UNBOUND){
        cout << "in unbound/infeasible - mode" << endl;
        if(!useFixedStageCount){
            THROWERROR("No optimal solution found (problem infeasible or unbounded!)");
        }
        else{
            REPORT(INFO, "No optimal solution for current stage count found (problem infeasible or unbound!).\n We will try it again with more stages");
            infeasible = true;
        }
    }
    else if(stat == ScaLP::status::OPTIMAL || stat == ScaLP::status::FEASIBLE || (stat == ScaLP::status::TIMEOUT && problemSolver->getResult().objectiveValue > 0)){
        infeasible = false;
    }

    ScaLP::Result result = problemSolver->getResult();

    if(infeasible && !useFixedStageCount){
        THROWERROR("No feasible solution found within the ILP timeout of " << timeout << " seconds");
    }
    else if(!infeasible){
        double area = problemSolver->getResult().objectiveValue;
        preReductionAreaCost += area;
        REPORT(DEBUG, "The size of the complete Solution is " << preReductionAreaCost << " LUTs");
        costOfCurrentSolution = area;
    }
#else
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

#endif //HAVE_SCALP
    cout << "finished setting infeasible" << endl;
    /*
    SCIP_Bool feasible;
    SCIP_CALL( SCIPcheckSol ( scip,sol,true,true,true,true,&feasible) );

    if(!feasible)
    {
            cerr << "No feasible solution found within the ILP timeout of " << timeout << " seconds" << endl;
            exit(-1);
    }
    */

#ifdef HAVE_SCALP
    if(!infeasible){
        if(!useHeuristic){
            solution.resize(compCountVars.size());
        }
        cout << "trying to read variables" << endl;

        //first try: go through compCountVars, get the values from each variable
        for(unsigned int s = 0; s < compCountVars.size(); s++){
            for(unsigned int e = 0; e < compCountVars[s].size(); e++){
                for(unsigned int c = 0; c < compCountVars[s][e].size(); c++){
                    double tempValue = result.values[compCountVars[s][e][c]];
                    tempValue += 0.00001;   //float to integer: add small value and cast
                    int integerValue = (int) tempValue;
                    if(integerValue > 0)
                    {
                        for(unsigned int k = 0; k < (unsigned int) integerValue; k++){
                            //cout << "in stages " << s << " and column " << c << " is compressor " << e << endl;
                            solution[s].push_back(pair<int,int>(e,c));
                        }
                    }

                }
            }
        }
    }
    cout << " finished writing scalp solution" << endl;
#else

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
#endif //HAVE_SCALP
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
#ifndef HAVE_SCALP
    SCIP_SOL* heuSol;

    //cout << "in passHeuristicSolution" << endl;

    for(unsigned int i = 0; i < heuristicSolutions.size(); i++){
        REPORT(DEBUG, "passing heuristic solution number " << i);
        SCIP_CALL( SCIPcreateSol(scip, &heuSol, NULL) );
        computeCompressorCount(i);      //computes also heuristicN

        //U: is the same for every solution
        for(unsigned int s = 0; s < newBits.size(); s++){
            for(unsigned int c = 0; c < newBits[s].size(); c++){

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

#endif //ifndef HAVE_SCALP

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
#ifndef HAVE_SCALP
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
#endif //HAVE_SCALP
    }

    int BitHeapILPCompression::cleanUp()
    {
        REPORT(DEBUG, "cleaning up SCIP... ");
        //clean up:
        cout << "!!!skipping cleanup for debug reasons!!!" << endl;
        return 0; //!!!
#ifdef HAVE_SCALP
        //deleting solvers
        for(int i = solvers.size() -1; i >= 0; i--){
            delete solvers[i];
        }

        return 0;
#else
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

#endif //HAVE_SCALP
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
