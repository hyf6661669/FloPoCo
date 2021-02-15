#include "IntMult/TilingStrategyOptimalILP.hpp"
#include "IntMultiplier.hpp"
#include "BaseMultiplierLUT.hpp"
#include "MultiplierTileCollection.hpp"

using namespace std;
namespace flopoco {

TilingStrategyOptimalILP::TilingStrategyOptimalILP(
		unsigned int wX_,
		unsigned int wY_,
		unsigned int wOut_,
		bool signedIO_,
		BaseMultiplierCollection* bmc,
		base_multiplier_id_t prefered_multiplier,
		float occupation_threshold,
		int maxPrefMult,
        MultiplierTileCollection mtc_,
        unsigned guardBits,
        unsigned keepBits,
        unsigned long long errorBudget,
        unsigned long long &centerErrConstant):TilingStrategy(
			wX_,
			wY_,
			wOut_,
			signedIO_,
			bmc),
		max_pref_mult_ {maxPrefMult},
		occupation_threshold_{occupation_threshold},
		tiles{mtc_.MultTileCollection},
        guardBits{guardBits},
        keepBits{keepBits},
        errorBudget{errorBudget},
        centerErrConstant{centerErrConstant}
	{
	    cout << errorBudget << endl;
	    cout << this->errorBudget << endl;
	    cout << "guardBits " << guardBits << " keepBits " << keepBits << endl;
        for(auto &p:tiles)
        {
                cout << p->getLUTCost(0, 0, wX, wY) << " " << p->getType() << endl;
        }
	}

void TilingStrategyOptimalILP::solve()
{

#ifndef HAVE_SCALP
    throw "Error, TilingStrategyOptimalILP::solve() was called but FloPoCo was not built with ScaLP library";
#else
    cout << "using ILP solver " << target->getILPSolver() << endl;
    solver = new ScaLP::Solver(ScaLP::newSolverDynamic({target->getILPSolver(),"Gurobi","CPLEX","SCIP","LPSolve"}));
    solver->timeout = target->getILPTimeout();

    constructProblem();

    // Try to solve
    cout << "starting solver, this might take a while..." << endl;
    solver->quiet = false;
    ScaLP::status stat = solver->solve();

    // print results
    cerr << "The result is " << stat << endl;
    //cerr << solver->getResult() << endl;
    ScaLP::Result res = solver->getResult();

    cout << "centerErrConstant was: " << centerErrConstant << endl;
    unsigned long long new_constant = 0;
    double total_cost = 0;
    int dsp_cost = 0, own_lut_cost=0;
    for(auto &p:res.values)
    {
        if(p.second > 0.5){     //parametrize all multipliers at a certain position, for which the solver returned 1 as solution, to flopoco solution structure
            std::string var_name = p.first->getName();
            if(var_name.substr(1,1) == "d"){
                int mult_id = stoi(var_name.substr(2,dpS));
                int x_negative = (var_name.substr(2+dpS,1).compare("m") == 0)?1:0;
                int m_x_pos = stoi(var_name.substr(2+dpS+x_negative,dpX)) * ((x_negative)?(-1):1);
                int y_negative = (var_name.substr(2+dpS+x_negative+dpX,1).compare("m") == 0)?1:0;
                int m_y_pos = stoi(var_name.substr(2+dpS+dpX+x_negative+y_negative,dpY)) * ((y_negative)?(-1):1);
                cout << "is true:  " << setfill(' ') << setw(dpY) << mult_id << " " << setfill(' ') << setw(dpY) << m_x_pos << " " << setfill(' ') << setw(dpY) << m_y_pos << " cost: " << setfill(' ') << setw(5) << tiles[mult_id]->getLUTCost(m_x_pos, m_y_pos, wX, wY) << std::endl;

                total_cost += (double)tiles[mult_id]->getLUTCost(m_x_pos, m_y_pos, wX, wY);
                own_lut_cost += tiles[mult_id]->ownLUTCost(m_x_pos, m_y_pos, wX, wY);
                dsp_cost += (double)tiles[mult_id]->getDSPCost();
                auto coord = make_pair(m_x_pos, m_y_pos);
                solution.push_back(make_pair(tiles[mult_id]->getParametrisation().tryDSPExpand(m_x_pos, m_y_pos, wX, wY, signedIO), coord));
            }
            if(var_name.substr(0,1) == "c"){
                int c_id = stoi(var_name.substr(1,dpC));
                new_constant |= (1ULL<<c_id);
                cout << var_name << " pos " << c_id << " dpC " << dpC << endl;
            }
        }
    }
    cout << "Total LUT cost:" << total_cost <<std::endl;
    cout << "Own LUT cost:" << own_lut_cost <<std::endl;
    cout << "Total DSP cost:" << dsp_cost <<std::endl;

    centerErrConstant = new_constant;
    cout << "centerErrConstant now is: " << centerErrConstant << endl;
/*
    solution.push_back(make_pair(tiles[1]->getParametrisation().tryDSPExpand(0, 0, wX, wY, signedIO), make_pair(0, 0)));
    solution.push_back(make_pair(tiles[0]->getParametrisation().tryDSPExpand(16, 0, wX, wY, signedIO), make_pair(16, 0)));
    solution.push_back(make_pair(tiles[0]->getParametrisation().tryDSPExpand(32, 0, wX, wY, signedIO), make_pair(32, 0)));
    solution.push_back(make_pair(tiles[0]->getParametrisation().tryDSPExpand(0, 24, wX, wY, signedIO), make_pair(0, 24)));
    solution.push_back(make_pair(tiles[0]->getParametrisation().tryDSPExpand(16, 24, wX, wY, signedIO), make_pair(16, 24)));
    solution.push_back(make_pair(tiles[0]->getParametrisation().tryDSPExpand(32, 24, wX, wY, signedIO), make_pair(32, 24)));
    solution.push_back(make_pair(tiles[0]->getParametrisation().tryDSPExpand(48, 24, wX, wY, signedIO), make_pair(48, 24)));
    solution.push_back(make_pair(tiles[0]->getParametrisation().tryDSPExpand(16, 48, wX, wY, signedIO), make_pair(16, 48)));
    solution.push_back(make_pair(tiles[0]->getParametrisation().tryDSPExpand(32, 48, wX, wY, signedIO), make_pair(32, 48)));
*/
//    solution.push_back(make_pair(tiles[1]->getParametrisation(), make_pair(0, 0)));

#endif
}

#ifdef HAVE_SCALP
void TilingStrategyOptimalILP::constructProblem()
{
	bool performOptimalTruncation = true;

    cout << "constructing problem formulation..." << endl;
    wS = tiles.size();


    //Assemble cost function, declare problem variables
    cout << "   assembling cost function, declaring problem variables..." << endl;
    ScaLP::Term obj;
    prodWidth = IntMultiplier::prodsize(wX, wY, signedIO, signedIO);
    int x_neg = 0, y_neg = 0;
    for(int s = 0; s < wS; s++){
        x_neg = (x_neg < (int)tiles[s]->wX())?tiles[s]->wX() - 1:x_neg;
        y_neg = (y_neg < (int)tiles[s]->wY())?tiles[s]->wY() - 1:y_neg;
    }
    int nx = wX-1, ny = wY-1, ns = wS-1, nc = prodWidth-wOut; dpX = 1; dpY = 1; dpS = 1; dpC = 1;   //calc number of decimal places, for var names
    nx = (x_neg > nx)?x_neg:nx;                                                                     //in case the extend in negative direction is larger
    ny = (y_neg > ny)?y_neg:ny;
    while (nx /= 10)
        dpX++;
    while (ny /= 10)
        dpY++;
    while (ns /= 10)
        dpS++;
    while (nc /= 10)
        dpC++;

    vector<vector<vector<ScaLP::Variable>>> solve_Vars(wS, vector<vector<ScaLP::Variable>>(wX+x_neg, vector<ScaLP::Variable>(wY+y_neg)));
    ScaLP::Term maxEpsTerm;
    __uint64_t sumOfPosEps = 0;
    // add the Constraints
    cout << "   adding the constraints to problem formulation..." << endl;
    for(int y = 0; y < wY; y++){
        for(int x = 0; x < wX; x++){
            stringstream consName;
            consName << "p" << setfill('0') << setw(dpX) << x << setfill('0') << setw(dpY) << y;            //one constraint for every position in the area to be tiled
            ScaLP::Term pxyTerm;
            for(int s = 0; s < wS; s++){					//for every available tile...
                for(int ys = 0 - tiles[s]->wY() + 1; ys <= y; ys++){					//...check if the position x,y gets covered by tile s located at position (xs, ys) = (x-wtile..x, y-htile..y)
                    for(int xs = 0 - tiles[s]->wX() + 1; xs <= x; xs++){
                        if(occupation_threshold_ == 1.0 && ((wX - xs) < (int)tiles[s]->wX() || (wY - ys) < (int)tiles[s]->wY())) break;
                        if(tiles[s]->shape_contribution(x, y, xs, ys, wX, wY, signedIO) == true){
                            if((wOut < (int)prodWidth) && ((xs+tiles[s]->wX()+ys+tiles[s]->wY()-2) < ((int)prodWidth-wOut-guardBits))) break;
                            if(tiles[s]->shape_utilisation(xs, ys, wX, wY, signedIO) >=  occupation_threshold_ ){
                                if(solve_Vars[s][xs+x_neg][ys+y_neg] == nullptr){
                                    stringstream nvarName;
                                    nvarName << " d" << setfill('0') << setw(dpS) << s << ((xs < 0)?"m":"") << setfill('0') << setw(dpX) << ((xs<0)?-xs:xs) << ((ys < 0)?"m":"")<< setfill('0') << setw(dpY) << ((ys<0)?-ys:ys) ;
                                    //std::cout << nvarName.str() << endl;
                                    ScaLP::Variable tempV = ScaLP::newBinaryVariable(nvarName.str());
                                    solve_Vars[s][xs+x_neg][ys+y_neg] = tempV;
                                    obj.add(tempV, (double)tiles[s]->getLUTCost(xs, ys, wX, wY));    //append variable to cost function
                                }
                                pxyTerm.add(solve_Vars[s][xs+x_neg][ys+y_neg], 1);
                            }
                        }
                    }
                }
            }
            ScaLP::Constraint c1Constraint;
            if(performOptimalTruncation == true && (wOut < (int)prodWidth) && ((x+y) < ((int)prodWidth-wOut))) {
                stringstream nvarName;
                nvarName << " b" << ((x < 0) ? "m" : "") << setfill('0') << setw(dpX) << ((x < 0) ? -x : x)
                         << ((y < 0) ? "m" : "") << setfill('0') << setw(dpY) << ((y < 0) ? -y : y);
                ScaLP::Variable tempV = ScaLP::newBinaryVariable(nvarName.str());
                maxEpsTerm.add(tempV, (-1) * ((long) 1 << (x + y)));
                maxEpsTerm.add((long) 1 << (x + y));

                c1Constraint = pxyTerm - tempV == 0;
            } else if(performOptimalTruncation == false && (wOut < (int)prodWidth) && ((x+y) < ((int)prodWidth-wOut-guardBits))){
                //c1Constraint = pxyTerm <= (bool)1;
            } else if(performOptimalTruncation == false && (wOut < (int)prodWidth) && ((x+y) == ((int)prodWidth-wOut-guardBits))){
                if((keepBits)?keepBits--:0){
                    c1Constraint = pxyTerm == (bool)1;
                    cout << "keepBit at" << x << "," << y << endl;
                } else {
                    cout << "NO keepBit at" << x << "," << y << endl;
                }
            } else {
                c1Constraint = pxyTerm == (bool)1;
            }

            c1Constraint.name = consName.str();
            solver->addConstraint(c1Constraint);
        }
    }

    //limit use of DSPs
    if(0 <= max_pref_mult_) {
        //check if DSP tiles are available
        bool nDSPTiles = false;
        for (int s = 0; s < wS; s++) {
            if (tiles[s]->getDSPCost()) {
                nDSPTiles = true;
                break;
            }
        }

        if (nDSPTiles) {
            cout << "   adding the constraint to limit the use of DSP-Blocks to " << max_pref_mult_ << " instances..."
                 << endl;
            stringstream consName;
            consName << "limDSP";
            ScaLP::Term pxyTerm;
            for (int y = 0 - 24 + 1; y < wY; y++) {
                for (int x = 0 - 24 + 1; x < wX; x++) {
                    for (int s = 0; s < wS; s++)
                        if (solve_Vars[s][x + x_neg][y + y_neg] != nullptr)
                            for (int c = 0; c < tiles[s]->getDSPCost(); c++)
                                pxyTerm.add(solve_Vars[s][x + x_neg][y + y_neg], 1);
                }
            }
            ScaLP::Constraint c1Constraint = pxyTerm <= max_pref_mult_;     //set max usage equ.
            c1Constraint.name = consName.str();
            solver->addConstraint(c1Constraint);
        }
    }

    //make shure the available precision is present in case of truncation
    if(performOptimalTruncation == true && (wOut < (int)prodWidth))
    {
        cout << "   multiplier is truncated by " << (int)prodWidth-wOut << " bits (err=" << (unsigned long)wX*(((unsigned long)1<<((int)wOut-guardBits))) << "), ensure sufficient precision..." << endl;
        cout << "   guardBits=" << guardBits << endl;
        /*int g, k;
        long long errorBudget;
        computeTruncMultParams((int)prodWidth-wOut, g, k, errorBudget);*/
        cout << "   g=" << guardBits << " k=" << keepBits << " errorBudget=" << errorBudget << " difference to conservative est: " << errorBudget-(long long)(((unsigned long)1)<<(prodWidth-(int)wOut-1)-1) << endl;
        //cout << sumOfPosEps << " " << ((unsigned long)1<<((int)prodWidth-wOut-1));
        //unsigned long maxErr = (prodWidth-(int)wOut-guardBits > prodWidth)?prodWidth:(prodWidth-(int)wOut-guardBits);
        //cout << "shift=" << maxErr << endl;
        //maxErr = ((unsigned long)((wX < wY) ? wX : wY)*(((unsigned long)1<<maxErr)));
        //unsigned long maxErr = ((unsigned long)1)<<(prodWidth-(int)wOut-1);

        stringstream nvarName;
        nvarName << "C";
        ScaLP::Variable Cvar = ScaLP::newIntegerVariable(nvarName.str());

        // C = 2^i*c_i + 2^(i+1)*c_(i+1)...
        ScaLP::Term cTerm;
        vector<ScaLP::Variable> cVars(guardBits-1);
        for(int i = 0; i < guardBits-1; i++){
            stringstream nvarName;
            nvarName << "c" << prodWidth-wOut-guardBits+i;
            cout << nvarName.str() << endl;
            cVars[i] = ScaLP::newBinaryVariable(nvarName.str());
            cTerm.add(cVars[i],  ( 1ULL << (prodWidth-wOut-guardBits+i)));
            obj.add(cVars[i], (double)0.65);    //append variable to cost function
        }
        ScaLP::Constraint cConstraint = cTerm - Cvar == 0;
        stringstream cName;
        cName << "cCons";
        cConstraint.name = cName.str();
        solver->addConstraint(cConstraint);

        //Limit value of constant to center the error around 0 to the error budget
        ScaLP::Constraint cLimConstraint = Cvar <= errorBudget;
        stringstream cLimName;
        cLimName << "cLimCons";
        cLimConstraint.name = cLimName.str();
        solver->addConstraint(cLimConstraint);

        //Limit the error budget
        cout << "  maxErr=" << errorBudget << endl;
        ScaLP::Constraint truncConstraint = maxEpsTerm - Cvar <= errorBudget;
        stringstream consName;
        consName << "maxEps";
        truncConstraint.name = consName.str();
        solver->addConstraint(truncConstraint);
    }

    // Set the Objective
    cout << "   setting objective (minimize cost function)..." << endl;
    solver->setObjective(ScaLP::minimize(obj));

    // Write Linear Program to file for debug purposes
    cout << "   writing LP-file for debuging..." << endl;
    solver->writeLP("tile.lp");
}


#endif

    void TilingStrategyOptimalILP::computeTruncMultParams(int w, int &g, int &k, long long &errorBudget){
        // first loop iterates over the columns, right to left
        int i = 0;
        long long weightedSumOfTruncatedBits = 0, constant;
        bool loop=true;
        while(loop){
            i++;
            weightedSumOfTruncatedBits += i * (1<<(i-1));
            constant = (1<<(w-1)) - (1<<(i-1));
            errorBudget = (1<<(w-1)) + constant;
            loop = (weightedSumOfTruncatedBits < errorBudget);
        } // when we exit the loop, we have found g
        g = w-(i-1);
        // Now add back bits in rigthtmost column, one by one
        k = 0;
        while(weightedSumOfTruncatedBits >= errorBudget) {
            weightedSumOfTruncatedBits -= (1<<(i-1));
            k++;
        }
    }


}   //end namespace flopoco
