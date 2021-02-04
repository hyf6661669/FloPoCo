#include "PseudoCompressionStrategyOptILP.hpp"
#include "ModReduction/RowAdder.hpp"

using namespace std;
namespace flopoco {

    PseudoCompressionptILP::PseudoCompressionptILP(BitHeap* bitheap):CompressionStrategy(bitheap)
	{

	}

    PseudoCompressionptILP::PseudoCompressionptILP(BitHeap* bitheap, unsigned wIn_, unsigned mod_):CompressionStrategy(bitheap), wIn(wIn_), mod(mod_)
    {

    }

void PseudoCompressionptILP::solve()
{

#ifndef HAVE_SCALP
    throw "Error, TilingAndCompressionOptILP::solve() was called but FloPoCo was not built with ScaLP library";
#else
    cout << "using ILP solver " << bitheap->getOp()->getTarget()->getILPSolver() << endl;
    solver = new ScaLP::Solver(ScaLP::newSolverDynamic({bitheap->getOp()->getTarget()->getILPSolver(),"Gurobi","CPLEX","SCIP","LPSolve"}));
    solver->timeout = bitheap->getOp()->getTarget()->getILPTimeout();

    //for debugging it might be better to order the compressors by efficiency
    orderCompressorsByCompressionEfficiency();

    ScaLP::status stat;
    s_max = 0;
    do{
        s_max++;
        //if(s_max > 4) break;
        solver->reset();
        cout << "constructing model with" << s_max << " stages..." << endl;
        constructProblem(s_max);

        // Try to solve
        cout << "starting solver, this might take a while..." << endl;
        solver->quiet = false;
        stat = solver->solve();

        // print results
        cerr << "The result is " << stat << endl;
        cout << solver->getResult() << endl;
    } while(stat == ScaLP::status::INFEASIBLE);

    ScaLP::Result res = solver->getResult();


#endif
}

void PseudoCompressionptILP::compressionAlgorithm() {
#ifndef HAVE_SCALP
		throw "Error, TilingAndCompressionOptILP::compressionAlgorithm() was called but FloPoCo was not built with ScaLP library";
#else
		//adds the Bits to stages and columns
		orderBitsByColumnAndStage();

		//populates bitAmount
		fillBitAmounts();

		//prints out how the inputBits of the bitheap looks like
		printBitAmounts();

        solve();

		//new solution
		CompressionStrategy::solution = BitHeapSolution();
		CompressionStrategy::solution.setSolutionStatus(BitheapSolutionStatus::OPTIMAL_PARTIAL);

		float compressor_cost = 0;
		vector<vector<int>> zeroInputsVector(s_max, vector<int>((int)wIn, 0));
        vector<vector<vector<int>>> row_adder(s_max, vector<vector<int>>((int)wIn, vector<int>(3, 0)));
		resizeBitAmount(s_max-1);
		ScaLP::Result res = solver->getResult();
		for(auto &p:res.values)
		{
			if(p.second > 0.5){     //parametrize all multipliers at a certain position, for which the solver returned 1 as solution, to flopoco solution structure
				std::string var_name = p.first->getName();
				//cout << var_name << "\t " << p.second << endl;
				//if(var_name.substr(0,1).compare("k") != 0) continue;
				switch(var_name.substr(0,1).at(0)) {
					case 'k':{      //decision variables 'k' for placed compressors
						int sta_id = stoi(var_name.substr(2, dpSt));
						int com_id = stoi(var_name.substr(2 + dpSt + 1, dpK));
						int col_id = stoi(var_name.substr(2 + dpSt + 1 + dpK + 1, dpC));
                        cout << "is compressor" << com_id << " stage " << sta_id << " column " << col_id
                             << " compressor type " << ((com_id<(int)possibleCompressors.size())?possibleCompressors[com_id]->getStringOfIO():"?") << endl;
                        if((int)possibleCompressors.size() <= com_id){
                            row_adder[sta_id][col_id][(com_id-possibleCompressors.size())%3] += 1;
                            //cout << "break" << endl;
                            break;
                        }
						compressor_cost += ((com_id < (int)possibleCompressors.size() && sta_id < s_max - 1)?possibleCompressors[com_id]->area:0) * p.second;
                        if (possibleCompressors[com_id] != flipflop && sta_id < s_max) {
							float instances = p.second;
							while (instances-- > 0.5) {
								CompressionStrategy::solution.addCompressor(sta_id, col_id, possibleCompressors[com_id]);
							}
						} else {
                            cout << "skipped " << possibleCompressors[com_id]  << "==" << flipflop << " " << (possibleCompressors[com_id] == flipflop) << endl;
                        }
						break;
					}
					case 'Z':{          //compressor input bits that are set zero have to be recorded for VHDL-Generation
						int sta_id = stoi(var_name.substr(2, dpSt));
						int col_id = stoi(var_name.substr(2 + dpSt + 1, dpC));
						zeroInputsVector[sta_id][col_id] = (int)(-1);
						break;
					}
					default:            //Other decision variables are not needed for VHDL-Generation
						break;
				}
				for(int stage = 0; stage < s_max; stage++){
					CompressionStrategy::solution.setEmptyInputsByRemainingBits(stage, zeroInputsVector[stage]);
				}
			}
		}

		cout << "Total compressor LUT-cost: " << compressor_cost << endl;

        replace_row_adders(CompressionStrategy::solution, row_adder);

        //reports the area in LUT-equivalents
		printSolutionStatistics();
        drawBitHeap();
		//here the VHDL-Code for the compressors as well as the bits->compressors->bits are being written.
		applyAllCompressorsFromSolution();
#endif
}


#ifdef HAVE_SCALP
void PseudoCompressionptILP::constructProblem(int s_max)
{
    cout << "constructing problem formulation..." << endl;

    //Assemble cost function, declare problem variables
    cout << "   assembling cost function, declaring problem variables..." << endl;
    ScaLP::Term obj;
    //wIn = 8;
    int ns = wS-1; dpS = 1;     //calc number of decimal places, for var names
    int nk = possibleCompressors.size()+4, nc = wIn + 1, nst = s_max; dpK = 1; dpC = 1; dpSt = 1;
    while (ns /= 10)
        dpS++;
    while (nk /= 10)
        dpK++;
    while (nc /= 10)
        dpC++;
    while (nst /= 10)
        dpSt++;

    vector<ScaLP::Term> bitsinColumn(wIn + 4);
    // add the Constraints
    cout << "   adding the constraints to problem formulation..." << endl;


    for(unsigned i = 0; i < wIn; i++){               //Fill array for the bits initially available on Bitheap
        bitsinColumn[i].add(bitAmount[0][i]);
    }

    // BitHeap compression part of ILP formulation
    bitheap->final_add_height = 2;   //Select if dual or ternary adder is used in final stage of bitheap compression
    //remove_all_but_Adders();    //Simplify Modell for test purposes -- REMOVE WHEN PROBLEMS SOLVED ------------------------
    addPseudocompressors(wIn, mod);     //Add Pseudocompressors to list of compressors
    addFlipFlop();              //Add FF to list of compressors
    vector<vector<ScaLP::Variable>> bitsInColAndStage(s_max+1, vector<ScaLP::Variable>(bitsinColumn.size()+1));
    vector<vector<ScaLP::Variable>> range_limits(s_max+1,vector<ScaLP::Variable>(2));
    vector<vector<ScaLP::Variable>> possibleConstBitsPos(s_max,vector<ScaLP::Variable>(wIn+1));
    vector<ScaLP::Variable> constBits(wIn+1);
    vector<ScaLP::Term> neg_range_change(s_max+1), pos_range_change(s_max+1);
    //ScaLP::Term selectLastStage;
    ScaLP::Term sign_ext_vect;
    init_cons_bit_vals(sign_ext_vect, possibleConstBitsPos, constBits);
    unsigned long max_sign_ext_val = 0;                                                                                 //TODO: extend precision
    for(int s = 0; s <= s_max; s++){
        vector<ScaLP::Term> bitsinNextColumn(wIn + 1);
        vector<ScaLP::Term> bitsinCurrentColumn(wIn + 1);
        vector<ScaLP::Term> rcdDependencies(wIn + 1);
        vector<ScaLP::Term> used_pseudocom_in_col(wIn);
        ScaLP::Term rangeMax, rangeMin;

        for(unsigned c = 0; c <= wIn; c++){        //one bit more for carry of ripple carry adder

            if(s < s_max){
                for(unsigned e = 0; e <= possibleCompressors.size() + 2; e++){                                          //place every possible compressor on current position and the following that get covered by it, index extended by 3 for RCA
                    if(e < possibleCompressors.size() - 1 && s < s_max-1){                                                                         //FFs and the row adders, are handled different since they currently can not be instantiated
                        if(possibleCompressors[e]->type.compare("pseudo") != 0 || c == 0){                               //pseudocompressors are only used right aligned on BitHeap
                            stringstream nvarName;
                            nvarName << "k_" << setfill('0') << setw(dpSt) << s << "_" << setfill('0') << setw(dpK) << e << "_" << setfill('0') << setw(dpC) << c;
                            //std::cout << nvarName.str() << " pscom: " << (used_pseudocom_in_col[c][0] != nullptr) << endl;
                            ScaLP::Variable tempV = ScaLP::newIntegerVariable(nvarName.str(), 0, ((possibleCompressors[e]->type.compare("pseudo") != 0)?ScaLP::INF():1));  //pseudocompressors can only be used once per column
                            obj.add(tempV, possibleCompressors[e]->area );                    //append variable to cost function, for r.c.a.-area (cost) is 1 6LUT, FFs cost 0.5LUT
                            for(int ce = 0; ce < (int) possibleCompressors[e]->getHeights() && ce < (int)bitsinCurrentColumn.size() - (int)c; ce++){                                //Bits that can be removed by compressor e in stage s in column c for constraint C1
                                //cout << possibleCompressors[e]->getHeightsAtColumn((unsigned) ce, false) << " c: " << c+ce << endl;
                                bitsinCurrentColumn[c+ce].add(tempV, possibleCompressors[e]->getHeightsAtColumn((unsigned) ce, false));
                                if(c+ce < wIn && possibleCompressors[e]->type.compare("pseudo") == 0){
                                    used_pseudocom_in_col[c+ce].add(tempV, -(int)possibleCompressors[e]->getHeightsAtColumn((unsigned) ce, false));
                                }
                                if(possibleCompressors[e]->getHeightsAtColumn((unsigned) ce, false)) {
                                    if(possibleCompressors[e]->type.compare("pseudo") == 0){
                                        pos_range_change[s].add(tempV, -(1 << ce));                                 //consider the range change of the bit removed by the pseudo-compressor
                                    }
                                }
                            }
                            for(int ce = 0; ce < (int) possibleCompressors[e]->getOutHeights() && ce < (int)bitsinNextColumn.size() - (int)c; ce++){   //Bits that are added by compressor e in stage s+1 in column c for constraint C2
                                //cout << possibleCompressors[e]->getOutHeightsAtColumn((unsigned) ce, false) << " c: " << c+ce << endl;
                                if(ce <= possibleCompressors[e]->ones_vector_start){
                                    bitsinNextColumn[c+ce].add(tempV, possibleCompressors[e]->getOutHeightsAtColumn((unsigned) ce, false));
                                } else {
                                    cout << "compr: " << e << " height: " << possibleCompressors[e]->getOutHeightsAtColumn((unsigned) ce, false) << " c: " << c+ce << endl;
                                }
                            }
                            if(possibleCompressors[e]->type.compare("pseudo") == 0){
                                //cout << "stage " << s << " col " << c << " comp " << e << " range chg " << possibleCompressors[e]->range_change << endl;
                                if(possibleCompressors[e]->range_change < 0){
                                    neg_range_change[s].add(tempV,possibleCompressors[e]->range_change);
                                    for(unsigned ci = possibleCompressors[e]->ones_vector_start; ci <= wIn; ci++){      //contributions to sign externsion vector of current neg. congruent pseudo-compressor
                                        sign_ext_vect.add(tempV, 1<<ci );
                                        max_sign_ext_val += (1<<ci);
                                    }
                                    cout << "pseudocompr. " << possibleCompressors[e]->range_change << " start of ones vector " << possibleCompressors[e]->ones_vector_start << endl;
                                } else {
                                    pos_range_change[s].add(tempV,possibleCompressors[e]->range_change);
                                }
                            }
                        }
                    } else {                                                                                            //Handling for FFs and row adders
                        if(possibleCompressors.size() <= e && e <= possibleCompressors.size()+2 && (c!=0 || possibleCompressors.size() == e) ) {             //MSB and middle element of row adder should not be put in the LSB column of the BitHeap
                            cout << "-------------------------- stage" << s << " col " << c << " compr" << e << " ---------------------------------- " << endl;                            stringstream nvarName;
                            nvarName << "k_" << setfill('0') << setw(dpSt) << s << "_" << setfill('0') << setw(dpK) << e << "_" << setfill('0') << setw(dpC) << c;
                            //cout << nvarName.str() << endl;
                            ScaLP::Variable tempV = ScaLP::newBinaryVariable(nvarName.str());
                            obj.add(tempV, (e == possibleCompressors.size()-1)?0.5:1);                                        //append variable to cost function, r.c.a.-area (cost) is 1 6LUT, FFs cost 0.5LUT
                            int takes_carry = (e == possibleCompressors.size() && bitheap->final_add_height == 2)?1:0;              //The dual input adder can process an additional carry input bit in its LSB column
                            bitsinCurrentColumn[c].add(tempV, ((e == possibleCompressors.size()-1)?1: (e == possibleCompressors.size()+2) ? bitheap->final_add_height : (bitheap->final_add_height + takes_carry) ));   //FFs remove only one bit from current stage, but the ripple carry adder two, the front element of ripple carry adder does not remove any bit but just provides the carry
                            bitsinNextColumn[c + 0].add(tempV, 1);
                            if(e == possibleCompressors.size()+2 && c < wIn) bitsinNextColumn[c + 1].add(tempV, 1);
                            if(0 < c && possibleCompressors.size()+1 <= e)                                                          //Middle and left element of RCA are counted negative in eq. for relations between RCA elements (C5)
                                rcdDependencies[c-1].add(tempV, -1);
                            if(c < bitsinColumn.size()-1 && possibleCompressors.size() <= e && e < possibleCompressors.size()+2)    //Middle and right element of RCA are counted positive in eq. for relations between RCA elements (C5)
                                rcdDependencies[c].add(tempV, 1);
                            if(c == bitsinColumn.size() && e == possibleCompressors.size()+2)                                       //only the left (and not the middle) element of the RCA should be put in the MSB column of the bitheap
                                rcdDependencies[c].add(tempV, 1);
                        }
                    }
                }
            }

            if(bitsInColAndStage[s][c] == nullptr){                                                                 //N_s_c: Bits that enter current compressor stage
                stringstream curBits;
                curBits << "N_" << s << "_" << c;
                //cout << curBits.str() << endl;
                bitsInColAndStage[s][c] = ScaLP::newIntegerVariable(curBits.str(), 0, ScaLP::INF());
            }
            if(s == 0 && bitsInColAndStage[s][c] != nullptr){
                C0_bithesp_input_bits(s, c, bitsinColumn, bitsInColAndStage);
            }
            if(s < s_max){
                bitsinNextColumn[c].add(possibleConstBitsPos[s][c], 1);                                     //possible bit form sign extension constant vector
                C1_compressor_input_bits(s, c, bitsinCurrentColumn, bitsInColAndStage);
                C2_compressor_output_bits(s, c, bitsinNextColumn, bitsInColAndStage);
            }

            if(s_max-1 <= s){                                               //limitation of the BitHeap height to final_add_height in s_max-1 and 1 in stage s_max
                C3_limit_BH_height(s, c, bitsInColAndStage);
            }
            if(0 < c && s < s_max){                                         //define constraint for the stucture of the RCA
                C5_RCA_dependencies(s, c, rcdDependencies);
            }
            if(s == 0){                                                     //append the weighted bits on bitheap to the equation describing the range
                rangeMax.add(bitsInColAndStage[s][c], 1<<c);
            } else {
                rangeMax.add(bitsInColAndStage[s-1][c], 1<<c);
            }
            if(c == wIn){                                                   //range constraint C6/C7
                C67_range_constraint(s, range_limits, rangeMin, rangeMax, neg_range_change, pos_range_change);
            }
        }
    }
    C8_sign_extension(sign_ext_vect, max_sign_ext_val);
    C9_sign_extension_bits(possibleConstBitsPos, constBits);

    // Set the Objective
    cout << "   setting objective (minimize cost function)..." << endl;
    solver->setObjective(ScaLP::minimize(obj));

    // Write Linear Program to file for debug purposes
    cout << "   writing LP-file for debuging..." << endl;
    solver->writeLP("pseudocompression.lp");
}


    void PseudoCompressionptILP::C0_bithesp_input_bits(int s, int c, vector<ScaLP::Term> &bitsinColumn, vector<vector<ScaLP::Variable>> &bitsInColAndStage){
        stringstream consName0;
        consName0 << "C0_" << s << "_" << c;
        bitsinColumn[c].add(bitsInColAndStage[s][c], -1);      //Output bits from sub-multipliers
        ScaLP::Constraint c0Constraint = bitsinColumn[c] == 0;     //C0_s_c
        c0Constraint.name = consName0.str();
        solver->addConstraint(c0Constraint);
    }

    void PseudoCompressionptILP::C1_compressor_input_bits(int s, int c, vector<ScaLP::Term> &bitsinCurrentColumn, vector<vector<ScaLP::Variable>> &bitsInColAndStage){
        stringstream consName1, zeroBits;
        consName1 << "C1_" << s << "_" << c;
        zeroBits << "Z_" << setfill('0') << setw(dpSt) << s << "_" << setfill('0') << setw(dpC) << c;
        //cout << zeroBits.str() << endl;
        bitsinCurrentColumn[c].add(ScaLP::newIntegerVariable(zeroBits.str(), 0, ScaLP::INF()), -1);      //Unused compressor input bits, that will be set zero
        bitsinCurrentColumn[c].add(bitsInColAndStage[s][c], -1);      //Bits arriving in current stage of the compressor tree
        ScaLP::Constraint c1Constraint = bitsinCurrentColumn[c] == 0;     //C1_s_c
        c1Constraint.name = consName1.str();
        solver->addConstraint(c1Constraint);
    }

    void PseudoCompressionptILP::C2_compressor_output_bits(int s, int c, vector<ScaLP::Term> &bitsinNextColumn, vector<vector<ScaLP::Variable>> &bitsInColAndStage){
        stringstream consName2, nextBits;
        consName2 << "C2_" << s << "_" << c;
        //cout << consName2.str() << endl;
        if(bitsInColAndStage[s+1][c] == nullptr){
            nextBits << "N_" << s+1 << "_" << c;
            //cout << nextBits.str() << endl;
            bitsInColAndStage[s+1][c] = ScaLP::newIntegerVariable(nextBits.str(), 0, ScaLP::INF());
        }
        bitsinNextColumn[c].add(bitsInColAndStage[s+1][c], -1); //Output Bits of compressors to next stage
        ScaLP::Constraint c2Constraint = bitsinNextColumn[c] == 0;     //C2_s_c
        c2Constraint.name = consName2.str();
        solver->addConstraint(c2Constraint);
    }

    void PseudoCompressionptILP::C3_limit_BH_height(int s, int c, vector<vector<ScaLP::Variable>> &bitsInColAndStage){
        stringstream consName3;
        consName3 << "C3_" << s << "_" << c;
        ScaLP::Constraint c3Constraint = bitsInColAndStage[s][c] <= ((s == s_max-1)?bitheap->final_add_height:1);     //C3_s_c
        c3Constraint.name = consName3.str();
        solver->addConstraint(c3Constraint);
        //cout << consName3.str() << " " << stage.str() << endl;
    }

    void PseudoCompressionptILP::C5_RCA_dependencies(int s, int c, vector<ScaLP::Term> &rcdDependencies){
        stringstream consName5;
        consName5 << "C5_" << s << "_" << c;
        //cout << consName5.str() << rcdDependencies[c-1] << endl;
        ScaLP::Constraint c5Constraint = rcdDependencies[c-1] == 0;     //C5_s_c
        c5Constraint.name = consName5.str();
        solver->addConstraint(c5Constraint);
    }

    void PseudoCompressionptILP::C67_range_constraint(int s, vector<vector<ScaLP::Variable>> &range_limits, ScaLP::Term &rangeMin, ScaLP::Term &rangeMax, vector<ScaLP::Term> &neg_range_change, vector<ScaLP::Term> &pos_range_change){
        for(int con = 0; con <= 1; con++){
            stringstream consName;
            consName << ((!con)?"C6_":"C7_") << s;
            ScaLP::Constraint cxConstraint;
            stringstream nvarName;
            nvarName << "s" << setfill('0') << setw(dpSt) << s << ((con)?"_min":"_max");
            range_limits[s][con] = ScaLP::newIntegerVariable(nvarName.str(), ((con)?-ScaLP::INF():0), ((!con)?ScaLP::INF():0) );
            if(!con && s != s_max) rangeMax.add(range_limits[s][con], -1);
            if( con && s != s_max) rangeMin.add(range_limits[s][con], -1);
            if(!con){
                if(s < s_max) cxConstraint = rangeMax + ((s != 0)?pos_range_change[s-1]:0) == 0;     //C6_s
                if(s ==s_max) cxConstraint = rangeMax + ((s != 0)?pos_range_change[s-1]:0)  < mod;     //C6_s
            } else {
                if(s < s_max) cxConstraint = rangeMin + ((s != 0)?neg_range_change[s-1]:0) == 0;     //C6_s
                if(s ==s_max) cxConstraint = rangeMin + ((s != 0)?neg_range_change[s-1]:0)  >= -(int)mod;     //C6_s
            }
            cxConstraint.name = consName.str();
            solver->addConstraint(cxConstraint);
        }
    }

    void PseudoCompressionptILP::init_cons_bit_vals(ScaLP::Term &sign_ext_vect, vector<vector<ScaLP::Variable>> &possibleConstBitsPos, vector<ScaLP::Variable> &constBits){
        for(unsigned ci = 0; ci <= wIn; ci++){      //generate decision variables for bits of the sign extension vector
            stringstream nvarName;
            nvarName << "b_" << setfill('0') << setw(dpC) << ci;
            cout << endl << " " << nvarName.str();
            constBits[ci] = ScaLP::newIntegerVariable(nvarName.str(), 0, 1);
            sign_ext_vect.add(constBits[ci], -(1<<ci) );
            for(unsigned si = 0; si < s_max; si++){
                stringstream bitName;
                bitName << "p_" << setfill('0') << setw(dpSt) << si << "_" << setfill('0') << setw(dpC) << ci;
                cout << " " << bitName.str();
                possibleConstBitsPos[si][ci] = ScaLP::newIntegerVariable(bitName.str(), 0, 1);
            }
        }
    }

    void PseudoCompressionptILP::C8_sign_extension(ScaLP::Term &sign_ext_vect, unsigned long max_sign_ext_val){
        stringstream consName8;
        consName8 << "C8";
        int wMax, wM = 1, wBitsMax = 1;
        while (max_sign_ext_val /= 2)
            wBitsMax++;
        wMax = wBitsMax;
        while (wMax /= 10)
            wM++;
        for(unsigned ci = wIn+1; ci <= wBitsMax; ci++){      //generate decision variables for bits of the sign extension vector
            stringstream nvarName;
            nvarName << "b_" << setfill('0') << setw(wM) << ci;
            ScaLP::Variable tempV = ScaLP::newIntegerVariable(nvarName.str(), 0, 1);
            sign_ext_vect.add(tempV, -(1<<ci) );
        }
        cout << consName8.str() << sign_ext_vect << endl;
        ScaLP::Constraint c8Constraint = sign_ext_vect == 0;     //C8_s
        c8Constraint.name = consName8.str();
        solver->addConstraint(c8Constraint);
    }

    void PseudoCompressionptILP::C9_sign_extension_bits(vector<vector<ScaLP::Variable>> &possibleConstBitsPos, vector<ScaLP::Variable> &constBits){
        for(unsigned ci = 0; ci <= wIn; ci++){
            ScaLP::Term bitPlacement;
            stringstream consName9;
            consName9 << "C9_" << setfill('0') << setw(dpC) << ci;
            for(unsigned si = 0; si < s_max; si++){
                bitPlacement.add(possibleConstBitsPos[si][ci], 1);
            }
            ScaLP::Constraint c9Constraint = bitPlacement - constBits[ci] == 0;     //C8_s
            c9Constraint.name = consName9.str();
            solver->addConstraint(c9Constraint);
        }

    }

    bool PseudoCompressionptILP::addFlipFlop(){
        //BasicCompressor* flipflop;
        bool foundFlipflop = false;
        for(unsigned int e = 0; e < possibleCompressors.size(); e++){
            //inputs
            if(possibleCompressors[e]->getHeights() == 1 && possibleCompressors[e]->getHeightsAtColumn(0) == 1){
                if(possibleCompressors[e]->getOutHeights() == 1 && possibleCompressors[e]->getOutHeightsAtColumn(0) == 1){
                    foundFlipflop = true;
                    flipflop = possibleCompressors[e];
                    break;
                }
            }
        }

        if(!foundFlipflop){
            //add flipflop at back of possibleCompressor
            //vector<int> newVect;
            //BasicCompressor *newCompressor;
            //int col0=1;
            //newVect.push_back(col0);
            BasicCompressor *newCompressor = new BasicCompressor(bitheap->getOp(), bitheap->getOp()->getTarget(), vector<int> {1}, 0.5, "combinatorial", true);
            possibleCompressors.push_back(newCompressor);

            flipflop = newCompressor;
        }

        return !foundFlipflop;
    }

    void PseudoCompressionptILP::remove_all_but_Adders(){
        for(unsigned int e = 0; e < possibleCompressors.size(); e++){
            cout << possibleCompressors[e]->getStringOfIO() << endl;
            if(possibleCompressors[e]->getHeights() != 1 || !(possibleCompressors[e]->getHeightsAtColumn(0) == 2 || possibleCompressors[e]->getHeightsAtColumn(0) == 3)){
                possibleCompressors.erase(possibleCompressors.begin()+e);
                e = 0;
                cout << "removed" << endl;
            }
        }
    }

    void PseudoCompressionptILP::resizeBitAmount(unsigned int stages){

        stages++;	//we need also one stage for the outputbits

        unsigned int columns = bitAmount[bitAmount.size() - 1].size();
        //we need one stage more for the
        while(bitAmount.size() < stages){
            bitAmount.resize(bitAmount.size() + 1);
            bitAmount[bitAmount.size() - 1].resize(columns, 0);
        }
    }

    void PseudoCompressionptILP::printIOHeights(void){
        for(int j=0; j < (int)possibleCompressors.size(); j++){
            cout << "Compressor input heights: ";
            for(int i = possibleCompressors[j]->getHeights()-1; 0 < i; i--){
                cout << possibleCompressors[j]->getHeightsAtColumn(i);
            }
            cout << endl << "Compressor output heights: ";
            for(int o = possibleCompressors[j]->getOutHeights()-1; 0 < o; o--){
                cout << possibleCompressors[j]->getOutHeightsAtColumn(o);
            }
            cout << endl;
        }

    }

    void PseudoCompressionptILP::drawBitHeap(){
        vector<vector<int>> bitsOnBitHeap(s_max+1, vector<int>((int)wIn, 0));
        ScaLP::Result res = solver->getResult();
        for(auto &p:res.values) {
            if (p.second >
                0.5) {     //parametrize all multipliers at a certain position, for which the solver returned 1 as solution, to flopoco solution structure
                std::string var_name = p.first->getName();
                //if(var_name.substr(0,1).compare("k") != 0) continue;
                switch (var_name.substr(0, 1).at(0)) {
 /*                   case 'k': {      //decision variables 'k' for placed compressors
                        int sta_id = stoi(var_name.substr(2, dpSt));
                        int com_id = stoi(var_name.substr(2 + dpSt + 1, dpK));
                        int col_id = stoi(var_name.substr(2 + dpSt + 1 + dpK + 1, dpC));
                        float instances = p.second;
                        while (instances-- > 0.5) {
                            if(com_id < possibleCompressors.size()){
                                for(int ci = 0; ci < possibleCompressors[com_id]->outHeights.size() && ci <= possibleCompressors[com_id]->ones_vector_start; ci++){
                                    bitsOnBitHeap[sta_id][col_id+ci] += possibleCompressors[com_id]->outHeights[ci];
                                }
                            } else {
                                bitsOnBitHeap[sta_id][col_id] += 1;
                                if(com_id == 26)
                                    bitsOnBitHeap[sta_id][col_id+1] += 1;
                            }

                        }
                        cout << var_name << "\t " << p.second << endl;
                        break;
                    }*/
                    case 'p': {          //constant bits from sign extension of negative congruent pseudo-compressors
                        int sta_id = stoi(var_name.substr(2, dpSt));
                        int col_id = stoi(var_name.substr(2 + dpSt + 1, dpC));
                        //bitsOnBitHeap[sta_id][col_id] += 1;
                        cout << var_name << "\t " << p.second << endl;
                        break;
                    }
                    case 'N': {          //bits present in a particular column an stage of BitHeap
                        int sta_id = stoi(var_name.substr(2, dpSt));
                        int col_id = stoi(var_name.substr(2 + dpSt + 1, dpC));
                        if(sta_id < s_max+1 && col_id < wIn)
                            bitsOnBitHeap[sta_id][col_id] = p.second;
                        //cout << var_name << "\t " << p.second << endl;
                        break;
                    }
                    default:            //Other decision variables are not needed for VHDL-Generation
                        break;
                }
            }
        }
        for(int s = 0; s < bitsOnBitHeap.size(); s++){
            for(int c = bitsOnBitHeap[0].size()-1; 0 <= c; c--){
                cout << bitsOnBitHeap[s][c];
            }
            cout << endl;
        }
    }

    void PseudoCompressionptILP::replace_row_adders(BitHeapSolution &solution, vector<vector<vector<int>>> &row_adders){
        cout << solution.getSolutionStatus() << endl;
        for(int s = 0; s < row_adders.size(); s++){
            for(int c = 0; c < row_adders[0].size(); c++){
                //cout << "at stage " << s << " col " << c << " " << row_adders[s][c][0] << row_adders[s][c][1] << row_adders[s][c][2] << endl;
                int ci; bool adder_started;
                while(0 < row_adders[s][c][0]){
                    ci = 0;
                    adder_started = false;
                    while(0 < row_adders[s][c+ci][0] || 0 < row_adders[s][c+ci][1] || 0 < row_adders[s][c+ci][2]){
                        if(0 < row_adders[s][c+ci][0] && adder_started == false){
                            row_adders[s][c+ci][0]--;
                            adder_started = true;
                        } else {
                            if(adder_started || 0 < row_adders[s][c+ci][1] || 0 < row_adders[s][c+ci][2]){
                                if(0 < row_adders[s][c+ci][1]){
                                    row_adders[s][c+ci][1]--;
                                } else {
                                    row_adders[s][c+ci][2]--;
                                    adder_started = false;
                                    cout << "adder in stage " << s <<  "from col " << c << " to " << c+ci << " width " << ci+1 << endl;
                                    BasicCompressor *newCompressor = new BasicRowAdder(bitheap->getOp(), bitheap->getOp()->getTarget(), ci+1);
                                    //possibleCompressors.push_back(newCompressor);
                                    cout << solution.getCompressorsAtPosition(s, c).size() << endl;
                                    solution.addCompressor(s, c, newCompressor);
                                    cout << solution.getCompressorsAtPosition(s, c).size() << " " << solution.getCompressorsAtPosition(s, c)[0].first->outHeights.size() << endl;
                                    cout << "ok" << endl;
                                }
                            }
                        }
                        ci++;
                    }

                }
            }
        }
    }

#endif


}   //end namespace flopoco
