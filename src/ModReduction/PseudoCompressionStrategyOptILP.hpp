#pragma once

#include "./../BitHeap/CompressionStrategy.hpp"

#ifdef HAVE_SCALP
#include <ScaLP/Solver.h>
#include <ScaLP/Exception.h>    // ScaLP::Exception
#include <ScaLP/SolverDynamic.h> // ScaLP::newSolverDynamic
#endif //HAVE_SCALP
#include <iomanip>


namespace flopoco {

/*!
 * The TilingAndCompressionOptILP class
 */
class PseudoCompressionptILP : public CompressionStrategy
{

public:
    PseudoCompressionptILP(BitHeap* bitheap);
    PseudoCompressionptILP(BitHeap* bitheap, unsigned wIn, unsigned mod);

    void solve();
	void compressionAlgorithm();

private:

    int dpS, wS, dpK, dpC, dpSt, s_max;
    unsigned wIn, mod;

#ifdef HAVE_SCALP
    BasicCompressor* flipflop;
    vector<BasicCompressor*> pseudocompressors;
    void constructProblem(int s_max);
    bool addFlipFlop();
    void remove_all_but_Adders();
    void addPseudocompressors();

    ScaLP::Solver *solver;

    void resizeBitAmount(unsigned int stages);
    void printIOHeights(void);

    void C0_bithesp_input_bits(int s, int c, vector<ScaLP::Term> &bitsinColumn,
                               vector<vector<ScaLP::Variable>> &bitsInColAndStage);
    void C1_compressor_input_bits(int s, int c, vector<ScaLP::Term> &bitsinCurrentColumn,
                                  vector<vector<ScaLP::Variable>> &bitsInColAndStage);
    void C2_compressor_output_bits(int s, int c, vector<ScaLP::Term> &bitsinNextColumn,
                                   vector<vector<ScaLP::Variable>> &bitsInColAndStage);
    void C3_limit_BH_height(int s, int c, vector<vector<ScaLP::Variable>> &bitsInColAndStage);
    void C5_RCA_dependencies(int s, int c, vector<ScaLP::Term> &rcdDependencies);
    void C67_range_constraint(int s, vector<vector<ScaLP::Variable>> &range_limits, ScaLP::Term &rangeMin,
                              ScaLP::Term &rangeMax, vector<ScaLP::Term> &neg_range_change,
                              vector<ScaLP::Term> &pos_range_change);
    void init_cons_bit_vals(ScaLP::Term &sign_ext_vect, vector<vector<ScaLP::Variable>> &possibleConstBitsPos, vector<ScaLP::Variable> &constBits);
    void C8_sign_extension(ScaLP::Term &sign_ext_vect, unsigned long max_sign_ext_val);
    void C9_sign_extension_bits(vector<vector<ScaLP::Variable>> &possibleConstBitsPos, vector<ScaLP::Variable> &constBits);
#endif



    };

}
