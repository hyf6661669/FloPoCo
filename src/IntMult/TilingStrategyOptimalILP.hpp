#pragma once

#include "TilingStrategy.hpp"

#ifdef HAVE_SCALP
#include <ScaLP/Solver.h>
#include <ScaLP/Exception.h>    // ScaLP::Exception
#include <ScaLP/SolverDynamic.h> // ScaLP::newSolverDynamic
#endif //HAVE_SCALP
#include <iomanip>
#include "BaseMultiplier.hpp"
#include "BaseMultiplierDSPSuperTilesXilinx.hpp"
#include "BaseMultiplierIrregularLUTXilinx.hpp"
#include "MultiplierTileCollection.hpp"

namespace flopoco {

/*!
 * The TilingStrategyOptimalILP class
 */
class TilingStrategyOptimalILP : public TilingStrategy
{

public:
    using TilingStrategy::TilingStrategy;
    TilingStrategyOptimalILP(
        unsigned int wX,
        unsigned int wY,
        unsigned int wOut,
        bool signedIO,
        BaseMultiplierCollection* bmc,
		base_multiplier_id_t prefered_multiplier,
		float occupation_threshold,
		int maxPrefMult,
        MultiplierTileCollection tiles_,
        unsigned guardBits,
        unsigned keepBits,
        unsigned long long errorBudget,
        unsigned long long &centerErrConstant,
        bool performOptimalTruncation);

    void solve() override;

private:
    float occupation_threshold_;
    int dpX, dpY, dpS, dpC, wS, max_pref_mult_;
    unsigned prodWidth, guardBits, keepBits;
    unsigned long long &centerErrConstant;
    unsigned long long errorBudget;
    vector<BaseMultiplierCategory*> tiles;
    bool performOptimalTruncation;
#ifdef HAVE_SCALP
    void constructProblem();

    ScaLP::Solver *solver;
#endif
    void computeTruncMultParams(int w, int &g, int &k, long long &errorBudget);
};

}
