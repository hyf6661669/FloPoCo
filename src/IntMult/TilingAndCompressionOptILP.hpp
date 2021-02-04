#pragma once

#include "TilingStrategy.hpp"
#include "./../BitHeap/CompressionStrategy.hpp"

#ifdef HAVE_SCALP
#include <ScaLP/Solver.h>
#include <ScaLP/Exception.h>    // ScaLP::Exception
#include <ScaLP/SolverDynamic.h> // ScaLP::newSolverDynamic
#endif //HAVE_SCALP
#include <iomanip>
#include "BaseMultiplier.hpp"
#include "BaseMultiplierDSPSuperTilesXilinx.hpp"
#include "MultiplierTileCollection.hpp"

namespace flopoco {

/*!
 * The TilingAndCompressionOptILP class
 */
class TilingAndCompressionOptILP : public TilingStrategy, public CompressionStrategy
{

public:
    using TilingStrategy::TilingStrategy;
    TilingAndCompressionOptILP(
        unsigned int wX,
        unsigned int wY,
        unsigned int wOut,
        bool signedIO,
        BaseMultiplierCollection* bmc,
		base_multiplier_id_t prefered_multiplier,
		float occupation_threshold,
		int maxPrefMult,
        MultiplierTileCollection tiles_,
        BitHeap* bitheap);

    void solve() override;
	void compressionAlgorithm() override;

private:
    float occupation_threshold_;
    int dpX, dpY, dpS, wS, dpK, dpC, dpSt, s_max, max_pref_mult_;
    unsigned prodWidth;
	vector<BaseMultiplierCategory*> tiles;

#ifdef HAVE_SCALP
    BasicCompressor* flipflop;
    void constructProblem(int s_max);
    bool addFlipFlop();

    ScaLP::Solver *solver;

    void resizeBitAmount(unsigned int stages);
#endif
    };

}
