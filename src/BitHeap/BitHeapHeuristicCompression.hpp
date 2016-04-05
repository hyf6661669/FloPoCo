#ifndef BITHEAPHEURISTICCOMPRESSION_H
#define BITHEAPHEURISTICCOMPRESSION_H

#ifdef HAVE_SCIP


#include "BitHeap.hpp"
#include "BitHeapILPCompression.hpp"

#include <scip/scip.h>
#include <scip/scipdefplugins.h>

#include <vector>
#include <list>
namespace flopoco
{



    class BitHeapHeuristicCompression
    {
    public:
        struct variableCompressor {
            int type;
            int column;
            int startCompressorWidth;
            int middleCompressorWidth;
            int endCompressorWidth;
        };


        BitHeapHeuristicCompression(BitHeap* bh, std::string mode);
        //BitHeapHeuristicCompression(std::string mode);
        //BitHeapHeuristicCompression();
        //BitHeapHeuristicCompression(int value);
        ~BitHeapHeuristicCompression();
        int generateProblem();
        int writeProblem(std::string filename = "");
        int solve();
        int plotSolution();
        int cleanUp();
        void setLowerBounds(string s);

/*
        //needed to simulate compressors with variable width, eg RCA
        struct variableBasicCompressor {
            vector<int> height;
            vector<int> outputs;
            double areaCost;
        }; */

        vector<list<pair<int,int> > > solution; //index is stage s, list contains pairs of compression elements (1st) and column (2nd), note that there may be several compression elements for the same column
        vector<list<pair<int,int> > > preSolution; //needed for the compressors which the preSolution adds
        vector<list<variableCompressor> > varCompSolution;  //this is the part of the solution which contains Compressors with variable length. index is stage s.



        BitHeap* bh_;
        std::string mode;

        double lowerBound;
        bool useHoles;
        bool generateSolution;
        bool passHeuristicSolution;
        bool passMultipleSolutions;
        bool usePreReduction;
        bool reduceILPStageCount;
        bool useMaxHeuristicStageCount;
        bool optimalAndPreReduction;
        bool differentStages;
        bool paSorting;
        bool buildSingleStages;
        bool useSmootherCompressorPlacing;
        bool dontUsePA;
        double lowerBounds[20];
        bool useVariableCompressors;
        bool useCompleteHeuristic;
        bool getLowerBoundsFromBitHeap;

        unsigned minimumHeightReduction;
        //double minimumEffReduction;
        unsigned minFixedStage;
        unsigned maxFixedStage;
        unsigned minSliceCount;
        unsigned numberOfBuildStages;

        struct bhCompressor {
            BasicCompressor * pointer;
            double maxEfficiency;
            unsigned originalPosition;
        };

    protected:

        BitHeapILPCompression bitHeapILPCompression;

        int noOfStages_;
        int noOfColumnsMax;

        int getMaxStageCount(int maxHeight);

        int printBitHeap();
        int printLowerBounds();
        double computeAreaofSolution();

        int noOfColumnsFor606(int s, int c, bool usePreviousStages, bool useSingleCompressor);
        bool checkColumnFor6Bits(int s, int column, bool usePreviousStages);
        int use606(int s, int start, int times, bool usePreviousStages, bool useSingleCompressor);
        int removeBitsIncludePreviousStages(int s, int c, int bitsToRemove);
        bool smootherCompressorPlacing(unsigned s, double efficiency);

        double computeEfficiency(BasicCompressor* comp);
        double computeCompressionRatioPA(BasicCompressor* comp);
        double compEffBitHeap(unsigned s, unsigned c, unsigned compPos);
        void useCompressor(unsigned s, unsigned column, unsigned newCompPos);
        void algorithm();
        void parandehAfshar();
        int generateAllHeu(bool allCombinations, double efficiency);
        void addHeuristicSolutionToILPSolver();
        void preReduction(unsigned minHeight, double minEff);
        void setUpILPForMoreStages(unsigned stages, bool firstOne);
        void buildVariableCompressorSolution();


        void setUpNewStageCountOfILP();

        BasicCompressor *flipflop;


     private:
        vector<vector<int> > newBits;
        vector<vector<int> > originalNewBits;       //stores the Bits before the heuristic Solution. Needed for passing the heuristic solution to the scip-solver. TODO: find better organisation

        vector<bhCompressor> compressors;

        //vector<variableBasicCompressor> variableBCompressors;            //contains the variableBasicCompressors. These are needed to represent variable Compressors

        //simple test for testing the passing of the heuristic solution
        vector<unsigned> originalCompressorPosition; //because we sorted the compressors for the parandeh-afshar heurisitic we need to remember the original position

    };
}



#endif // HAVE_SCIP

#endif // BITHEAPHEURISTICCOMPRESSION_H
