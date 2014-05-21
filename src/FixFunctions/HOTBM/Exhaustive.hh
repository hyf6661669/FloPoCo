#ifndef _EXHAUSTIVE_HH_
#define _EXHAUSTIVE_HH_

#include <list>

#include "Param.hh"
#include "HOTBMInstance.hh"

using namespace std;

namespace flopoco{


	class Exhaustive {
	public:
		Exhaustive(Function &f_, Param &p_);
		~Exhaustive();

		HOTBMInstance *getInstance();
		static double score(const HOTBMInstance &inst);
	
		enum ScoreType{
			ScoreAreaSquaredDelay, //! score = A^2 * D (default)
			ScoreAreaDelay, //! score = A*D
			ScoreArea	//! score = A
		};
		
		//! Set the method used to score different implementations
		static ScoreType setScoreType(ScoreType score);
		
		//! Get the method used to score different implementations
		static ScoreType getScoreType();
	
		//! Set the minimum alpha to be considered (i.e. always use at least 2^alpha segments)
		/*! By default this is zero, so a single segment is fine. */
		static int setMinAlpha(int alpha);
		
		//! Set the minimum alpha to be considered (i.e. always use at least 2^alpha segments)
		/*! By default this is unbounded, any number of segments is allowed. */
		static int setMaxAlpha(int x);
		
		static int getMinAlpha();
		static int getMaxAlpha();

	private:
		typedef pair<HOTBMInstance *, double> tInstance;
		struct ltInstance {
			bool operator()(const tInstance &i1, const tInstance &i2) const;
		};
		typedef multiset<tInstance, ltInstance> tInstSet;


		int process(list<Param> &pList, tInstSet &instSet, int nMax = -1);
		void skim(tInstSet &instSet, int nMax);

		Function &f;
		Param p;

		HOTBMInstance *instance;
		
		static ScoreType scoreType;
		static int minAlpha;
		static int maxAlpha;
	};
}
#endif // _EXHAUSTIVE_HH_
