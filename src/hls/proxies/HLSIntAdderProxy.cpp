
#include "hls/HLSOperatorProxy.hpp"
#include "hls/HLSContext.hpp"
#include "hls/HLSExpr.hpp"
#include "IntAdder.hpp"

namespace flopoco
{
    class IntAdderHLSProxy
      : public HLSOperatorProxy<IntAdder>
    {

    public:
        IntAdderHLSProxy(
            const IntAdder *adder
        )
            : HLSOperatorProxy<IntAdder>(adder)
        {}

      virtual HLSOperator *clone() const override
      {
	return new IntAdderHLSProxy(m_op);
      }
    
        void emitHLSBody(HLSContext &ctxt, HLSScope &scope) const
        {
	  hls_get("R").assign( hls_get("X") + hls_get("Y") );
	    }
    };

    HLSOperator *makeIntAdderHLSOperator(const IntAdder *adder)
    {
        return new IntAdderHLSProxy(adder);
    }
    
}; // flopoco
