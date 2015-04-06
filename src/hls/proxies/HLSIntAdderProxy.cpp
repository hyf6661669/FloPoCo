
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
    
        void emitHLSBody(HLSContext &ctxt, HLSScope &scope) const
        {
            int wIn=m_op->getWIn();
            
            hls_get("R") = hls_get("X") + hls_get("Y");
	    }
    };

    HLSOperator *makeIntAdderHLSOperator(const IntAdder *adder)
    {
        return new IntAdderHLSProxy(adder);
    }
    
}; // flopoco
