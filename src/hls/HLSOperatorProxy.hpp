#ifndef flopoco_random_hls_HLSOperatorProxy_hpp
#define flopoco_random_hls_HLSOperatorProxy_hpp

#include "HLSOperator.hpp"

namespace flopoco
{

    template<class TOp>
    class HLSOperatorProxy
        : public HLSOperator
    {
    protected:
        const TOp *m_op;
        
        HLSOperatorProxy(const TOp *op)
            : m_op(op)
        {}
    public:       
        virtual const Operator &getOperator() const
        {
            return *m_op;
        }
        
        virtual void releaseHLS() const
        {
            // We are the point of lifetime management
            delete this;
        }
    };
    

}; // flopoco
    
#endif
