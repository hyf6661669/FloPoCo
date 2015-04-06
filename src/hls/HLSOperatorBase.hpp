#ifndef HLSOperatorBase_hpp
#define HLSOperatorBase_hpp

#include "HLSOperator.hpp"

namespace flopoco
{
    class HLSOperatorBase
        : public HLSOperator
    {
    public:       
        virtual const Operator &getOperator() const
        {
            const Operator *op=dynamic_cast<const Operator*>(this);
            if(op==NULL)
                throw std::logic_error("HLSOperatorBase::getOperator - This object does not seem to be an Operator.");
            return *op;
        }
        
        virtual void releaseHLS()
        {
            // NOP : Lifetime is managed on the Operator, not this object
        }
    };
}; // flopoco
    
#endif
