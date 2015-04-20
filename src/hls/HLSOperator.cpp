#include "HLSOperator.hpp"
#include "HLSScope.hpp"

#include "Table.hpp"
#include "IntAdder.hpp"

#include <cassert>

namespace flopoco
{
    

    HLSOperator *makeTableHLSOperator(const Table *table);
    HLSOperator *makeIntAdderHLSOperator(const IntAdder *adder);
        
    HLSOperator *getHLSOperator(const Operator *op)
    {
        Operator *dop=const_cast<Operator*>(op);


        HLSOperator *hlsOp=dynamic_cast<HLSOperator*>(dop);
        if(hlsOp)
            return hlsOp;
        
        if(dynamic_cast<const Table*>(op)){
            return makeTableHLSOperator(dynamic_cast<const Table*>(op));
        }
        if(dynamic_cast<const IntAdder*>(op)){
            return makeIntAdderHLSOperator(dynamic_cast<const IntAdder*>(op));
        }
        
        std::stringstream acc;
        acc<<"getHLSOperator - Operator '"<<op->getName()<<"' does not support HLS.";
        throw std::runtime_error(acc.str());
    }
};
