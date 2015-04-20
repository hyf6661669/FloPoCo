
#include "hls/HLSOperatorProxy.hpp"
#include "hls/HLSContext.hpp"
#include "Table.hpp"

#include "hls/HLSTypes.hpp"
#include "hls/HLSExpr.hpp"

namespace flopoco
{
    class TableHLSProxy
      : public HLSOperatorProxy<Table>
    {

    public:
        TableHLSProxy(
            const Table *table
        )
            : HLSOperatorProxy<Table>(table)
        {}
    
        void emitHLSBody(HLSContext &ctxt, HLSScope &scope) const
        {
            int minIn=m_op->minIn, maxIn=m_op->maxIn;
            int wOut=m_op->wOut, wIn=m_op->wIn;
            
            if(minIn!=0){
                throw std::runtime_error("Haven't specialsied for minIn!=0 yet.");
            }
            
            HLSTypePtr addrType=HLSTypeInt::create(false, wIn);
            HLSTypePtr eltType=HLSTypeInt::create(false, wOut);
            
            if(ctxt.isTargetTool("c++11")){
                ctxt.writeLine("static const %s romData[%u] = {", ctxt.strRep(eltType).c_str(), maxIn-minIn);
                ctxt.indent();
                for(int i=0; i<=maxIn; i++){
                	ctxt.write("%s", ctxt.strRep(hls_cg(wOut, const_cast<Table*>(m_op)->function(i))).c_str());
                    if(i!=maxIn)
                        ctxt.write(",");
                    ctxt.writeLine();
                }
                ctxt.unindent();
                ctxt.writeLine("};");

                hls_get("Y") = HLSExpr(HLSNodeOpaque::create(eltType, "romData[Y]", "c++11"));
            }else{
                throw std::runtime_error("Can't generate code for this HLS target.");
            }
        }
    };

    HLSOperator *makeTableHLSOperator(const Table *table)
    {
        return new TableHLSProxy(table);
    }
    
}; // flopoco
