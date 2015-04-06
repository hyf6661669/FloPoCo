#ifndef flopoco_hls_scope_hpp
#define flopoco_hls_scope_hpp

#include "HLSContext.hpp"
#include "HLSOperator.hpp"

#include "HLSExpr.hpp"

#include <map>

namespace flopoco
{
    
    class HLSScope
    {
    private:
        std::map<std::string,HLSNodeDeclPtr> m_decls;
        HLSScope *m_pParent;
    
        static thread_local HLSScope *sm_pCurr;

        HLSScope &operator=(const HLSScope &o) = delete;
        HLSScope(const HLSScope &o) = delete;
    public:
        HLSScope(const Operator *op, HLSContext &ctxt);
        
        static HLSScope &getAmbientScope();

        HLSNodeDeclPtr get(const std::string &x);
        
        HLSNodeDeclPtr declare(const std::string &x, unsigned w);
    };

};

#endif
