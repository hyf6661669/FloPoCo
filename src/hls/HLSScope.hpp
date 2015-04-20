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
        HLSContext &m_ctxt;
        bool m_flushed;
    
        static thread_local HLSScope *sm_pCurr;

        HLSScope &operator=(const HLSScope &o) = delete;
        HLSScope(const HLSScope &o) = delete;

    public:
        HLSScope(const Operator *op, HLSContext &ctxt);
        ~HLSScope();
        
        // Must be called before the object goes out of scope
        void flush();

        static HLSScope &getAmbientScope();

        static HLSContext &getAmbientContext();

      HLSNodeDeclPtr call(const HLSOperator *op, std::string name, const std::map<std::string,HLSExpr> &inputs, const std::map<std::string,std::string> &outputs);

        HLSNodeDeclPtr get(const std::string &x);
        
      HLSNodeDeclPtr declare(const std::string &x, unsigned w);
    };

  HLSExpr hls_declare(const std::string &name, int width);
  HLSExpr hls_get(const std::string &name);

  HLSExpr hls_call(const HLSOperator *op, std::string instName,  HLSExpr arg0, const std::string &resName);
  HLSExpr hls_call(const HLSOperator *op, std::string instName, const std::map<std::string,HLSExpr> &args, const std::string &resName);

  /*! Create a call to a sub-operator
    \param op The operator to instantiate
    \param instName The name of the instance in the code (not used by all targets)
    \param inputs A mapping of the operator argument names to (valid) expressions
    \param outputs A mapping of output output argument names to a set of unique external names

    The call will appear as a HLSNodeCall declaration with name "instName" and a Void type. The
    outputs will appear as HLSNodeCallOutput declarations with type inferred from the operator
  */

  void hls_call(const HLSOperator *op, std::string instName, const std::map<std::string,HLSExpr> &inputs, const std::map<std::string,std::string> &outputs);


};

#endif
