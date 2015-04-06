#ifndef HLSContext_hpp
#define HLSContext_hpp

#include "Operator.hpp"

#include <set>
#include <string>
#include <iostream>
#include <memory>

#include <stdarg.h>

#include "hls/HLSTypes.hpp"

namespace flopoco
{
    class HLSExpr;
    
    class HLSOperator;

    class HLSContext
    {
    private:    
        FILE *m_dst;
        std::string m_prefix;
        bool m_sol;
    
        std::set<std::string> m_functionDecls;
        std::set<std::string> m_functionDefs;

      void writeImpl(const char *msg, va_list args, bool eol);


      // This is wierd, as it's a reference to *this. It models
      // some older behaviour, and looks nicer for internal writes
      HLSContext &dst;
    public:
        HLSContext(
            FILE *dst
        );

        bool isTargetTool(const std::string &name) const;
        
        /*! Returns the HLS specific type for whatever the signal is */
        HLSTypePtr makeType(
            const Signal &sig
        );

        HLSExpr makeIntConstant(
				  bool isSigned,
				  unsigned width,
				  const mpz_class &val
				  );
        
        void emitInputParameter(
            const Operator &op,
            const Signal &sig
        );
        
        void emitOutputParameter(
            const Operator &op,
            const Signal &sig
        );
       
        //! Emit the HLS prototype for the given function (used for either decl or def)
        void emitSignature(
            const Operator &op
        );
        
 
        void emitDeclaration(
            const Operator &op
        );

        void emitDefinition(
            const HLSOperator &op
        );
        

        void indent();
        void unindent();

        void write(const std::string &msg, ...);

        void writeLine(const std::string &msg, ...);

        void writeLine();

        template<class T>
        HLSContext &operator<<(const T &x)
        {
            std::stringstream acc;
            acc<<x;
            write("%s", acc.str().c_str());
            return *this;
        }
        
        //! Get an input/output, or a previously declared variable
        HLSExpr get(const std::string &name);
        
        //! Declare an unsigned integer variable with given width
        HLSExpr declare(const std::string &name, const HLSType &type);

        std::string strRep(const HLSExpr &e);

        std::string strRep(const HLSTypePtr &t);
    };
    

};

#endif
