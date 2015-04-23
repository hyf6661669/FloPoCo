#include "HLSContext.hpp"
#include "HLSOperator.hpp"
#include "HLSScope.hpp"

#include "HLSExpr.hpp"

#include <cassert>

namespace flopoco
{
    
    HLSContext::HLSContext(
        FILE *dst
    ) : m_dst(dst)
      , m_prefix("")
      , m_sol(true)
      , dst(*this)
    {
    }

    bool HLSContext::isTargetTool(const std::string &msg) const
    {
    return msg=="c++11";
    }
    
    void HLSContext::indent()
    {
      if(!m_sol)
        writeLine();
      m_prefix=m_prefix+"    ";
    }
    
    void HLSContext::unindent()
    {
      if(!m_sol)
        writeLine();
      m_prefix=m_prefix.substr(0, m_prefix.size()-4);
    }

  void HLSContext::writeImpl(const char *msg, va_list args, bool eol)
  {
    if(m_sol){
      fputs(m_prefix.c_str(), m_dst);
    }
    m_sol=false;
    vfprintf(m_dst, msg, args);
    if(eol){
      fputs("\n", m_dst);
      m_sol=true;
    }
  }

  void HLSContext::write(const std::string &msg, ...)
  {
    va_list va;
    va_start(va, msg);
    writeImpl(msg.c_str(), va, false);
    va_end(va);
  }

  void HLSContext::writeLine(const std::string &msg, ...)
  {
    va_list va;
    va_start(va, msg);
    writeImpl(msg.c_str(), va, true);
    va_end(va);
  }

  void HLSContext::writeLine()
  {
    if(m_sol){
      fputs(m_prefix.c_str(), m_dst);
    }
    fputs("\n", m_dst);
    m_sol=true;
  }
    
    HLSTypePtr HLSContext::makeType(
        const Signal &sig
    ){
      return makeHLSType(sig);
    };

  void HLSContext::emitCall(const HLSNodeCall &x)
  {
	    writeLine("// instance %s\n", x.getName().c_str());
	    const HLSOperator *hOp=x.getOperator();
	    const Operator &op=hOp->getOperator();
	    writeLine("// Declaration of outputs");
	    indent();
	    for(int i=0; i<op.getNumberOfOutputs(); i++){
	      std::string name=op.getOutputSignal(i)->getName();
	      HLSNodePtr node=x.getOutput(name);
	      writeLine("%s %s;", strRep(node->getType()).c_str(), name.c_str());
	    }
	    unindent();
	    writeLine("// Actual call");
	    writeLine("%s(", op.getName().c_str());
	    indent();
	    writeLine("// Input bindings:");
	    indent();
	    for(int i=0; i<op.getNumberOfInputs(); i++){
	      std::string name=op.getInputSignal(i)->getName();
	      HLSNodePtr val=x.getInput(name);
	      writeLine("// Input %u, argName=%s", i, name.c_str());
	      write("%s,", strRep(val).c_str());
	    }
	    unindent();
	    writeLine("// Output bindings:");
	    indent();
	    for(int i=0; i<op.getNumberOfOutputs(); i++){
	      std::string name=op.getOutputSignal(i)->getName();
	      std::string decl=x.getOutput(name)->getName();
	      writeLine("// Output %u, argName=%s", i, name.c_str());
	      write("& %s", decl.c_str());
	      if(i+1<op.getNumberOfOutputs()){
		writeLine(",");
	      }else{
		writeLine();
	      }
	    }
	    unindent();
	    unindent();
	    writeLine(");");
  }

    void HLSContext::emitDeclareAndAssign(const HLSNodeVar &v)
    {
    	dst<<strRep(v.getType())<<" "<<v.getName()<<" = "<<strRep(v.getSrc());
    	writeLine();
    }

    void HLSContext::emitAssignOutput(const HLSNodeOutput &o)
    {
    	dst<<o.getName()<<" = "<<strRep(o.getSrc());
    	writeLine();
    }

    std::string HLSContext::strRep(const HLSTypePtr &t)
    {
    	std::stringstream acc;

    	auto pI=std::dynamic_pointer_cast<HLSTypeInt>(t);
    	if(pI){
    		acc<<(pI->isSigned()?"ap_int":"ap_uint")<<"<"<<pI->getWidth()<<">";
    		return acc.str();
    	}

    	auto pB=std::dynamic_pointer_cast<HLSTypeBool>(t);
    	if(pB){
    		return "bool";
    	}

    	auto pFP=std::dynamic_pointer_cast<HLSTypeFloat>(t);
    	if(pFP){
    		acc<<"ap_uint<"<<pFP->getWidth()<<">";
    		return acc.str();
    	}

    	throw std::runtime_error("HLSContext::strRep - Unsupported type.");
    }

    class OutputExpr
		: public HLSNodeVisitor
	{
	private:
    	HLSContext &m_ctxt;
		std::stringstream dst;
	public:
		OutputExpr(HLSContext &ctxt)
			: m_ctxt(ctxt)
		{}

		std::string str()
		{ return dst.str(); }

		bool visitGeneric(const HLSNode &c)
		{
			throw std::runtime_error(std::string("HLSContext/OutputExpr - Haven't handled node type")+typeid(&c).name());
		}

		void visit(const HLSNodeConstantInt &x)
		 {
			if(x.getValue()<0)
				throw std::runtime_error("OutputExpr/HLSNodeConstantInt - Haven't thought about negative consts for HLS yet.");

			dst<<m_ctxt.strRep(x.getType());
			if(x.getValue()>0x7FFFFFFFul){
			  dst<<"(\"";
			  dst<<x.getValue().get_str(16);
			  dst<<"\",16)";
			}else{
			  dst<<"(0x"<<x.getValue().get_str(16)<<")";
			}
		 }

		 void visit(const HLSNodeOutput &x)
		 {
			 // Is this ok?
			throw std::runtime_error("Shouldn't be possible to use an expr in an output (?).");
		 }

		 void visit(const HLSNodeVar &x)
		 {
			 dst<<x.getName();
		 }

		 void visit(const HLSNodeInput &x)
		 {
			 dst<<x.getName();
		 }

		 void visit(const HLSNodeAdd &x)
		 {
			 dst<<"(";
			 x.getLeft()->accept(*this);
			 dst<<"+";
			 x.getRight()->accept(*this);
			 dst<<")";
		 }

		 void visit(const HLSNodeEquals &x)
		 {
			 dst<<"(";
			 x.getLeft()->accept(*this);
			 dst<<"==";
			 x.getRight()->accept(*this);
			 dst<<")";
		 }

		 void visit(const HLSNodeNotEquals &x)
		 {
			 dst<<"(";
			 x.getLeft()->accept(*this);
			 dst<<"!=";
			 x.getRight()->accept(*this);
			 dst<<")";
		 }

		 void visit(const HLSNodeLogicalOr &x)
		 {
			 dst<<"(";
			 x.getLeft()->accept(*this);
			 dst<<"||";
			 x.getRight()->accept(*this);
			 dst<<")";
		 }

		 void visit(const HLSNodeSub &x)
		 {
			 dst<<"(";
			 x.getLeft()->accept(*this);
			 dst<<"-";
			 x.getRight()->accept(*this);
			 dst<<")";
		 }

		 void visit(const HLSNodeReinterpretBits &x)
		 {
			 const HLSTypePtr srcT=x.getSrc()->getType();
			 const HLSTypePtr dstT=x.getType();

			 if(
					 (srcT->isInt()|| srcT->isBool() || srcT->isFloat())
					 &&
					 (dstT->isInt()|| dstT->isBool() || dstT->isFloat())
			 ){
				 x.getSrc()->accept(*this);
			 }else{
				 // Hrmm. Need to think more about this.
				 assert(0); // Not tested
				 dst<<"(reinterpret_cast<";
				 dst<<m_ctxt.strRep(dstT);
				 dst<<">(";
				 x.getSrc()->accept(*this);
				 dst<<"))";
			 }
		 }

		 void visit(const HLSNodeSelectBits &x)
		 {
			 dst<<"(apint_get_range(";
			 x.getSrc()->accept(*this);
			 dst<<",";
			 dst<<x.getHi();
			 dst<<",";
			 dst<<x.getLo();
			 dst<<"))";
		 }

		 void visit(const HLSNodeSelect &x)
		 {
			 dst<<"(";
			 x.getCond()->accept(*this);
			 dst<<" ? ";
			 x.getTrueValue()->accept(*this);
			 dst<<" : ";
			 x.getFalseValue()->accept(*this);
			 dst<<")";
		 }

			void visit(const HLSNodeCat &x)
			 {
				dst<<"apint_concatenate(";
				x.getLeft()->accept(*this);
				dst<<",";
				x.getRight()->accept(*this);
				dst<<")";
			 }

	  void visit(const HLSNodeCallOutput &x)
	  {
	    dst<<x.getName();
	  }

	  void visit(const HLSNodeCall &x)
	  {
	    throw std::runtime_error("VisitExpr - Should be visiting HLSNodeCall here.");
	  }

		void visit(const HLSNodeOpaque &x)
	  {
			if(!m_ctxt.isTargetTool(x.getTool()))
				throw std::runtime_error("Received HLSNodeOpaque for wrong tool.");
			dst<<x.getValue();
		}

	};


    std::string HLSContext::strRep(const HLSExpr &w)
	{
    	OutputExpr oe(*this);
    	w.getNode()->accept(oe);
		return oe.str();
	}

    
    void HLSContext::emitInputParameter(
        const Operator &op,
        const Signal &sig
    ){
      dst<<"const "<<strRep(makeType(sig))<<" &"<<sig.getName();
    }
    
    void HLSContext::emitOutputParameter(
        const Operator &op,
	const Signal &sig
    ){
        *this<<strRep(makeType(sig))<<" &"<<sig.getName();
    }
   
    void HLSContext::emitSignature(
        const Operator &op
    ){
        *this<<"void "<<op.getUniqueName()<<"(";
        indent();
        
        auto signals=*op.getIOList();
        
        writeLine(" //inputs");
        for (unsigned i=0; i<signals.size(); i++){
            Signal* s = signals[i];
            if(s->type()==Signal::in){
                emitInputParameter(op, *s);
		if(i!=signals.size()-1)
		  write(",");
		writeLine();
            }
        }
        
        writeLine(" //outputs");
        for (unsigned i=0; i<signals.size(); i++){
            Signal* s = signals[i];
            if(s->type()==Signal::out){
                emitOutputParameter(op, *s);
		if(i!=signals.size()-1)
		  write(",");
		writeLine();
            }
        }
            
        
        unindent();
        dst<<")";
    }
    
    void HLSContext::emitDeclaration(
        const Operator &op
    ){
        if(m_functionDecls.find(op.getName())!=m_functionDecls.end())
            return;
        
        emitSignature(op);
        writeLine(";");
        writeLine();
        
        m_functionDecls.insert(op.getName());
    }
    
    void HLSContext::emitDefinition(
        const HLSOperator &op
    ){
        const Operator &bop=op.getOperator();
        
        if(m_functionDefs.find(bop.getName())!=m_functionDefs.end())
            return;
        
        writeLine(" // Sub-components count = %u", bop.getSubComponents().size());
	std::set<Operator*> suppress=op.suppressedHLSDefinitions();
        for(auto sc : bop.getSubComponents()){
            writeLine(" // decl %s", sc.second->getName().c_str());
	    if(suppress.find(sc.second)!=suppress.end()){
	      writeLine(" // ^^^ Supressed for HLS output");
	      continue;
	    }
            HLSOperator *hlsDecl=getHLSOperator(sc.second);
            emitDefinition(*hlsDecl);
            hlsDecl->releaseHLS();
        }
        
        if(bop.isRecirculatory()){
            throw std::runtime_error("Operator "+op.getOperator().getName()+" is recirculatory, cannot generate HLS.");
        }
        
        emitSignature(bop);
        writeLine("{");
        indent();
  
        {
			HLSScope scope(&bop, *this);
			op.emitHLSBody(*this, scope);
			unindent();
			scope.flush();
        }

        writeLine("}");
        writeLine();
    }   

};
