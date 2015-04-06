#include "HLSContext.hpp"
#include "HLSOperator.hpp"
#include "HLSScope.hpp"

#include "HLSExpr.hpp"

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
        if(sig.isFP()){
            return HLSTypeFloat::create(sig.wE(), sig.wF());
        }else if( sig.isBus()){
            return HLSTypeInt::create(false, sig.width());
        }else if(sig.width()==1 && !sig.isBus()){
            return HLSTypeBool::create();
        }else{
            throw std::runtime_error("Unknown signal type for signal "+sig.getName());
        }
    };

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

    std::string HLSContext::strRep(const HLSExpr &w)
	{
		throw std::runtime_error("HLSContext::strRep - Unsupported expression.");
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
        for(auto sc : bop.getSubComponents()){
            writeLine(" // decl %s", sc.second->getName().c_str());
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
  
        HLSScope scope(&bop, *this);
        op.emitHLSBody(*this, scope);
        unindent();
        writeLine("}");
        writeLine();
    }   

};
