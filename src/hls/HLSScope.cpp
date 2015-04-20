#include "HLSContext.hpp"
#include "HLSOperator.hpp"
#include "HLSScope.hpp"

#include <set>
#include <cassert>

namespace flopoco
{
    
	thread_local HLSScope *HLSScope::sm_pCurr=NULL;

	HLSScope::HLSScope(const Operator *op, HLSContext &ctxt)
		: m_ctxt(ctxt)
		, m_flushed(false)
	{
		if(sm_pCurr!=NULL) // No idea if this is bad, can't be bothered thinking about it.
			throw std::runtime_error("HLSScope - Nested scopes not supported.");
		m_pParent=sm_pCurr;
		sm_pCurr=this;

		auto signals=op->getIOList();

		for (unsigned i=0; i<signals->size(); i++){
			Signal* s = signals->at(i);
			HLSNodeDeclPtr node;
			if(s->type()==Signal::in){
				node=std::make_shared<HLSNodeInput>(s->getName(), ctxt.makeType(*s));
			}else if(s->type()==Signal::out){
				node=std::make_shared<HLSNodeOutput>(s->getName(), ctxt.makeType(*s));
			}
			if(node){
				m_decls[s->getName()]=node;
			}
		}
	}

	class OutputDecls
		: public HLSNodeVisitor
	{
	private:
		std::set<std::string> m_done;
		std::set<std::string> m_pending;
	public:
		OutputDecls()
		{}

		void define(const HLSNodeDecl& n)
		{
			if(dynamic_cast<const HLSNodeInput*>(&n))
				return; // always defined

			if(m_done.find(n.getName())!=m_done.end())
				return;
			if(m_pending.find(n.getName())!=m_pending.end())
				throw std::runtime_error("There appears to be recursive definition for "+n.getName());

			m_pending.insert(n.getName());

			if(!n.isDefined())
				throw std::runtime_error("Declaration "+n.getName()+" has not been given a definition.");

			auto pSrc=n.getSrc();
			assert(pSrc);
			pSrc->accept(*this);

			auto pC=dynamic_cast<const HLSNodeCall*>(&n);
			auto pV=dynamic_cast<const HLSNodeVar*>(&n);
			auto pO=dynamic_cast<const HLSNodeOutput*>(&n);

			if(pC){
			  HLSScope::getAmbientContext().emitCall(*pC);
			}else if(pV){
				HLSScope::getAmbientContext().emitDeclareAndAssign(*pV);

			}else if(pO){
				HLSScope::getAmbientContext().emitAssignOutput(*pO);

			}else{
				throw std::runtime_error("Unknown decl type.");
			}

			auto it=m_pending.find(n.getName());
			assert(it!=m_pending.end());
			m_pending.erase(it);
			m_done.insert(n.getName());
		}

		 void visit(const HLSNodeOutput &x)
		 {
			 // Maybe people use outputs in expressions? Should they?
			throw std::runtime_error("Shouldn't be possible to visit an output (?).");
		 }

		 void visit(const HLSNodeVar &x)
		 {
			 //std::cerr<<"// Visiting : "<<x.getName()<<"\n";
			 define(x);
		 }

	  void visit(const HLSNodeCall &x)
	  {
	    define(x);
	  }

	  void visit(const HLSNodeCallOutput &x)
	  {
	    x.getCall()->accept(*this);
	  }
	};

	void HLSScope::flush()
	{
		OutputDecls visitor;

		// Have to flush everything out to the context
		for(auto x : m_decls){
			//std::cerr<<"// define : "<<x.first<<"\n";
			visitor.define(*x.second);
		}

		m_flushed=true;
	}

	HLSScope::~HLSScope()
	{
		// Can't throw here
		assert(m_flushed);

		sm_pCurr=m_pParent;
	}

	HLSScope &HLSScope::getAmbientScope()
	{
		HLSScope *pCurr=sm_pCurr;
		if(pCurr==NULL)
			throw std::runtime_error("HLSScope::getAmbientScope - no scope is active.");
		return *pCurr;
	}

	HLSContext &HLSScope::getAmbientContext()
	{
		return getAmbientScope().m_ctxt;
	}

	HLSNodeDeclPtr HLSScope::get(const std::string &x)
	{
		auto it=m_decls.find(x);
		if(it==m_decls.end())
			throw std::runtime_error("No node called '"+x+"'");
		return it->second;
	}

	HLSNodeDeclPtr HLSScope::declare(const std::string &x, unsigned w)
	{
		auto it=m_decls.find(x);
		if(it!=m_decls.end())
			throw std::runtime_error("A node called '"+x+"' already exists.");

		auto res=HLSNodeVar::create(x, HLSTypeInt::create(false, w));
		m_decls[x]=res;
		return res;
	}

  HLSNodeDeclPtr HLSScope::call(const HLSOperator *hop, std::string name, const std::map<std::string,HLSExpr> &inputs, const std::map<std::string,std::string> &outputs)
  {
    const Operator &op=hop->getOperator();

    // Make sure the instance name is unique
		auto it=m_decls.find(name);
		if(it!=m_decls.end())
			throw std::runtime_error("A node called '"+name+"' already exists.");

		// Walk through the operator inputs, make sure all have a binding
		for(int i=0;i<op.getNumberOfInputs();i++){
		  const Signal *s=op.getInputSignal(i);
		  if(inputs.find(s->getName())==inputs.end())
		    throw std::runtime_error("HLSScope::call - No binding for input "+s->getName());
		}

		// Walk through the operator ouputs, make sure they all have a declaration to go to
		for(int i=0;i<op.getNumberOfOutputs();i++){
		  const Signal *s=op.getOutputSignal(i);
		  if(inputs.find(s->getName())==inputs.end())
		    throw std::runtime_error("HLSScope::call - No binding for output "+s->getName());
		}

		// Check all given inputs exist on operator, and convert to nodes
		std::map<std::string,HLSNodePtr> inputsN;
		for(auto x : inputs){
		  const Signal *s=op.getInputSignal(x.first);
		  if(!s)
		    throw std::runtime_error("HLSScope::call - No input called "+x.first);
		  inputsN[x.first]=x.second.getNode();
		}

		// Create the declaration representing the instance
		auto res=HLSNodeCall::create(hop, name, inputsN, outputs);
		m_decls[name]=res;
		
		// Bind the outputs to declarations, also checking the signals are unique
		for(auto x : outputs){
		  const Signal *s=op.getOutputSignal(x.first);
		  if(!s)
		    throw std::runtime_error("HLSScope::call - No output called "+x.first);
		  		  
		  m_decls[x.second]=res->getOutput(x.first);
		}
		return res;
  }

    HLSExpr hls_declare(const std::string &name, int width)
    {
    	return HLSExpr(HLSScope::getAmbientScope().declare(name, width));
    }

    HLSExpr hls_get(const std::string &name)
    {
    	return HLSExpr(HLSScope::getAmbientScope().get(name));
	}

  HLSExpr hls_call(const HLSOperator *hop, HLSExpr arg0, const std::string &resName)
  {
    const Operator &op=hop->getOperator();

    if(op.getNumberOfInputs()!=1)
      throw std::runtime_error("hls_call - There is not exactly one input.");
    std::string innerArgName=op.getInputSignal(0)->getName();

    std::map<std::string,HLSExpr> inputs;
    inputs[innerArgName]=arg0;
    return hls_call(hop, inputs, resName);
  }

  HLSExpr hls_call(const HLSOperator *hop, const std::map<std::string,HLSExpr> &args, const std::string &resName)
  {
    const Operator &op=hop->getOperator();

    if(op.getNumberOfOutputs()!=1)
      throw std::runtime_error("hls_call - There is not exactly one output.");
    std::string innerResName=op.getOutputSignal(0)->getName();
    
    std::map<std::string,std::string> outputs;
    outputs[innerResName]=resName;
    hls_call(hop, args, outputs);
    return hls_get(resName);
  }

  void hls_call(const HLSOperator *op, const std::map<std::string,HLSExpr> &inputs, const std::map<std::string,std::string> &outputs)
  {    
    throw std::runtime_error("hls_call - Not implemented.");
  }


};
