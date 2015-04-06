#include "HLSContext.hpp"
#include "HLSOperator.hpp"
#include "HLSScope.hpp"

namespace flopoco
{
    
	thread_local HLSScope *HLSScope::sm_pCurr=NULL;

	HLSScope::HLSScope(const Operator *op, HLSContext &ctxt)
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

	HLSScope &HLSScope::getAmbientScope()
	{
		HLSScope *pCurr=sm_pCurr;
		if(pCurr==NULL)
			throw std::runtime_error("HLSScope::getAmbientScope - no scope is active.");
		return *pCurr;
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

    HLSExpr hls_declare(const std::string &name, int width)
    {
    	return HLSExpr(HLSScope::getAmbientScope().declare(name, width));
    }

    HLSExpr hls_get(const std::string &name)
    {
    	return HLSExpr(HLSScope::getAmbientScope().get(name));
	}

};
