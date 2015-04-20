#include "HLSExpr.hpp"

#include "HLSOperator.hpp"

namespace flopoco{

	bool HLSNodeVisitor::visitGeneric(const HLSNode &c)
	{
		return true;
	}

	 void HLSNodeVisitor::visit(const HLSNodeConstantInt &x)
	{
		 visitGeneric(x);
		// do nothing
	}

	 void HLSNodeVisitor::visit(const HLSNodeOpaque &x)
	{
		 visitGeneric(x);
	}

	 void HLSNodeVisitor::visit(const HLSNodeAdd &x)
	{
		 if(!visitGeneric(x))
			 return ;
		x.getLeft()->accept(*this);
		x.getRight()->accept(*this);
	}

	void HLSNodeVisitor::visit(const HLSNodeSub &x)
	{
		if(!visitGeneric(x))
					 return ;
		x.getLeft()->accept(*this);
		x.getRight()->accept(*this);
	}

	 void HLSNodeVisitor::visit(const HLSNodeLessThan &x)
	 {
		 if(!visitGeneric(x))
		 			 return ;
		x.getLeft()->accept(*this);
		x.getRight()->accept(*this);
	}

	 void HLSNodeVisitor::visit(const HLSNodeLessThanEquals &x)
	 {
		 if(!visitGeneric(x))
		 			 return ;
		x.getLeft()->accept(*this);
		x.getRight()->accept(*this);
	}

	 void HLSNodeVisitor::visit(const HLSNodeEquals &x)
	 {
		 if(!visitGeneric(x))
		 			 return ;
		x.getLeft()->accept(*this);
		x.getRight()->accept(*this);
	 }

	 void HLSNodeVisitor::visit(const HLSNodeLogicalOr &x)
	 {
		 if(!visitGeneric(x))
			 return ;
		x.getLeft()->accept(*this);
		x.getRight()->accept(*this);
	 }

	 void HLSNodeVisitor::visit(const HLSNodeSelect &x)
	 {
		 if(!visitGeneric(x))
		 			 return ;
		x.getCond()->accept(*this);
		x.getTrueValue()->accept(*this);
		x.getFalseValue()->accept(*this);
	 }

	 void HLSNodeVisitor::visit(const HLSNodeCat &x)
	 {
		 if(!visitGeneric(x))
		 			 return ;
		x.getLeft()->accept(*this);
		x.getRight()->accept(*this);
	 }

	 void HLSNodeVisitor::visit(const HLSNodeSelectBits &x)
	 {
		 if(!visitGeneric(x))
		 			 return ;
		x.getSrc()->accept(*this);
	 }

	 void HLSNodeVisitor::visit(const HLSNodeReinterpretBits &x)
	 {
		 if(!visitGeneric(x))
		 			 return ;
		x.getSrc()->accept(*this);
	 }

	 void HLSNodeVisitor::visit(const HLSNodeInput &x)
	 {
		 if(!visitGeneric(x))
		 			 return ;
		// do nothing
	 }

	 void HLSNodeVisitor::visit(const HLSNodeOutput &x)
	 {
		 if(!visitGeneric(x))
		 			 return ;
		if(x.getSrc())
			x.getSrc()->accept(*this);
	 }

	 void HLSNodeVisitor::visit(const HLSNodeVar &x)
	 {
		 if(!visitGeneric(x))
		 			 return ;
		if(x.getSrc())
			x.getSrc()->accept(*this);
	 }

  HLSNodeCall::HLSNodeCall(const HLSOperator *op, std::string name, std::map<std::string,HLSNodePtr> inputs, std::map<std::string,std::string> outputs)
	: HLSNodeDecl(name, HLSTypeVoid::create())
	, m_op(op)
	, m_inputs(inputs)
      {
	std::shared_ptr<HLSNodeCall> pMe=std::dynamic_pointer_cast<HLSNodeCall>(shared_from_this());

	const Operator &pop=op->getOperator();

	for(auto x : inputs){
	  const Signal *s=pop.getInputSignal(x.first);
	  if(!s)
	    throw std::runtime_error("HLSNodeCall::HLSNodeCall - No input called "+x.first+" on operator "+pop.getName());
	}
	for(auto x : outputs){
	  const Signal *s=pop.getOutputSignal(x.first);
	  if(!s)
	    throw std::runtime_error("HLSNodeCall::HLSNodeCall - No output called "+x.first+" on operator "+pop.getName());
	  HLSTypePtr type=makeHLSType(*s);
	  m_outputs[x.first]=std::make_shared<HLSNodeCallOutput>(x.second, type, pMe, x.first);
	}
	
      }

	 //////////////////////////////////////////////////
	 // HLSNodeSelectBits

		std::shared_ptr<HLSNodeSelectBits> HLSNodeSelectBits::create(const HLSNodePtr &a, int hi, int lo)
		{
			// Unlike the constructor, this will insert a cast to
			// an integer type if necessary
			auto pa=std::dynamic_pointer_cast<HLSTypeInt>(a);
			if(pa)
				return std::make_shared<HLSNodeSelectBits>(a, hi, lo);

			auto rt=HLSTypeInt::create(false, a->getType()->getWidth());
			return std::make_shared<HLSNodeSelectBits>(
					HLSNodeReinterpretBits::create(a,rt), hi, lo
			);
		}


		std::shared_ptr<HLSNodeSelectBits> HLSNodeSelectBits::create(const HLSNodePtr &a, const std::string &x)
		{
			std::stringstream src(x);
			int hi, lo;
			char pre, post;
			std::string mid;
			src>>pre>>hi>>mid>>lo>>post;
			if(pre!='(' || mid!="downto" || post!=')')
				throw std::runtime_error("HLSNodeSelectBits - String range must be of form '(<num> downto <num>)'");

			return create(a, hi, lo);
		}



};

