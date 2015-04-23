#ifndef flopoco_hls_expr_hpp
#define flopoco_hls_expr_hpp

#include "hls/HLSTypes.hpp"

namespace flopoco
{
    class HLSNode;
    typedef std::shared_ptr<HLSNode> HLSNodePtr;

    class HLSOperator; // Declared elsewhere, needed for HLSNodeCall

    class HLSNodeCall;
    class HLSNodeCallOutput;
    class HLSNodeOpaque;
    class HLSNodeConstantInt;
    class HLSNodeAdd;
    class HLSNodeSub;
    class HLSNodeLessThan;
    class HLSNodeLessThanEquals;
    class HLSNodeEquals;
  class HLSNodeNotEquals; 
    class HLSNodeLogicalOr;
    class HLSNodeSelect;
    class HLSNodeCat;
    class HLSNodeSelectBits;
    class HLSNodeReinterpretBits;
    class HLSNodeInput;
    class HLSNodeOutput;
    class HLSNodeVar;

    class HLSNodeVisitor
    {
    public:
    	virtual ~HLSNodeVisitor()
    	{}

    	/*! Used if you want to provide some special case dynamic
    		stuff, like throwing on unhandled node types.
    		Return false to stop the default handler executing, true
    		to do whatever it would normally do for that node.
    		*/
    	virtual bool visitGeneric(const HLSNode &);

        virtual void visit(const HLSNodeCall &);
        virtual void visit(const HLSNodeCallOutput &);
    	virtual void visit(const HLSNodeOpaque &);
    	virtual void visit(const HLSNodeConstantInt &);
    	virtual void visit(const HLSNodeAdd &);
    	virtual void visit(const HLSNodeSub &);
    	virtual void visit(const HLSNodeLessThan &);
    	virtual void visit(const HLSNodeLessThanEquals &);
    	virtual void visit(const HLSNodeEquals &);
    	virtual void visit(const HLSNodeNotEquals &);
    	virtual void visit(const HLSNodeLogicalOr &);
    	virtual void visit(const HLSNodeSelect &);
    	virtual void visit(const HLSNodeCat &);
    	virtual void visit(const HLSNodeSelectBits &);
    	virtual void visit(const HLSNodeReinterpretBits &);
    	virtual void visit(const HLSNodeInput &);
    	virtual void visit(const HLSNodeOutput &);
    	virtual void visit(const HLSNodeVar &);
    };

    class HLSNode
      : public std::enable_shared_from_this<HLSNode>
    {
    private:
        HLSTypePtr m_type;
    protected:
        HLSNode(const HLSTypePtr &type)
            : m_type(type)
        {}
    public:
        virtual ~HLSNode()
        {}

        const HLSTypePtr &getType() const
        { return m_type; }

        virtual void assign(const HLSNodePtr &src)
        {
            throw std::runtime_error("HLSNode - Cannot assign this node type.");
        }

        virtual void accept(HLSNodeVisitor &visitor) const =0;
    };

  //! Represents outputs of invoked sub-operators
    class HLSNodeCallOutput;
  typedef std::shared_ptr<HLSNodeCallOutput> HLSNodeCallOutputPtr;

  class HLSNodeCall;
  typedef std::shared_ptr<HLSNodeCall> HLSNodeCallPtr;


    /*! This represents some opaque expression specific to a particular
     * output tool. This will only be used by operators that are some kind
     * of primitive that has to be specialised for the tool (e.g. ROMs)
     */
    class HLSNodeOpaque
        : public HLSNode
    {
    private:
        std::string m_value;
        std::string m_tool;
    public:
        HLSNodeOpaque(const HLSTypePtr &type, std::string value, std::string tool)
			: HLSNode(type)
			, m_value(value)
    		, m_tool(tool)
		{
		}

        const std::string &getValue() const
        { return m_value; }

        const std::string &getTool() const
        { return m_tool; }

        virtual void accept(HLSNodeVisitor &visitor) const
        { visitor.visit(*this); }

        static std::shared_ptr<HLSNodeOpaque> create(const HLSTypePtr &type, std::string value, std::string tool)
		{ return std::make_shared<HLSNodeOpaque>(type, value, tool); }
    };

    class HLSNodeConstantInt
        : public HLSNode
    {
    private:
        mpz_class m_value;
    public:
        HLSNodeConstantInt(const HLSTypePtr &type, const mpz_class &value)
			: HLSNode(type)
			, m_value(value)
		{
			auto p=std::dynamic_pointer_cast<HLSTypeInt>(type);
			if(!p)
				throw std::runtime_error("HLSNodeConstantInt - Type does not appear to be an int.");

			if(value < p->minValue() || p->maxValue()< value)
				throw std::runtime_error("HLSNodeConstantInt - Value is out of range for type.");
		}

        const mpz_class &getValue() const
        { return m_value; }

        virtual void accept(HLSNodeVisitor &visitor) const
        { visitor.visit(*this); }

        static std::shared_ptr<HLSNodeConstantInt> create(unsigned w, const mpz_class &value)
        { return create(HLSTypeInt::create(false,w), value); }

        static std::shared_ptr<HLSNodeConstantInt> create(const HLSTypePtr &type, const mpz_class &value)
		{ return std::make_shared<HLSNodeConstantInt>(type, value); }
    };

    class HLSNodeBinOp
        : public HLSNode
    {
    private:
        HLSNodePtr m_left;
        HLSNodePtr m_right;
    protected:
        HLSNodeBinOp(const HLSNodePtr &left, const HLSNodePtr &right, const HLSTypePtr &type)
            : HLSNode(type)
    		, m_left(left)
            , m_right(right)
        {}

        HLSNodeBinOp(const HLSNodePtr &left, const HLSNodePtr &right)
                    : HLSNode(same_type(left->getType(),right->getType()))
            		, m_left(left)
                    , m_right(right)
                {}

        HLSTypePtr same_type(const HLSTypePtr &left, const HLSTypePtr &right)
        {
        	if(!left->equals(right))
        		throw std::runtime_error("HLSNodeBinOp::same_type - Types must be the same.");
        	return left;
        }
    public:
        const HLSNodePtr &getLeft() const
        { return m_left; }

        const HLSNodePtr &getRight() const
        { return m_right; }
    };

    class HLSNodeAdd
        : public HLSNodeBinOp
    {
    public:
        HLSNodeAdd(const HLSNodePtr &left, const HLSNodePtr &right)
            : HLSNodeBinOp(left, right)
        {}

        virtual void accept(HLSNodeVisitor &visitor)  const
        { visitor.visit(*this); }

        static std::shared_ptr<HLSNodeAdd> create(const HLSNodePtr &left, const HLSNodePtr &right)
        { return std::make_shared<HLSNodeAdd>(left, right); }
    };

    class HLSNodeSub
        : public HLSNodeBinOp
    {
    public:
        HLSNodeSub(const HLSNodePtr &left, const HLSNodePtr &right)
            : HLSNodeBinOp(left, right)
        {}

        virtual void accept(HLSNodeVisitor &visitor)  const
        { visitor.visit(*this); }

        static std::shared_ptr<HLSNodeSub> create(const HLSNodePtr &left, const HLSNodePtr &right)
        { return std::make_shared<HLSNodeSub>(left, right); }
    };

    class HLSNodeCmp
       : public HLSNodeBinOp
	{
	public:
		HLSNodeCmp(const HLSNodePtr &left, const HLSNodePtr &right)
			: HLSNodeBinOp(left, right, HLSTypeBool::create())
		{}
	};

    class HLSNodeLessThan
            : public HLSNodeCmp
	{
	public:
		HLSNodeLessThan(const HLSNodePtr &left, const HLSNodePtr &right)
			: HLSNodeCmp(left, right)
		{}

        virtual void accept(HLSNodeVisitor &visitor)  const
        { visitor.visit(*this); }
	};


    class HLSNodeLessThanEquals
            : public HLSNodeCmp
	{
	public:
		HLSNodeLessThanEquals(const HLSNodePtr &left, const HLSNodePtr &right)
			: HLSNodeCmp(left, right)
		{}

        virtual void accept(HLSNodeVisitor &visitor)  const
        { visitor.visit(*this); }
	};

    class HLSNodeEquals
    		: public HLSNodeCmp
    {
	public:
		HLSNodeEquals(const HLSNodePtr &left, const HLSNodePtr &right)
			: HLSNodeCmp(left, right)
		{}

        virtual void accept(HLSNodeVisitor &visitor)  const
        { visitor.visit(*this); }
    };

    class HLSNodeNotEquals
    		: public HLSNodeCmp
    {
	public:
		HLSNodeNotEquals(const HLSNodePtr &left, const HLSNodePtr &right)
			: HLSNodeCmp(left, right)
		{}

        virtual void accept(HLSNodeVisitor &visitor)  const
        { visitor.visit(*this); }
    };

    class HLSNodeLogicalOr
			: public HLSNodeCmp
	{
	public:
		HLSNodeLogicalOr(const HLSNodePtr &left, const HLSNodePtr &right)
			: HLSNodeCmp(left, right)
		{}

        virtual void accept(HLSNodeVisitor &visitor)  const
        { visitor.visit(*this); }
	};

    class HLSNodeSelect
        : public HLSNode
    {
    private:
    	HLSNodePtr m_cond, m_trueValue, m_falseValue;
    public:
        HLSNodeSelect(const HLSNodePtr &cond, const HLSNodePtr &trueValue, const HLSNodePtr &falseValue)
            : HLSNode(trueValue->getType())
    		, m_cond(cond)
    		, m_trueValue(trueValue)
    		, m_falseValue(falseValue)
        {
        	if(!trueValue->getType()->equals(falseValue->getType()))
        		throw std::runtime_error("HLSNodeSelect - true and false have different types.");

        }

        const HLSNodePtr &getCond() const
        { return m_cond; }

        const HLSNodePtr &getTrueValue() const
        { return m_trueValue; }

        const HLSNodePtr &getFalseValue() const
        { return m_falseValue; }

        virtual void accept(HLSNodeVisitor &visitor)  const
        { visitor.visit(*this); }

        static std::shared_ptr<HLSNodeSelect> create(const HLSNodePtr &cond, const HLSNodePtr &trueValue, const HLSNodePtr &falseValue)
        { return std::make_shared<HLSNodeSelect>(cond,trueValue,falseValue); }
    };

    class HLSNodeCat
			: public HLSNodeBinOp
	{
    private:
    	HLSTypePtr cat_type(const HLSTypePtr &a, const HLSTypePtr &b)
    	{
    		auto pa=std::dynamic_pointer_cast<HLSTypeInt>(a);
    		auto pb=std::dynamic_pointer_cast<HLSTypeInt>(b);
    		if(!pa || !pb)
    			throw std::runtime_error("HLSNodeCat - Both arguments must have integer type.");
    		if(pa->isSigned() || pb->isSigned()) // Is this reasonable? What should semantics be?
    			throw std::runtime_error("HLSNodeCat - Both arguments must be unsigned.");
    		return HLSTypeInt::create(false, pa->getWidth()+pb->getWidth());
    	}
	public:
		HLSNodeCat(const HLSNodePtr &left, const HLSNodePtr &right)
			: HLSNodeBinOp(left, right, cat_type(left->getType(), right->getType()))
		{}

        virtual void accept(HLSNodeVisitor &visitor)  const
        { visitor.visit(*this); }

		static std::shared_ptr<HLSNodeCat> create(const HLSNodePtr &left, const HLSNodePtr &right)
		{ return std::make_shared<HLSNodeCat>(left, right); }
	};

    class HLSNodeSelectBits
			: public HLSNode
	{
    private:
    	HLSNodePtr m_src;
    	int m_hi, m_lo;

    	HLSTypePtr select_type(const HLSTypePtr &a, int hi, int lo)
    	{
    		auto pa=std::dynamic_pointer_cast<HLSTypeInt>(a);
    		if(!pa) // TODO This should/must be relaxed
    			throw std::runtime_error(std::string("HLSNodeSelectBits - Argument must have integer type, it has type.")+a->getName());

    		if(hi>=a->getWidth())
    			throw std::runtime_error("HLSNodeSelectBits - hi index is out of range.");
    		if(lo<0)
    		    throw std::runtime_error("HLSNodeSelectBits - lo index is out of range.");
    		if(hi<lo)
    			throw std::runtime_error("HLSNodeSelectBits - Have hi<lo. Range must be at least one bit.");

    		// The return type is always unsigned... just because
    		return HLSTypeInt::create(false, hi-lo+1);
    	}
	public:
		HLSNodeSelectBits(const HLSNodePtr &a, int hi, int lo)
			: HLSNode(select_type(a->getType(), hi, lo))
			, m_src(a)
			, m_hi(hi)
			, m_lo(lo)
		{}

		HLSNodePtr getSrc() const
		{ return m_src; }

		int getHi() const
		{ return m_hi; }

		int getLo() const
		{ return m_lo; }

        virtual void accept(HLSNodeVisitor &visitor)  const
        { visitor.visit(*this); }

		static std::shared_ptr<HLSNodeSelectBits> create(const HLSNodePtr &a, int hi, int lo);

		static std::shared_ptr<HLSNodeSelectBits> create(const HLSNodePtr &a, int index);

		static std::shared_ptr<HLSNodeSelectBits> create(const HLSNodePtr &a, const std::string &x);
	};

    class HLSNodeReinterpretBits
			: public HLSNode
	{
	private:
		HLSNodePtr m_src;

	public:
		HLSNodeReinterpretBits(const HLSNodePtr &a, const HLSTypePtr &target)
			: HLSNode(target)
			, m_src(a)
		{
			if(a->getType()->getWidth() != target->getWidth())
				throw std::runtime_error("Cannot reinterpret to different bit widths");
		}

		HLSNodePtr getSrc() const
		{ return m_src; }

		virtual void accept(HLSNodeVisitor &visitor)  const
		{ visitor.visit(*this); }

		static std::shared_ptr<HLSNodeReinterpretBits> create(const HLSNodePtr &a, const HLSTypePtr &target)
		{ return std::make_shared<HLSNodeReinterpretBits>(a, target); }

	};

    class HLSNodeDecl
        : public HLSNode
    {
    private:
        std::string m_name;
    protected:
        HLSNodeDecl(const std::string &name, const HLSTypePtr &type)
            : HLSNode(type)
            , m_name(name)
        {}
    public:
        const std::string &getName() const
        { return m_name; }

        virtual HLSNodePtr getSrc() const =0;

        virtual bool isDefined() const =0;
    };
    typedef std::shared_ptr<HLSNodeDecl> HLSNodeDeclPtr;

    class HLSNodeInput
        : public HLSNodeDecl
    {
    public:
        HLSNodeInput(const std::string &name, const HLSTypePtr &type)
            : HLSNodeDecl(name, type)
        {}

        virtual bool isDefined() const
        { return true; }

        virtual void accept(HLSNodeVisitor &visitor)  const
        { visitor.visit(*this); }

        HLSNodePtr getSrc() const
        { return HLSNodePtr(); }
    };

    class HLSNodeOutput
        : public HLSNodeDecl
    {
    protected:
        HLSNodePtr m_src;
    public:
        HLSNodeOutput(const std::string &name, const HLSTypePtr &type)
			: HLSNodeDecl(name, type)
		{}

        virtual bool isDefined() const
        { return !!m_src; }

        HLSNodePtr getSrc() const
        { return m_src; }

        void assign(const HLSNodePtr &o)
        {
            if(m_src)
                throw std::runtime_error("Output "+getName()+" has already been assigned.");
            m_src=o;
        }

        virtual void accept(HLSNodeVisitor &visitor)  const
        { visitor.visit(*this); }
    };

    class HLSNodeCallOutput
        : public HLSNodeDecl
    {
    protected:
      HLSNodeCallPtr m_call;
      std::string m_innerName;
    public:
      HLSNodeCallOutput(const std::string &name, const HLSTypePtr &type, HLSNodeCallPtr call, std::string innerName)
			: HLSNodeDecl(name, type)
			, m_call(call)
			, m_innerName(innerName)
		{}

        virtual bool isDefined() const
        { return true; }

      virtual HLSNodePtr getSrc() const
      { return HLSNodePtr(); }

        HLSNodeCallPtr getCall() const
        { return m_call; }

      std::string getInnerName() const
      { return m_innerName; }

        virtual void accept(HLSNodeVisitor &visitor)  const
        { visitor.visit(*this); }
    };

    class HLSNodeVar
        : public HLSNodeDecl
    {
    protected:
        bool m_defined;
        HLSNodePtr m_src;
    public:
        HLSNodeVar(const std::string &name, const HLSTypePtr &type)
			: HLSNodeDecl(name, type)
		{}

        virtual bool isDefined() const
        { return !!m_src; }

        HLSNodePtr getSrc() const
        { return m_src; }

        void assign(const HLSNodePtr &o)
        {
            if(m_src)
                throw std::runtime_error("Variable "+getName()+" has already been assigned.");
            m_src=o;
        }

        virtual void accept(HLSNodeVisitor &visitor)  const
        { visitor.visit(*this); }

        static std::shared_ptr<HLSNodeVar> create(const std::string &name, const HLSTypePtr &type)
		{ return std::make_shared<HLSNodeVar>(name, type); }
    };

    //! Represents a call to a sub-operator
    /*! This is kind of wierd, as it exposes zero or more sub-nodes representing the outputs */
    class HLSNodeCall
      : public HLSNodeDecl
    {
    private:
      const HLSOperator *m_op;
      std::map<std::string,HLSNodePtr> m_inputs;
      std::map<std::string,HLSNodeCallOutputPtr> m_outputs;
    public:
      HLSNodeCall(const HLSOperator *op, std::string name, std::map<std::string,HLSNodePtr> inputs, std::map<std::string,std::string> outputs);

      ~HLSNodeCall();

      const HLSOperator *getOperator() const
      { return m_op; }

      unsigned getInputCount() const
      { return m_inputs.size(); }

      HLSNodePtr getInput(std::string name) const
      { return m_inputs.at(name); }

      //! Get input as defined by the underlying Operator's view
      HLSNodePtr getInput(unsigned index) const;

      unsigned getOutputCount() const
      { return m_outputs.size(); }

      HLSNodeCallOutputPtr getOutput(std::string name) const
      { return m_outputs.at(name); }

        virtual bool isDefined() const
        { return true; }

        HLSNodePtr getSrc() const
        { return HLSNodePtr(); }

      virtual void accept(HLSNodeVisitor &visitor) const
      {
	for(auto x : m_inputs){
	  HLSNodePtr i=m_inputs.at(x.first);
	  i->accept(visitor);
	}
	for(auto x : m_outputs){
	  HLSNodePtr o=m_outputs.at(x.first);
	  o->accept(visitor);
	}
	visitor.visit(*this);
      }

      static std::shared_ptr<HLSNodeCall> create(const HLSOperator *op, std::string name, std::map<std::string,HLSNodePtr> inputs, std::map<std::string,std::string> outputs)
      { return std::make_shared<HLSNodeCall>(op, name, inputs, outputs); }
    };


    class HLSExpr
    {
    protected:
        HLSNodePtr m_base;
    public:
      //! "Empty" expression. Any use will throw an error
      /*! Mainly here to make things like STL containers easier */
      HLSExpr()
      {}

        HLSExpr(const HLSNodePtr &b)
    		: m_base(b)
    	{}

        const HLSNodePtr getNode() const
        { return m_base; }

        const HLSTypePtr getType() const
        {
	  if(!m_base)
	    throw std::runtime_error("HLSExpr::getType - expression is null.");
	  return m_base->getType();
        }

        int getWidth() const
        { return getType()->getWidth(); }

        HLSExpr operator +(const HLSExpr &o) const
        { return HLSExpr(std::make_shared<HLSNodeAdd>(this->m_base, o.m_base)); }

        HLSExpr operator +(int o) const
        {
        	return HLSExpr(
        			HLSNodeAdd::create(
        					this->m_base,
        					HLSNodeConstantInt::create(m_base->getType(), mpz_class(o))
        			)
        	);
        }

        HLSExpr operator -(const HLSExpr &o) const
        { return HLSExpr(std::make_shared<HLSNodeSub>(this->m_base, o.m_base)); }

        HLSExpr operator -(int o) const
        { return HLSExpr(HLSNodeSub::create(this->m_base, HLSNodeConstantInt::create(m_base->getType(), mpz_class(o)))); }

        HLSExpr operator <(const HLSExpr &o) const
        { return HLSExpr(std::make_shared<HLSNodeLessThan>(this->m_base, o.m_base)); }

        HLSExpr operator <=(const HLSExpr &o) const
        { return HLSExpr(std::make_shared<HLSNodeLessThanEquals>(this->m_base, o.m_base)); }

        HLSExpr operator ==(const HLSExpr &o) const
        { return HLSExpr(std::make_shared<HLSNodeEquals>(this->m_base, o.m_base)); }

        HLSExpr operator !=(const HLSExpr &o) const
        { return HLSExpr(std::make_shared<HLSNodeNotEquals>(this->m_base, o.m_base)); }

        HLSExpr operator ||(const HLSExpr &o) const
        { return HLSExpr(std::make_shared<HLSNodeLogicalOr>(this->m_base, o.m_base)); }

        HLSExpr operator[](const std::string &idx) const
        {
        	return HLSExpr(HLSNodeSelectBits::create(this->m_base, idx));
        }

        HLSExpr operator[](int idx) const
        {
        	return HLSExpr(HLSNodeSelectBits::create(this->m_base, idx));
        }


        //! Performs an assignment to a previously declared variable.
        HLSExpr operator=(const HLSExpr &o)
        {
	  if(!m_base)
	    throw std::runtime_error("Attempt to assign to uninitialised (empty) HLSExpr.");
            m_base->assign(o.m_base);
            return *this;
        }
    };

  HLSExpr select_if(const HLSExpr &cond0, const HLSExpr &expr0,
		    const HLSExpr &exprFalse);

  HLSExpr select_if(const HLSExpr &cond0, const HLSExpr &expr0,
		    const HLSExpr &cond1, const HLSExpr &expr1,
		    const HLSExpr &exprFalse);

  HLSExpr select_if(const HLSExpr &cond0, const HLSExpr &expr0,
		    const HLSExpr &cond1, const HLSExpr &expr1,
		    const HLSExpr &cond2, const HLSExpr &expr2,
		    const HLSExpr &exprFalse);

    inline HLSExpr reinterpret_bits(const HLSExpr &x, const HLSTypePtr &t)
	{
		return HLSExpr(HLSNodeReinterpretBits::create(x.getNode(), t));
	}

    inline HLSExpr cat(const HLSExpr &a, const HLSExpr &b)
    {
    	return HLSExpr(HLSNodeCat::create(a.getNode(),b.getNode()));
    }

    inline HLSExpr cat(const HLSExpr &a, const HLSExpr &b, const HLSExpr &c)
    {
    	return cat(a.getNode(),cat(b.getNode(),c.getNode()));
    }

    inline HLSExpr hls_cg(unsigned width, const mpz_class &value)
    {
    	return HLSExpr(HLSNodeConstantInt::create(width, value));
    }

        inline HLSExpr hls_og(unsigned width)
	{
		mpz_class tmp(1);
		tmp<<width;
		return hls_cg(width, tmp-1);
	}

	inline HLSExpr hls_zg(unsigned width)
	{
		return hls_cg(width, mpz_class(0));
	}

	inline HLSExpr hls_xg(unsigned width)
	{
		return hls_og(width);
	}

  HLSExpr hls_call(const Operator *op, std::string name, HLSExpr arg0, const std::string &resName);

  HLSExpr hls_call(const HLSOperator *op, std::string name, HLSExpr arg0, const std::string &resName);

  HLSExpr hls_call(const HLSOperator *op, std::string name, const std::map<std::string,HLSExpr> &args, const std::string &resName);

  HLSExpr hls_call(const HLSOperator *op, std::string name, const std::map<std::string,HLSExpr> &args, const std::string &resName);

  void hls_call(const HLSOperator *op, std::string name, const std::map<std::string,HLSExpr> &inputs, const std::map<std::string,std::string> &outputs);

  void hls_call(const HLSOperator *op, std::string name, const std::map<std::string,HLSExpr> &inputs, const std::map<std::string,std::string> &outputs);

}; // flopoco

#endif
