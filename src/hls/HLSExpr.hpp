#ifndef flopoco_hls_expr_hpp
#define flopoco_hls_expr_hpp

#include "HLSContext.hpp"
#include "HLSOperator.hpp"

namespace flopoco
{
    class HLSNode;
    typedef std::shared_ptr<HLSNode> HLSNodePtr;

    class HLSNode
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
	};


    class HLSNodeLessThanEquals
            : public HLSNodeCmp
	{
	public:
		HLSNodeLessThanEquals(const HLSNodePtr &left, const HLSNodePtr &right)
			: HLSNodeCmp(left, right)
		{}
	};

    class HLSNodeEquals
    		: public HLSNodeCmp
    {
	public:
		HLSNodeEquals(const HLSNodePtr &left, const HLSNodePtr &right)
			: HLSNodeCmp(left, right)
		{}
    };

    class HLSNodeLogicalOr
			: public HLSNodeCmp
	{
	public:
		HLSNodeLogicalOr(const HLSNodePtr &left, const HLSNodePtr &right)
			: HLSNodeCmp(left, right)
		{}
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

		static std::shared_ptr<HLSNodeCat> create(const HLSNodePtr &left, const HLSNodePtr &right)
		{ return std::make_shared<HLSNodeCat>(left, right); }
	};

    class HLSNodeSelectBits
			: public HLSNode
	{
    private:
    	HLSTypePtr select_type(const HLSTypePtr &a, int hi, int lo)
    	{
    		auto pa=std::dynamic_pointer_cast<HLSTypeInt>(a);
    		if(!pa) // TODO This should/must be relaxed
    			throw std::runtime_error("HLSNodeSelectBits - Argument must have integer type.");
    		if(hi>=a->getWidth())
    			throw std::runtime_error("HLSNodeSelectBits - hi index is out of range.");
    		if(lo<0)
    		    throw std::runtime_error("HLSNodeSelectBits - lo index is out of range.");
    		if(hi<lo)
    			throw std::runtime_error("HLSNodeSelectBits - Have hi<lo. Range must be at least one bit.");

    		return HLSTypeInt::create(false, hi-lo+1);
    	}
	public:
		HLSNodeSelectBits(const HLSNodePtr &a, int hi, int lo)
			: HLSNode(select_type(a->getType(), hi, lo))
		{}

		static std::shared_ptr<HLSNodeSelectBits> create(const HLSNodePtr &a, int hi, int lo)
		{ return std::make_shared<HLSNodeSelectBits>(a, hi, lo); }

		static std::shared_ptr<HLSNodeSelectBits> create(const HLSNodePtr &a, const std::string &x)
		{
			std::stringstream src(x);
			int hi, lo;
			std::string mid;
			src>>hi>>mid>>lo;
			if(mid!="downto")
				throw std::runtime_error("HLSNodeSelectBits - String range must be of form '<num> downto <num>'");

			return create(a, hi, lo);
		}
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

        void assign(const HLSNodePtr &o)
        {
            if(m_src)
                throw std::runtime_error("Output "+getName()+" has already been assigned.");
            m_src=o;
        }
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

        void assign(const HLSNodePtr &o)
        {
            if(m_src)
                throw std::runtime_error("Variable "+getName()+" has already been assigned.");
            m_src=o;
        }

        static std::shared_ptr<HLSNodeVar> create(const std::string &name, const HLSTypePtr &type)
		{ return std::make_shared<HLSNodeVar>(name, type); }
    };

    class HLSExpr
    {
    protected:
        HLSNodePtr m_base;
    public:
        HLSExpr(const HLSNodePtr &b)
    		: m_base(b)
    	{}

        const HLSNodePtr getNode() const
        { return m_base; }

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

        HLSExpr operator ||(const HLSExpr &o) const
        { return HLSExpr(std::make_shared<HLSNodeLogicalOr>(this->m_base, o.m_base)); }

        HLSExpr operator[](const std::string &idx) const
        {

        }


        //! Performs an assignment to a previously declared variable.
        HLSExpr operator=(const HLSExpr &o)
        {
            m_base->assign(o.m_base);
            return *this;
        }
    };

    inline HLSExpr select_if(const HLSExpr &cond, const HLSExpr &trueExpr, const HLSExpr &falseExpr)
    {
    	return HLSExpr(HLSNodeSelect::create(cond.getNode(), trueExpr.getNode(), falseExpr.getNode()));
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


}; // flopoco

#endif
