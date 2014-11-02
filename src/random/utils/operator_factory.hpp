#ifndef flopoco_random_utils_operator_factory_hpp
#define flopoco_random_utils_operator_factory_hpp

#include "Operator.hpp"

#include <boost/smart_ptr.hpp>

namespace flopoco
{
namespace random
{

class OperatorFactory;
typedef boost::shared_ptr<OperatorFactory> OperatorFactoryPtr;

class OperatorFactory
{
private:
	static std::vector<OperatorFactoryPtr> sm_factoriesByIndex;
	static std::map<std::string,OperatorFactoryPtr> sm_factoriesByName;
public:
	virtual const std::string &Name() const =0;
	
	virtual const std::vector<std::string> &Categories() const=0;
	
	virtual void Usage(std::ostream &dst) const=0;
	
	/*! Consumes zero or more string arguments, and creates an operator
		\param args The offered arguments start at index 0 of the vector, and it is up to the
				factory to check the types and whether there are enough.
		\param consumed On exit, the factory indicates how many of the arguments are used up.
		\param inputDelayMap Delays used during construction
		\note I would prefer this to be shared_ptr<Operator>, but this is more in the
		flopoco style. If desired, the creator can immediately wrap it in shared_ptr, as it would
		be a unique reference.
	*/
	virtual Operator *Parse(
		Target *target,
		const std::vector<std::string> &args,
		int &consumed
	) const =0;

	// Returns a set of valid parameters. For each list of parameters the component should build and be testable
	virtual std::vector<std::vector<std::string> > GetValidTestParameters() const
	{ return std::vector<std::vector<std::string> >(); }

	static void RegisterOperatorFactory(OperatorFactoryPtr factory);

	static unsigned GetFactoryCount();
	static OperatorFactoryPtr GetFactoryByIndex(unsigned i);

	static OperatorFactoryPtr FindFactory(std::string operatorName);

	//! For compatibility with flopoco command line interface
	static void classic_usage(char *name, string opName);
	
	//! For compatibility with flopco command line interface
	static bool classic_parseCommandLine(int argc, char* argv[], Target *target, std::string opname, int &i);
	
	//! Emulates the OP macro from main.cpp, with appropriate colouring
	static void classic_OP(std::ostream &dst, std::string name, std::string args, bool newOperator);
};

class DefaultOperatorFactory
{
public:
	// Note: I'm not using boost::function here, as it's likely to scare people, and also drags in quite a few header dependencies
	typedef void (*usage_func_t)(std::ostream &);
	typedef Operator * (*parser_func_t)(Target *,const std::vector<std::string> &,int &);	
private:
	class SimpleOperatorFactory
		: public OperatorFactory
	{
	private:
		std::string m_name;
		std::vector<std::string> m_categories;
		usage_func_t m_usage;
		parser_func_t m_parser;
		std::vector<std::vector<std::string> > m_testParams;
	public:
		SimpleOperatorFactory(
			std::string name,			
			std::string categories,
			usage_func_t usage,
			parser_func_t parser,
			const std::vector<std::vector<std::string> > &testParams
		);
		
		virtual const std::string &Name() const
		{ return m_name; }
		
		virtual const std::vector<std::string> &Categories() const
		{ return m_categories; }
		
		virtual void Usage(std::ostream &dst) const
		{ m_usage(dst); }
		
		virtual Operator *Parse(
			Target *target,
			const std::vector<std::string> &args,
			int &consumed
		)const
		{ return m_parser(target, args, consumed); }
	};
public:	
	static std::vector<std::string> Parameters(std::string a)
	{ return std::vector<std::string>(1, a); }
	
	static std::vector<std::string> Parameters(std::string a, std::string b)
	{ std::vector<std::string> res; res.push_back(a); res.push_back(b); return res; }
	
	static std::vector<std::string> Parameters(std::string a, std::string b, std::string c)
	{ std::vector<std::string> res; res.push_back(a); res.push_back(b); res.push_back(c); return res; }
	
	static std::vector<std::string> Parameters(std::string a, std::string b, std::string c, std::string d)
	{ std::vector<std::string> res; res.push_back(a); res.push_back(b); res.push_back(c); res.push_back(d); return res; }
	
	static std::vector<std::string> Parameters(std::string a, std::string b, std::string c, std::string d, std::string e)
	{ std::vector<std::string> res; res.push_back(a); res.push_back(b); res.push_back(c); res.push_back(d); res.push_back(e); return res; }
	
	
	static std::vector<std::vector<std::string> > ParameterList(const std::vector<std::string> &a)
	{ std::vector<std::vector<std::string> > res; res.push_back(a); return res; }

	static std::vector<std::vector<std::string> > ParameterList(const std::vector<std::string> &a, const std::vector<std::string> &b)
	{ std::vector<std::vector<std::string> > res; res.push_back(a); res.push_back(b); return res; }
	
	static std::vector<std::vector<std::string> > ParameterList(const std::vector<std::string> &a, const std::vector<std::string> &b, const std::vector<std::string> &c)
	{ std::vector<std::vector<std::string> > res; res.push_back(a); res.push_back(b); res.push_back(c); return res; }
	
	/** Implements a no-frills factory, given a usage function and a parser function
		\param name Name for the operator. The factory will only respond for this name (case sensitive)
		\param categories A semi-colon seperated list of categories that the operator belongs to.
		\param usage A function that can print the usage for the function to a std::ostream
		\param parser A function that can parse a vector of string arguments into an Operator instance
		\param testParameters Zero or more sets of arguments that the operator can be tested with.
	**/
	static void Register(
		std::string name,			
		std::string categories,	// semi-colon seperated list of categories
		usage_func_t usage,
		parser_func_t parser,
		const std::vector<std::vector<std::string> > &testParameters=std::vector<std::vector<std::string> >()
	);
};

}; // random
}; // flopoco

#endif
