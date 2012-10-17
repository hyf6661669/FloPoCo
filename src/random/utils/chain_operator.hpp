#ifndef flopoco_random_utils_chain_operator_hpp
#define flopoco_random_utils_chain_operator_hpp

#include "Operator.hpp"

namespace flopoco
{
namespace random
{

class ChainOperator
	: public Operator
{
public:
	/* The mapping info defines what happens to the output ports of the first operator,
		and the input ports of the second operator (all inputs of the first operator are
		exposed, and all outputs of the second operator are exposed).
		(X,Y) - Connect output X of the first operator to input Y of the second operator
		(X,-) - Output X of the first operator is unconnected
		(X,+) - Output X of the first operator is brought out to top level
		(+,Y) - Input Y of the first operator is brought out to the top level

		It is currently an error to have an input of the second operator left dangling, it must
		be connected to a top-level input or an output of the first operator, so the default
		is to map to (+,Y) if it is not specified.

		If an output of the first operator is not mentioned, then it is equivalent to (X,-), i.e.
		it is not connected.
	*/
	typedef std::pair<std::string,std::string> mapping_t;
	typedef std::list<std::pair<std::string,std::string> > mapping_list_t;
	
private:	
	Operator *m_a;
	mapping_list_t m_mappings;
	Operator *m_b;

	// Get mapping for output of first operator
	mapping_t getOutputMapping(std::string name)
	{
		mapping_list_t::iterator it=m_mappings.begin();
		while(it!=m_mappings.end()){
			if(it->first==name)
				return *it;
			++it;
		}
		return mapping_t(name,"-");
	}
	
	// Get mapping for input of second operator
	mapping_t getInputMapping(std::string name)
	{
		mapping_list_t::iterator it=m_mappings.begin();
		while(it!=m_mappings.end()){
			if(it->second==name)
				return *it;
			++it;
		}
		return mapping_t("+",name);
	}


	ChainOperator(
		std::string name,
		Operator *a,
		const mapping_list_t &mapping,
		Operator *b
	)
		: Operator(a->getTarget())
		, m_a(a)
		, m_mappings(mapping)
		, m_b(b)
	{
		setName(name);
		
		oplist.push_back(a);
		oplist.push_back(b);
		
		// All inputs of first operator to first level
		for(int i=0;i<a->getNumberOfInputs();i++){
			Signal *sig=a->getInputSignal(i);
			addInput(sig->getName(), sig->width());
			inPortMap(a, sig->getName(), sig->getName());
		}
		
		// Outputs of first operator go to various places
		for(int i=0;i<a->getNumberOfOutputs();i++){
			Signal *sig=a->getOutputSignal(i);
			mapping_t m=getOutputMapping(sig->getName());
			if(m.second=="-"){
				// do nothing
			}else if(m.second=="+"){
				addOutput(sig->getName(), sig->width());
				outPortMap(a, sig->getName(), sig->getName()+"_out");
			}else{
				outPortMap(a, sig->getName(), m.first+"_to_"+m.second);
			}
		}
		vhdl << instance(a, "a");
		
		// Do some synchronisation - I'm not absolutely sure if this is needed
		for(int i=0;i<a->getNumberOfOutputs();i++){
			Signal *sig=a->getOutputSignal(i);
			mapping_t m=getOutputMapping(sig->getName());
			if(m.second=="+"){
				syncCycleFromSignal(sig->getName()+"_out");
			}else if(m.second!="-"){
				syncCycleFromSignal(m.first+"_to_"+m.second);
			}
		}
		
		// We need to get a out of the way before starting on b, as in principle a==b might happen
		
		// Output's second operator come from various places
		for(int i=0;i<b->getNumberOfInputs();i++){
			Signal *sig=b->getInputSignal(i);
			mapping_t m=getInputMapping(sig->getName());
			if(m.first=="-"){
				throw std::string("ChainOperator - inputs of second operator can't be unbound.");
			}else if(m.first=="+"){
				addInput(sig->getName(), sig->width());
				inPortMap(b, sig->getName(), sig->getName());
			}else{
				inPortMap(b, sig->getName(), m.first+"_to_"+m.second);
			}
		}
		// All outputs go to top level
		for(int i=0;i<b->getNumberOfOutputs();i++){
			Signal *sig=b->getOutputSignal(i);
			addOutput(sig->getName(), sig->width());
			outPortMap(b, sig->getName(), join(sig->getName(),"_out"));
		}
		vhdl<<instance(b, "b");
		
		// More synchronisation - again, not sure if this is needed
		for(int i=0;i<b->getNumberOfOutputs();i++){
			Signal *sig=b->getOutputSignal(i);
			syncCycleFromSignal(sig->getName()+"_out");
		}
		
		// Handle any outputs of the first operator going to the top level
		for(int i=0;i<a->getNumberOfOutputs();i++){
			Signal *sig=a->getOutputSignal(i);
			if(getOutputMapping(sig->getName()).second=="+"){
				vhdl<<tab<<sig->getName()<<" <= "<<sig->getName()<<"_out;\n";
			}
		}
		
		// Handle all outputs of the second operator (all go to top level)
		for(int i=0;i<b->getNumberOfOutputs();i++){
			Signal *sig=b->getOutputSignal(i);
			vhdl<<tab<<sig->getName()<<" <= "<<sig->getName()<<"_out;\n";
		}
	}
	
	TestCase *MapTestCase(TestCase *src)
	{
		if(m_a->getNumberOfInputs() < getNumberOfInputs())
			throw std::string("ChainOperator - can't generate test cases if inputs of second operator come to top level.");
		
		TestCase *res=new TestCase(this);
		res->setSetupCycle(src->getSetupCycle());
		
		for(int i=0;i<getNumberOfInputs();i++){
			Signal *sig=getInputSignal(i);
			res->setInputValue(sig->getName(), src->getInputValue(sig->getName()));
		}
		
		return res;
	}
public:
	static Operator *Create(std::string name,Operator *a,const mapping_list_t &mapping,Operator *b)
	{
		return new ChainOperator(name, a, mapping, b);
	}
	
	void emulate(TestCase *t)
	{
		if(m_a->getNumberOfInputs() < getNumberOfInputs())
			throw std::string("ChainOperator - can't emulate if inputs of second operator come to top level.");
		
		TestCase *ta=new TestCase(m_a);
		for(int i=0;i<m_a->getNumberOfInputs();i++){
			Signal *sig=m_a->getInputSignal(i);
			ta->setInputValue(sig->getName(), t->getInputValue(sig->getName()));
		}
		
		m_a->emulate(ta);
		
		TestCase *tb=new TestCase(m_b);
		
		for(int i=0;i<m_a->getNumberOfOutputs();i++){
			Signal *sig=m_a->getOutputSignal(i);
			mapping_t m=getOutputMapping(sig->getName());
			std::vector<mpz_class> ee=ta->getExpectedOutputValues(sig->getName());
			if(m.second=="-"){
				// do nothing
			}else if(m.second=="+"){
				for(unsigned j=0;j<ee.size();j++){
					t->addExpectedOutput(sig->getName(), ee[j]);
				}
			}else{
				if(ee.size()!=1)
					throw std::string("ChainOperator - Can only emulate with first operator that produces single expected value.");
				tb->setInputValue(m.second, ee[0]);
			}
		}
		
		delete ta;
		
		m_b->emulate(tb);
		
		for(int i=0;i<m_b->getNumberOfOutputs();i++){
			Signal *sig=m_b->getOutputSignal(i);
			std::vector<mpz_class> ee=tb->getExpectedOutputValues(sig->getName());
			for(unsigned j=0;j<ee.size();j++){
				t->addExpectedOutput(sig->getName(), ee[j]);
			}
		}
		
		delete tb;
	}
	
	void buildStandardTestCases(TestCaseList* tcl)
	{
		fprintf(stderr, "Here\n");
		
		TestCaseList *firstLevel=new TestCaseList();
		
		m_a->buildStandardTestCases(firstLevel);
		
		for(int i=0;i<firstLevel->getNumberOfTestCases();i++){
			
			// This will give us the inputs for the firstlevel
			TestCase *src=firstLevel->getTestCase(i);
			
			TestCase *actual=MapTestCase(src);
			emulate(actual);
			tcl->add(actual);
		}
		
		fprintf(stderr, "There\n");
		
		delete firstLevel;
	}
};

}; // random
}; // flopoco

#endif
