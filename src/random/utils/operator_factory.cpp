#include "operator_factory.hpp"

void addOperator(flopoco::Operator *op);

namespace flopoco
{	
namespace random
{

std::vector<OperatorFactoryPtr> OperatorFactory::sm_factoriesByIndex;
std::map<std::string,OperatorFactoryPtr> OperatorFactory:: sm_factoriesByName;

void OperatorFactory::RegisterOperatorFactory(OperatorFactoryPtr factory)
{
	if(sm_factoriesByName.find(factory->Name())!=sm_factoriesByName.end())
		throw std::string("OperatorFactory - Factory with name '"+factory->Name()+" has already been registered.");
	
	sm_factoriesByIndex.push_back(factory);
	sm_factoriesByName.insert(std::make_pair(factory->Name(),factory));
}

unsigned OperatorFactory::GetFactoryCount()
{
	return sm_factoriesByIndex.size();
}

OperatorFactoryPtr OperatorFactory::GetFactoryByIndex(unsigned i)
{
	return sm_factoriesByIndex.at(i);
}

OperatorFactoryPtr OperatorFactory::FindFactory(std::string operatorName)
{
	return sm_factoriesByName[operatorName];
}

void OperatorFactory::classic_usage(char *name, string opName)
{
	std::cerr<<"NumFactories="<<sm_factoriesByIndex.size()<<"\n";
	if(opName==""){
		for(unsigned i=0;i<sm_factoriesByIndex.size();i++){
			std::cerr<<"Name = "<<sm_factoriesByIndex[i]->Name()<<"\n";
			sm_factoriesByIndex[i]->Usage(std::cerr);
		}
		
	}else{
		OperatorFactoryPtr factory=sm_factoriesByName[opName];
		if(factory)
			factory->Usage(std::cerr);
	}
}

bool OperatorFactory::classic_parseCommandLine(int argc, char* argv[], Target *target, std::string opname, int &i)
{
	OperatorFactoryPtr factory=sm_factoriesByName[opname];
	if(!factory)
		return false;
	
	// Collect only the arguments that are left
	std::vector<std::string> args;
	for(int ii=i;ii<argc;ii++){
		args.push_back(std::string(argv[ii]));
	}
	
	try{
		int consumed=0;
		Operator *res=factory->Parse(target, args, consumed);
		if(res!=NULL)	// Some factories don't actually create an operator
			addOperator(res);
		assert((consumed>=0) && (consumed<=(int)args.size()));
		i+=consumed;
		
		return true;
	}catch(std::string &s){
		std::cerr<<"Error : "<<s<<"\n";
		factory->Usage(std::cerr);
		exit(1);
	}catch(std::exception &s){
		std::cerr<<"Exception : "<<s.what()<<"\n";
		factory->Usage(std::cerr);
		exit(1);	
	}
}

void OperatorFactory::classic_OP(std::ostream &dst, std::string name, std::string args, bool newOperator)
{
	//static const int BRIGHT=1;
	//static const int RED=31;
	static const int OPER=32;
	static const int NEWOPER=32;
	static const int PARAM=34;
	
	int OP=newOperator?NEWOPER:OPER;
	
	dst << "    ";
		dst<<"\033[1;"<<OP<<"m";
	dst<<name;
		dst<<"\033[0m";
	dst<< " ";
		dst<<"\033[1;"<<PARAM<<"m";
	dst << args;
		dst<<"\033[0m";
	dst<<"\n";
}

DefaultOperatorFactory::SimpleOperatorFactory::SimpleOperatorFactory(
		std::string name,			
		std::string categories,
		usage_func_t usage,
		parser_func_t parser,
		const std::vector<std::vector<std::string> > &testParameters=std::vector<std::vector<std::string> >()
)
	: m_name(name)
	, m_usage(usage)
	, m_parser(parser)
	, m_testParams(testParameters)
{
	int start=0;
	while(start<(int)categories.size()){
		int end=categories.find(';', start);
		std::string part;
		if(end==-1)
			part=categories.substr(start, end);
		else
			part=categories.substr(start, end-start);
		if(part.size()!=0)
			m_categories.push_back(part);
		if(end==-1)
			break;
		start=end+1;
	}
}
		
void DefaultOperatorFactory::Register(
	std::string name,			
	std::string categories,	// semi-colon seperated list of categories
	usage_func_t usage,
	parser_func_t parser,
	const std::vector<std::vector<std::string> > &testParameters
){
	OperatorFactoryPtr factory(new SimpleOperatorFactory(name, categories, usage, parser, testParameters));
	OperatorFactory::RegisterOperatorFactory(factory);
}


}; // random
}; // flopoco
