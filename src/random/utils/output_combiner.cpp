#include "Operator.hpp"

#include "random/utils/operator_factory.hpp"

namespace flopoco
{
namespace random
{

/*! This operator takes every single output of the operator and xors it into a giant LFSR-type
	structure, which the tools can't optimise out, and which takes precisely one LUT-FF per output bit.
*/
class OutputCombinerOperator
	: public Operator
{
	OutputCombinerOperator(
		Operator *op
	)
		: Operator(op->getTarget())
	{
		setName(op->getName()+"combiner");
		
		// We have recurrences of one cycle, anything more is an error
		setHasDelay1Feedbacks();
		
		oplist.push_back(op);
		
		std::map<std::string,std::pair<int,int> > mappings;
		
		for(int i=0;i<op->getNumberOfInputs();i++){
			Signal *sig=op->getInputSignal(i);
			addInput(sig->getName(), sig->width(), sig->width()>1);
			
			inPortMap(op, sig->getName(), sig->getName());
		}
		for(int i=0;i<op->getNumberOfOutputs();i++){
			Signal *sig=op->getOutputSignal(i);
			outPortMap(op, sig->getName(), sig->getName());
		}
		
		vhdl << instance(op, "op");
		for(int i=0;i<op->getNumberOfOutputs();i++){
			Signal *sig=op->getOutputSignal(i);
			syncCycleFromSignal(sig->getName());
		}
		
		// Ok, now we are synced up to every output	
		
		int total=0;
		for(int i=0;i<op->getNumberOfOutputs();i++){
			Signal *sig=op->getOutputSignal(i);
			for(int j=0;j<sig->width();j++){
				declare(join("reg_",total),1,false, Signal::registeredWithZeroInitialiser);
				vhdl << tab << use(join("reg_",total)) << "<=" << join("val_",total) << ";" << endl;
				total++;
			}
		}
		vhdl<<declare(join("reg_",total))<<" <= '0';\n";
		
		nextCycle();
		
		total=0;
		for(int i=0;i<op->getNumberOfOutputs();i++){
			Signal *sig=op->getOutputSignal(i);
			if(sig->width()==1){
				vhdl << declare(join("val_",total))<<" <= "<< sig->getName() << " xor "<<join("reg_",total)<<" xor "<<join("reg_",total+1)<<";\n";
				total++;
			}else{
				for(int j=0;j<sig->width();j++){
					vhdl << declare(join("val_",total))<<" <= "<< sig->getName() << "("<<j<<") xor "<<join("reg_",total)<<" xor "<<join("reg_",total+1)<<";\n";
					total++;
				}
			}
		}
		
		addOutput("O");
		vhdl << " O <= reg_0;\n";
	}
public:
	static OutputCombinerOperator *Create(Operator *op)
	{
		return new OutputCombinerOperator(op);
	}
};


static void OutputCombinerUsage(std::ostream &dst)
{
	OperatorFactory::classic_OP(dst, "OutputCombiner", "<no arguments>", false);
	dst << "    Take whatever was the last operator, and combine all n bits of its output\n";
	dst << "  into a single bit, using n LUT-FF pairs.\n";
	dst << "\n";
}



static Operator *OutputCombinerParser(Target *target ,const std::vector<std::string> &args,int &consumed)
{
	consumed=0;
	
	std::vector<Operator*> * ops=target->getGlobalOpListRef();
	if(ops->size()==0)
		throw std::string("OutputCombinerParser - No operator to combine outputs of.");
	
	Operator *op=ops->back();
	return OutputCombinerOperator::Create(op);
}

void OutputCombiner_registerFactory()
{
	DefaultOperatorFactory::Register(
		"OutputCombiner",
		"operator",
		OutputCombinerUsage,
		OutputCombinerParser
	);
}

}; // random
}; // flopoco
