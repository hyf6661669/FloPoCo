#include <iostream>
#include <sstream>

#include "gmp.h"
#include "mpfr.h"
#include "GenericAddSub.hpp"
#include "Primitive.hpp"

using namespace std;
namespace flopoco {
    GenericAddSub::GenericAddSub(Target* target, const uint32_t &wIn, const uint32_t &flags) : Operator(target), flags_(flags) {
        setCopyrightString( UniKs::getAuthorsString( UniKs::AUTHOR_MKLEINLEIN ) );
        this->useNumericStd();
        srcFileName="GenericAddSub";
		ostringstream name;
        name << "GenericAddSub_w" << wIn << "_" << printFlags();
        setName(name.str());

        addInput("iL",wIn);
        addInput("iR",wIn);
        if(hasFlags(TERNARY))
            addInput("IM",wIn);

        if( hasFlags(CONF_LEFT) )
            addInput("iL_c",1,false);
        if( hasFlags(CONF_RIGHT) )
            addInput("iR_c",1,false);
        if( hasFlags(TERNARY|CONF_MID) )
            addInput("iM_c",1,false);

        addOutput("sum_o",wIn);

        if( UserInterface::useTargetSpecificOptimization && target->getVendor() == "Xilinx" )
            buildXilinx(target,wIn);
        else if( UserInterface::useTargetSpecificOptimization && target->getVendor()=="Altera" )
            buildAltera(target,wIn);
        else
            buildCommon(target,wIn);

    }

    void GenericAddSub::buildXilinx(Target* target, const uint32_t &wIn){
        REPORT(LIST,"Xilinx junction not fully implemented, fall back to common.");
        buildCommon(target,wIn);

    }

    void GenericAddSub::buildAltera(Target *target, const uint32_t &wIn){
        REPORT(LIST,"Altera junction not fully implemented, fall back to common.");
        buildCommon(target,wIn);

    }

    void GenericAddSub::buildCommon(Target* target, const uint32_t &wIn){
        const uint16_t c_count = (hasFlags(CONF_LEFT)?1:0) + (hasFlags(CONF_RIGHT)?1:0) + (hasFlags(TERNARY|CONF_MID)?1:0);
        if(c_count >0){
            vhdl << declare( "CONF", c_count ) << " <= ";
            if( hasFlags(CONF_LEFT) ){
                vhdl << "iL_c";
                if( hasFlags(TERNARY|CONF_MID) ){
                    vhdl << "& iM_c";
                }
                if( hasFlags(CONF_RIGHT) ){
                    vhdl << "& iR_c";
                }
            }else if( hasFlags(TERNARY|CONF_MID) ){
                vhdl << "iM_c";
                if( hasFlags(CONF_RIGHT) ){
                    vhdl << "& iM_c";
                }
            }else if( hasFlags(CONF_RIGHT) ){
                vhdl << "& iM_c";
            }
            vhdl << std::endl;
            vhdl << "case CONF is" << std::endl;
            for( int i = 0;i<(1<<c_count);++i ){
                vhdl << "\t" << "when \"";
                for( int j=c_count-1;j>=0;--j )
                    vhdl << (i&(1<<j)?"1":"0");
                vhdl << "\"\t=> std_logic_vector(";
                uint16_t mask = 1<<(c_count-1);
                if( hasFlags(CONF_LEFT) ){
                    vhdl << (i&(mask)?"-":"+") << "unsigned(iL)";
                    mask >>= 1;
                    if( hasFlags(TERNARY|CONF_MID) ){
                        vhdl << (i&mask?"-":"+") << "unsigned(iM)";
                        mask >>= 1;
                    }
                    if( hasFlags(CONF_RIGHT) ){
                        vhdl << (i&mask?"-":"+") << "unsigned(iR)";
                    }
                }else if( hasFlags(TERNARY|CONF_MID) ){
                    vhdl << "unsigned(iL)";
                    vhdl << (i&mask?"-":"+") << "unsigned(iM)";
                    mask >>= 1;
                    if( hasFlags(CONF_RIGHT) ){
                        vhdl << (i&mask?"-":"+") << "unsigned(iR)";
                    }
                }else if( hasFlags(CONF_RIGHT) ){
                    vhdl << (i&mask?"-":"+") << "unsigned(iR)";
                }
                vhdl << ");" << std::endl;
            }
            vhdl << "\t" << "when others => sum_o <= (others=>'X');" << std::endl;
            vhdl << "end case;" << std::endl;
        }
        else
        {
            vhdl << "\tsum_o <= std_logic_vector(" << (hasFlags(SUB_LEFT)?"-":"");
            if(hasFlags(TERNARY))
            {
                vhdl << "signed(iM)" << (hasFlags(SUB_MID)?"-":"+") << " signed(iM)";
            }
            vhdl << "signed(iL)" << (hasFlags(SUB_RIGHT)?"-":"+") << " signed(iR));" << endl;
        }
    }

    string GenericAddSub::getInputName(const uint32_t &index, const bool &c_input) const{
        switch(index){
        case 0: return "iL" + std::string(c_input?"_c":"");
        case 1: return std::string(hasFlags(TERNARY|CONF_MID)?"iM":"iR") + std::string(c_input?"_c":"");
        case 2: return std::string(hasFlags(TERNARY|CONF_MID)?"iR" + std::string(c_input?"_c":""):"");
        default: return "";
        }
    }

    string GenericAddSub::getOutputName() const{
        return "sum_o";
    }

    bool GenericAddSub::hasFlags(const uint32_t &flag) const{
        return flags_ & flag;
    }

    const uint32_t GenericAddSub::getInputCount() const{
        uint32_t c = (hasFlags(TERNARY)?3:2);
        if(hasFlags(CONF_LEFT)) c++;
        if(hasFlags(TERNARY|CONF_MID)) c++;
        if(hasFlags(CONF_RIGHT)) c++;
        return c;
    }

    string GenericAddSub::printFlags() const{
        std::stringstream o;
        o << (hasFlags(SUB_LEFT)?"s":"");
        o << (hasFlags(CONF_LEFT)?"c":"");
        o << "L";

        if( hasFlags(TERNARY) ){
            o << (hasFlags(SUB_MID)?"s":"");
            o << (hasFlags(CONF_MID)?"c":"");
            o << "M";
        }

        o << (hasFlags(SUB_RIGHT)?"s":"");
        o << (hasFlags(CONF_RIGHT)?"c":"");
        o << "R";
        return o.str();
    }


    void GenericAddSub::emulate(TestCase * tc) {

	}


    void GenericAddSub::buildStandardTestCases(TestCaseList * tcl) {
		// please fill me with regression tests or corner case tests!
	}




#if 0
    OperatorPtr GenericAddSub::parseArguments(Target *target, vector<string> &args) {
        /*int param0, param1;
		UserInterface::parseInt(args, "param0", &param0); // param0 has a default value, this method will recover it if it doesnt't find it in args, 
		UserInterface::parseInt(args, "param1", &param1); 
        return new UserDefinedOperator(target, param0, param1);*/
        return nullptr;
	}
	
    void GenericAddSub::registerFactory(){
        /*UserInterface::add("UserDefinedOperator", // name
											 "An heavily commented example operator to start with FloPoCo.", // description, string
											 "Miscellaneous", // category, from the list defined in UserInterface.cpp
											 "", //seeAlso
											 // Now comes the parameter description string.
											 // Respect its syntax because it will be used to generate the parser and the docs
											 // Syntax is: a semicolon-separated list of parameterDescription;
											 // where parameterDescription is parameterName (parameterType)[=defaultValue]: parameterDescriptionString 
											 "param0(int)=16: A first parameter, here used as the input size; \
                        param1(int): A second parameter, here used as the output size",
											 // More documentation for the HTML pages. If you want to link to your blog, it is here.
											 "Feel free to experiment with its code, it will not break anything in FloPoCo. <br> Also see the developper manual in the doc/ directory of FloPoCo.",
											 UserDefinedOperator::parseArguments
                                             ) ;*/
	}
#endif
}//namespace
