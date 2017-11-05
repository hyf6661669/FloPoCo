// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "CoreTE.hpp"
#include "../ShiftReg.hpp"
#include "ConvolutionalCore.hpp"
#include "SOP_KCM.hpp"

#include <math.h>

using namespace std;
namespace flopoco {


CoreTE::CoreTE(Target *target, int input_bit_width_, int const_bit_width_, int no_of_products_) : Operator(target), input_bit_width(input_bit_width_), const_bit_width(const_bit_width_),no_of_products(no_of_products_)
{
    srcFileName="CoreTE";
    // definition of the name of the operator
    ostringstream name;
    name << "CoreTE_" << input_bit_width << "_" << const_bit_width << "_" << no_of_products;
    setName(name.str());
    // Copyright
    setCopyrightString("UNIVERSITY of Kassel 2017");
    int method= 1; // 1=SOPKCM, 2 = ConvolutionalCore

    addInput("X",input_bit_width);
    addInput("LUT_Config_clk_enable");



    int LUT_bit_width=4;
    int LUT_per_stage = const_bit_width+LUT_bit_width;
    int No_of_stages = ceil((float)const_bit_width / (float)LUT_bit_width);
    {
        int max_delay=no_of_products;
        int delay_bit_width = ceil(((float)LUT_per_stage)/2)*No_of_stages;

        ShiftReg* my_ShiftReg = new ShiftReg(target,delay_bit_width,max_delay);
        addToGlobalOpList(my_ShiftReg);
        addInput("cdi", delay_bit_width);
        inPortMap(my_ShiftReg,"X","cdi");
        for(int i=0; i< max_delay;++i)
        {
            //declare(join("Xd",i),delay_bit_width);
            outPortMap(my_ShiftReg,join("Xd",i),join("cdi_no_",i));
        }
        vhdl << instance(my_ShiftReg,"ShiftReg_cdi_inst")<< std::endl;
        //vhdl << my_ShiftReg->instance(this,"ShiftReg_cdi_inst") << std::endl;
    }

    {
        int max_delay=no_of_products;
        ShiftReg* my_ShiftReg = new ShiftReg(target,input_bit_width,max_delay);
        addToGlobalOpList(my_ShiftReg);
        inPortMap(my_ShiftReg,"X","X");
        for(int i=0; i< max_delay;++i)
        {
            //declare(join("Xd",i),delay_bit_width);
            outPortMap(my_ShiftReg,join("Xd",i),join("Xd", i));
        }


        vhdl << instance(my_ShiftReg,"ShiftReg_X_inst") << std::endl;
        //vhdl << my_ShiftReg->instance(this,"ShiftReg_X_inst") << std::endl;
    }

    SOP_KCM *my_SOP_KCM = new SOP_KCM(target,input_bit_width, const_bit_width,no_of_products,1);
    addToGlobalOpList(my_SOP_KCM);
    inPortMap(my_SOP_KCM,"LUT_Config_clk_enable","LUT_Config_clk_enable");
    for(int i=0; i < no_of_products; ++i)
    {
        inPortMap(my_SOP_KCM,join("X",i),join("Xd", i));
    }
    for(int i=0; i < no_of_products; ++i)
    {
        inPortMap(my_SOP_KCM,join("cdi_no_",i),join("cdi_no_",i));
    }

    outPortMap(my_SOP_KCM, "result","result");

    addOutput("Y",input_bit_width);
    vhdl <<  "Y <= result;" << std::endl;

    //input_No_0_cdi
    //LUT_Config_clk_enable



    vhdl << instance(my_SOP_KCM, "Core_inst") << std::endl;
}


OperatorPtr CoreTE::parseArguments(Target *target, vector<string> &args)
{
    int param0, param1, param2;


    UserInterface::parseInt(args, "input_bit_width", &param0); // param0 has a default value, this method will recover it if it doesnt't find it in args,
    UserInterface::parseInt(args, "const_bit_width", &param1);
    UserInterface::parseInt(args, "no_of_products", &param2);
    if(param1 == -1)
    {
        param1 = param0;
    }


    return new CoreTE(target, param0, param1, param2);
}

void CoreTE::registerFactory()
{
    UserInterface::add("CoreTE", // name
                                         "My first CoreTE.", // description, string
                                         "NeuralNetworks", // category, from the list defined in UserInterface.cpp
                                         "", //seeAlso
                                         // Now comes the parameter description string.
                                         // Respect its syntax because it will be used to generate the parser and the docs
                                         // Syntax is: a semicolon-separated list of parameterDescription;
                                         // where parameterDescription is parameterName (parameterType)[=defaultValue]: parameterDescriptionString
                                         "input_bit_width(int)=16: input word size; \
                                         const_bit_width(int)=-1: coefficient word size, per default the same as the input bit width;\
                                         no_of_products(int)=1: the nomber of products to accumulate",
                                         // More documentation for the HTML pages. If you want to link to your blog, it is here.
                                         "Feel free to experiment with its code, it will not break anything in FloPoCo. <br> Also see the developper manual in the doc/ directory of FloPoCo.",
                                         CoreTE::parseArguments
                                         ) ;

}

};
