// general c++ library for manipulating streams
#include <iostream>
#include <sstream>

/* header of libraries to manipulate multiprecision numbers
   There will be used in the emulate function to manipulate arbitraly large
   entries */
#include "gmp.h"
#include "mpfr.h"

// include the header of the Operator
#include "ConfigurationController.hpp"

using namespace std;
namespace flopoco
{




    ConfigurationController::ConfigurationController(Target* target) : Operator(target)
    {
        srcFileName="ConfigurationController";
        this->useNumericStd();
		// definition of the name of the operator
		ostringstream name;
        name << "ConfigurationController";
        setName(name.str());
        // Copyright
        setCopyrightString("UNIVERSITY of Kassel 2017");


				addInput("init_conv_i");
				addOutput("ce_o");
				addOutput("ctrl_int",3,true);
				
				addType("states","(idle,count)");

				addType("rom","array ( 0 to 31) of std_logic_vector(2 downto 0)");
                //the order of the two 4 input LUT configurations has to switch because of the exchange of O5 and O6 (enforced by a mapping Error)
                //addConstant("ctrl_lut","rom","(\"001\",\"000\",\"000\",\"000\",\"000\",\"000\",\"000\",\"000\",\"010\",\"000\",\"000\",\"000\",\"000\",\"000\",\"000\",\"000\",\"101\",\"100\",\"100\",\"100\",\"100\",\"100\",\"100\",\"100\",\"110\",\"100\",\"100\",\"100\",\"100\",\"100\",\"100\",\"100\");\nsignal state : states");
                addConstant("ctrl_lut","rom","(\"101\",\"100\",\"100\",\"100\",\"100\",\"100\",\"100\",\"100\",\"110\",\"100\",\"100\",\"100\",\"100\",\"100\",\"100\",\"100\",\"001\",\"000\",\"000\",\"000\",\"000\",\"000\",\"000\",\"000\",\"010\",\"000\",\"000\",\"000\",\"000\",\"000\",\"000\",\"000\");\nsignal state : states");
   // Attention, the line above contains a really dirty hack to use an own type for the state machine, which should be fix in Operator and Signal
				vhdl << "coeff_switch: process(clk,init_conv_i) -- state machine to replace coefficients  " << endl;
				vhdl << "  variable counter: integer range 0 to 31;                                         " << endl;
				vhdl << "begin                                                                              " << endl;
				vhdl << "  if clk'event and clk = '1' then                                              " << endl;
				vhdl << "      case state is                                                                " << endl;
				vhdl << "      when idle =>                                                                 " << endl;
				vhdl << "        if init_conv_i = '1' then                                                  " << endl;
				vhdl << "          state <= count;                                                          " << endl;
				vhdl << "          counter := 0;                                                            " << endl;
				vhdl << "          ce_o <= '1';                                                             " << endl;
				vhdl << "        end if;                                                                    " << endl;
				vhdl << "      when count =>                                                                " << endl;
				vhdl << "        if counter = 31 then                                                       " << endl;
				vhdl << "          state <= idle;                                                           " << endl;
				vhdl << "          ce_o <= '0';                                                             " << endl;
				vhdl << "        else                                                                       " << endl;
				vhdl << "          counter := counter + 1;                                                  " << endl;
				vhdl << "          state <= count;                                                          " << endl;
				vhdl << "          ce_o <= '1';                                                             " << endl;
				vhdl << "        end if;                                                                    " << endl;
				vhdl << "      end case;                                                                    " << endl;
				vhdl << "  end if;                                                                          " << endl;
				vhdl << "  ctrl_int <= ctrl_lut(counter);                                                   " << endl;
				vhdl << "                                                                                   " << endl;
				vhdl << "end process;                                                                       " << endl;
    }

    void ConfigurationController::emulate(TestCase * tc)
    {
	}

    void ConfigurationController::buildStandardTestCases(TestCaseList * tcl)
    {
		// please fill me with regression tests or corner case tests!
	}

    OperatorPtr ConfigurationController::parseArguments(Target *target, vector<string> &args)
    {
        return new ConfigurationController(target);
	}
	
    void ConfigurationController::registerFactory()
    {
        UserInterface::add("ConfigurationController", // name
                       "ConfigurationController for NN reconfiguration in 32 clock cycles", // description, string
											 "NeuralNetworks", // category, from the list defined in UserInterface.cpp
											 "", //seeAlso
											 // Now comes the parameter description string.
											 // Respect its syntax because it will be used to generate the parser and the docs
											 // Syntax is: a semicolon-separated list of parameterDescription;
											 // where parameterDescription is parameterName (parameterType)[=defaultValue]: parameterDescriptionString 
                      "",
                      "",                  
											 // More documentation for the HTML pages. If you want to link to your blog, it is here.
                      ConfigurationController::parseArguments
											 ) ;
	}

}//namespace
