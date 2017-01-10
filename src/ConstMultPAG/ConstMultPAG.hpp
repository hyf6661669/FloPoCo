#ifndef  CONSTMULTPAG_HPP
#define CONSTMULTPAG_HPP

/*  CONSTMULTPAG.hpp
 *  Version:
 *  Datum: 20.11.2014
 */

#ifdef HAVE_PAGLIB

#include "Operator.hpp"
#include "utils.hpp"
#include "pag_lib/adder_graph.h"
#include "ConstMultPAG_types.hpp"

using namespace std;

namespace flopoco {
  class ConstMultPAG : public Operator {
    public:
     static ostream nostream;
     int noOfPipelineStages;

     ConstMultPAG(Target* target,int wIn_, char* pipelined_realization_str, bool pipelined_=true, bool syncInOut_=true, int syncEveryN_=1,bool syncMux_=true);
     ConstMultPAG(Target* target, int wIn_, vector<vector<int64_t> > &coefficients, bool pipelined_=true, bool syncInOut_=true, int syncEveryN_=1, bool syncMux_=true);

      ~ConstMultPAG() {}

     void emulate(TestCase * tc);
     void buildStandardTestCases(TestCaseList* tcl);
     struct output_signal_info{
         string signal_name;
         vector<vector<int64_t> > output_factors;
         int wordsize;};
     list<output_signal_info>& GetOutputList();

  protected:
     int wIn;
     bool syncInOut;
     bool pipelined;
     int syncEveryN;
     bool syncMux;

     bool RPAGused;
     int emu_conf;
     vector<vector<int64_t> > output_coefficients;
     adder_graph_t pipelined_adder_graph;
     list<string> input_signals;

     list<output_signal_info> output_signals;

     int noOfInputs;
     int noOfConfigurations;
     bool needs_unisim;

     void ProcessConstMultPAG(Target* target, string pipelined_realization_str);


     string generateSignalName(adder_graph_base_node_t* node);
     ConstMultPAG_TYPES::ConstMultPAG_BASE* identifyNodeType(adder_graph_base_node_t* node);
     bool TryRunRPAG(string pipelined_realization_str,string& out);

     void identifyOutputConnections(adder_graph_base_node_t* node, map<adder_graph_base_node_t *, ConstMultPAG_TYPES::ConstMultPAG_BASE *> &infoMap);
     void printAdditionalNodeInfo(map<adder_graph_base_node_t *, ConstMultPAG_TYPES::ConstMultPAG_BASE *> &infoMap );
     string getShiftAndResizeString(string signalName, int outputWordsize, int inputShift, bool signedConversion=true);
     string getBinary(int value, int wordsize);


	OperatorPtr parseArguments( Target *target, vector<string> &args ) {
        int wIn, sync_every = 0;
        std::string graph;
		bool pipeline,sync_inout,sync_muxes;

        UserInterface::parseInt( args, "wIn", &wIn );
        UserInterface::parseString( args, "graph", &graph );
		UserInterface::parseBool( args, "pipeline", &pipeline );
		UserInterface::parseBool( args, "sync_inout", &sync_inout );
		UserInterface::parseBool( args, "sync_muxes", &sync_muxes );
        UserInterface::parseInt( args, "sync_every", &sync_every );

        return new ConstMultPAG( target, wIn, graph.c_str(), pipeline,sync_inout,sync_every,sync_mux );
    }

    void registerFactory() {
        UserInterface::add( "ConstMultPAG", // name
                            "A component for building constant multipliers based on pipelined adder graphs(pag).", // description, string
                            "BasicInteger", // category, from the list defined in UserInterface.cpp
                            "", //seeAlso
                            "wIn (int): Wordsize of pag inputs; \
                            graph (string): Realization string of the pag; \
                            pipeline (bool)=true: Enable pipelining of the pag; \
                            sync_inout (bool)=true: Enable pipeline registers for input and output stage; \
							sync_muxes (bool)=true: Enable counting mux-only stages as full stage; \
                            sync_every (int)=1: Count of stages after which will be pipelined",
                            "Nope.",
                            ConstMultPAG::parseArguments
                          ) ;
  };


}//namespace
#endif // HAVE_PAGLIB
#endif
