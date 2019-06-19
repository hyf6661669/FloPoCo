#ifndef  CONSTMULTPAG_HPP
#define CONSTMULTPAG_HPP

/*  CONSTMULTPAG.hpp
 *  Version:
 *  Datum: 20.11.2014
 */
#ifdef HAVE_PAGLIB

#include "Operator.hpp"
#include "utils.hpp"
#include "pagsuite/adder_graph.h"
#include "ConstMultPAG_types.hpp"

using namespace std;

namespace flopoco {
  class ConstMultPAG : public Operator {
    public:
     static ostream nostream;
     int noOfPipelineStages;

     ConstMultPAG(Target* target,int wIn_, string pipelined_realization_str, bool pipelined_=true, bool syncInOut_=true, int syncEveryN_=1,bool syncMux_=true);
     ConstMultPAG(Target* target, int wIn_, vector<vector<int64_t> > &coefficients, bool pipelined_=true, bool syncInOut_=true, int syncEveryN_=1, bool syncMux_=true);

      ~ConstMultPAG() {}

     void emulate(TestCase * tc);
     void buildStandardTestCases(TestCaseList* tcl);
     struct output_signal_info{
         string signal_name;
         vector<vector<int64_t> > output_factors;
         int wordsize;};
     list<output_signal_info>& GetOutputList();
     static OperatorPtr parseArguments( Target *target, vector<string> &args );

     static void registerFactory();

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


 };
}//namespace
#endif // HAVE_PAGLIB
#endif
