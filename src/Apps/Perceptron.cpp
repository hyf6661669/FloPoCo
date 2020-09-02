#include <iostream>
#include <sstream>
#include "Perceptron.hpp"
#include "BitHeap/BitHeap.hpp"
#include "FixFunctions/FixFunctionByTable.hpp"

using namespace std;
namespace flopoco {
    Perceptron::Perceptron(OperatorPtr parentOp, Target* target, int inputMsb_, int inputLsb_, int expLsb_, int sumMsb_,
        int sumLsb_) :
    Operator(parentOp, target), inputMsb(inputMsb_), inputLsb(inputLsb_), expLsb(expLsb_), sumMsb(sumMsb_), sumLsb(sumLsb_)
    {
        srcFileName="Perceptron";

        ostringstream name;
        name << "Perceptron_" << inputMsb << "_" << (inputLsb < 0 ? "N":"") << abs(inputLsb);
        setName(name.str());
        setCopyrightString("Maxime Christ, Florent de Dinechin 2020");

        // IO defintions
        // Inputs X are unsigned
        // Inputs W have 1 sign bit + coded in 2's complement
        int inputPrec = inputMsb - inputLsb + 1;
        int prevLayerWidth = 10;
        // int layerWidth = 10;
        for (int i = 0; i < prevLayerWidth; i++) {
            string inputSignal = join("X", i);
            string weightSignal = join("W", i);
            addInput(inputSignal, inputMsb - inputLsb + 1, true); // zero flag + unsigned (always positive)
            addInput(weightSignal, inputMsb - inputLsb + 2, true); // zero flag + sign bit + signed int

            vhdl << tab << declare(join("zeroW", i)) << " <= " << weightSignal << "(" << inputPrec << ");" << endl;
            vhdl << tab << declare(join("signW", i)) << " <= " << weightSignal << "(" << inputPrec - 1 << ");" << endl;
            vhdl << tab << declareFixPoint(join("valueW", i), true, inputMsb - 2, inputLsb) << " <= signed(" << weightSignal;
            vhdl<< "(" << inputPrec - 2 << " downto 0));" << endl;

            vhdl << tab << declare(join("zeroX", i)) << " <= " << inputSignal << "(" << inputPrec << ");" << endl;
            vhdl << tab << declareFixPoint(join("valueX", i), true, inputMsb - 1, inputLsb) << " <= signed(" << inputSignal;
            vhdl<< "(" << inputPrec - 1 << " downto 0));" << endl;
        }

        addOutput("A", inputPrec);

        // LNS adder
        for (int i = 0; i < prevLayerWidth; i++) {
            vhdl << tab << declare(join("z", i)) << " <= " << join("zeroW", i) << " or " << join("zeroX", i) << ";" << endl;
            vhdl << tab << declareFixPoint(join("LogProduct", i), true, inputMsb + 1, inputLsb) << " <= ";
            vhdl << join("valueX", i) << " + " << join("valueW", i) << ";" << endl;
        }

        // 2^x fction
        REPORT(0, "inputLSB is: " << inputLsb << " and expLsb is: " << expLsb);
        REPORT(0, "EXP VECTOR IS: ");
        vector<mpz_class> exp_vector;
        for(int i = 0; i < (1<<-inputLsb + 1); i++) {
            int s = i >> -inputLsb;
            int x = i - (s << -inputLsb);

            int input_scale = 1 << -inputLsb;
            double xd = (double)x / input_scale; // exact

            double rd;
            if (s == 0)
                rd = exp2(xd); // rounding @ lsb=-53
            else
                rd = -exp2(xd);

            int output_scale = 1 << -expLsb;
            int r = (int)(rd * output_scale); // exact
            mpz_class exp_x = r;
            exp_x = signedToBitVector(exp_x, -expLsb + 2);
            exp_vector.push_back(exp_x);
            REPORT(0, "; X = " << x << " exp = " << exp_x);
        }
        REPORT(0, endl);

        for (int i = 0; i < prevLayerWidth; i++) {
            string expInput = join("expInput", i);
            string expOutput = join("expOutput", i);
            declare(expInput, -inputLsb + 2);
            declare(expOutput, -expLsb + 2);
            vhdl << tab << expInput << " <= std_logic_vector(" << join("signW", i);
            vhdl << " & " << join("LogProduct", i) << "(" << -inputLsb << "downto 0));" << endl;

            Table::newUniqueInstance(this, expInput, expOutput, exp_vector,
                                        join("ExpTable", i), -inputLsb + 2,
                                        -expLsb + 2);

            vhdl << tab << declareFixPoint(join("summand", i), true, 1, expLsb) << " <= signed(";
            vhdl << zg(2 - expLsb) << " when " << join("z", i) << " = '1' eles "<< expOutput << " when " << join("signW", i) << " = '0' else not " << expOutput << ");";
            // vhdl << rangeAssign(1 - expLsb, 0,join("signW", i)) << ";" << endl;
        }

        // Accumulator
        BitHeap* bh = new BitHeap(this, sumMsb, sumLsb);

        for (int i = 0; i < prevLayerWidth; i++)
            bh->addSignal(join("summand", i), 0);
        bh->startCompression();


        vhdl << tab << "A <= " << bh->getSumName(1, -4) << ";" << endl;
        // Trunc, act + log
    };

    OperatorPtr Perceptron::parseArguments(OperatorPtr parentOp, Target* target, vector<string> &args) {
        int inputMsb_, inputLsb_, expLsb_, sumMsb_, sumLsb_;
        UserInterface::parseInt(args, "inputMsb", &inputMsb_);
        UserInterface::parseInt(args, "inputLsb", &inputLsb_);
        UserInterface::parseInt(args, "expLsb", &expLsb_);
        UserInterface::parseInt(args, "sumMsb", &sumMsb_);
        UserInterface::parseInt(args, "sumLsb", &sumLsb_);
        return new Perceptron(parentOp, target, inputMsb_, inputLsb_, expLsb_, sumMsb_, sumLsb_);
    };

    void Perceptron::registerFactory() {
        UserInterface::add("Perceptron",
                           "Emulate a neuron using LNS to avoid FP multipliers",
                           "Primitives",
                           "",
                           "inputMsb(int): input exponent size; \
                           inputLsb(int): input fractional size; \
                           expLsb(int): exponential func output lsb; \
                           sumMsb(int): final summation msb; \
                           sumLsb(int): final summation lsb;",
                           "",
                           Perceptron::parseArguments);
    };

    void Perceptron::emulate(TestCase * tc) {};
    void Perceptron::buildStandardTestCases(TestCaseList* tcl) {};

}
