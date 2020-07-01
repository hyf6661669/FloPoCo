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
            string inputSignal = join("inputX", i);
            string weightSignal = join("W", i);
            addInput(inputSignal, inputMsb - inputLsb + 1); // unsigned (always positive)
            addInput(weightSignal, inputMsb - inputLsb + 2); // signed + sign bit
            vhdl << tab << declare(join("signW", i)) << " <= " << weightSignal << "(" << inputPrec << ");" << endl;
            vhdl << tab << declareFixPoint(join("valueW", i), true, inputMsb, inputLsb) << " <= signed (";
            vhdl << weightSignal << "(" << inputPrec - 1 << " downto 0));" << endl;
            vhdl << tab << declareFixPoint(join("X", i), true, inputMsb, inputLsb) << " <= signed(";
            vhdl << join("inputX", i) << ");" << endl;
        }

        addOutput("A", inputPrec);

        // Mult (add with LNS representation)
        for (int i = 0; i < prevLayerWidth; i++) {
            vhdl << tab << declare(join("LogProduct", i), -inputLsb + inputMsb + 2) << " <= ";
            vhdl << join("X", i) << " + " << join("valueW", i) << ";" << endl;
        }

        // 2^x table
        cout << "inputLSB is: " << inputLsb << " and expLsb is: " << expLsb << endl;

        cout << "EXP VECTOR IS: " << endl;
        vector<mpz_class> exp_vector;
        for(int i=0; i<(1<<-inputLsb + 1); i++) {
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
            exp_vector.push_back(exp_x);
            cout << "; X = " << x << " exp = " << exp_x;
        }
        cout << endl;

        for (int i = 0; i < prevLayerWidth; i++) {
            string expInput = join("expInput", i);
            string expOutput = join("expOutput", i);
            declare(expInput, -inputLsb + 1);
            declare(expOutput, -expLsb + 1);
            vhdl << tab << expInput << " <= " << join("signW", i);
            vhdl << " & " << join("LogProduct", i) << "(" << -inputLsb + 1 << "downto 0);" << endl;

            Table::newUniqueInstance(this, expInput, expOutput, exp_vector,
                                        join("ExpTable", i), -inputLsb + 2,
                                        -expLsb + 1);
            // vhdl << tab << join("expX", i) << " <= signed(" << join("output", i) << ");" << endl;
            // vhdl << instance(expTable, join("expTable", i));
        }

        // Accumulator
        BitHeap* bh = new BitHeap(this, sumMsb, sumLsb);

        for (int i = 0; i < prevLayerWidth; i++)
            bh->addSignal(join("expOutput", i), 0);
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
