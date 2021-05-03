//
// Created by Philipp Keller on 18.03.21.
//
#include "IntelRCCM.hpp"
#include "IntelLcellTestOperator.hpp"
#include "Intel_LCELL.hpp"

#define S0(x) "S0(" + to_string(x) + ")"
#define S1(x) "S1(" + to_string(x) + ")"
#define A1(x) "A1(" + to_string(x) + ")"
#define B1(x) "B1(" + to_string(x) + ")"
#define E(x) "E(" + to_string(x) + ")"
#define A2(x) "A2(" + to_string(x) + ")"
#define Y(x) "Y(" + to_string(x) + ")"
#define CARRY(x) "carry(" + to_string(x) + ")"
#define A1_BIT "0101"
#define NOT_A1 "1010"
#define A2_BIT "0011"
#define NOT_A2 "1100"
#define B1_BIT "1100"
#define NOT_B1 "0011"
#define ZERO "0000"
#define ONE "1111"

using namespace std;

namespace flopoco
{
    IntelRCCM::IntelRCCM(Operator *parentOp, Target *target, int wIn, string type) : Operator(parentOp, target)
    {
        this->type = type;
        this->wIn = wIn;
        srcFileName="IntelRCCM";

        // definition of the name of the operator
        ostringstream name;
        name << "IntelRCCM_" << wIn << "_" << type;
        setName(name.str()); // See also setNameWithFrequencyAndUID()

        declare("carry", wIn + 1);

        addInput("S0", 1);
        addInput("S1", 1);
        addInput("A1", wIn);
        addInput("B1", wIn);
        addInput("E", wIn);
        addInput("A2", wIn);
        addOutput("Y", wIn + 1);

        string lut_mask_parts_a[4];
        string lut_mask_parts_b[4];
        string init_carry_parts[4];

        for(int i = 0; i < type.length(); i++){
            switch(type[i]){
                case '0':
                    lut_mask_parts_a[i] = A1_BIT;
                    lut_mask_parts_b[i] = B1_BIT;
                    init_carry_parts[i] = ZERO;
                    break;
                case '1':
                    lut_mask_parts_a[i] = A1_BIT;
                    lut_mask_parts_b[i] = NOT_B1;
                    init_carry_parts[i] = ONE;
                    break;
                case '2':
                    lut_mask_parts_a[i] = A1_BIT;
                    lut_mask_parts_b[i] = ZERO;
                    init_carry_parts[i] = ZERO;
                    break;
                case '3':
                    lut_mask_parts_a[i] = NOT_A1;
                    lut_mask_parts_b[i] = B1_BIT;
                    init_carry_parts[i] = ONE;
                    break;
                case '4':
                    lut_mask_parts_a[i] = NOT_A1;
                    lut_mask_parts_b[i] = ZERO;
                    init_carry_parts[i] = ONE;
                    break;
                case '5':
                    lut_mask_parts_a[i] = A2_BIT;
                    lut_mask_parts_b[i] = B1_BIT;
                    init_carry_parts[i] = ZERO;
                    break;
                case '6':
                    lut_mask_parts_a[i] = A2_BIT;
                    lut_mask_parts_b[i] = NOT_B1;
                    init_carry_parts[i] = ONE;
                    break;
                case '7':
                    lut_mask_parts_a[i] = A2_BIT;
                    lut_mask_parts_b[i] = ZERO;
                    init_carry_parts[i] = ZERO;
                    break;
                case '8':
                    lut_mask_parts_a[i] = NOT_A2;
                    lut_mask_parts_b[i] = B1_BIT;
                    init_carry_parts[i] = ONE;
                    break;
                case '9':
                    lut_mask_parts_a[i] = NOT_A2;
                    lut_mask_parts_b[i] = ZERO;
                    init_carry_parts[i] = ONE;
                    break;
                case 'a':
                    lut_mask_parts_a[i] = ONE;
                    lut_mask_parts_b[i] = B1_BIT;
                    init_carry_parts[i] = ZERO;
                    break;
                case 'b':
                    lut_mask_parts_a[i] = ONE;
                    lut_mask_parts_b[i] = NOT_B1;
                    init_carry_parts[i] = ONE;
                    break;
                case 'c':
                    lut_mask_parts_a[i] = ONE;
                    lut_mask_parts_b[i] = ZERO;
                    init_carry_parts[i] = ZERO;
                    break;
                default:
                    lut_mask_parts_a[i] = ZERO;
                    lut_mask_parts_b[i] = ZERO;
                    init_carry_parts[i] = ZERO;
            }
        }

        string lut_mask_a;
        string lut_mask_b;
        string init_carry;
        for(int i = 0; i < 4; i++){
            for(int j = 3; j >= 0; j--){
                lut_mask_a += lut_mask_parts_a[j][i];
                lut_mask_b += lut_mask_parts_b[j][i];
                init_carry += init_carry_parts[j][i];

            }
        }

        string init_carry_inv;
        for(char bit : init_carry){
            init_carry_inv += (bit == '0') ? '1' : '0';
        }

        string lut_mask = "\"0000000000000000" + lut_mask_a + "0000000000000000" + lut_mask_b + "\"";
        string lut_mask_init_carry = "\"0000000000000000" + init_carry_inv + "0000000000000000" + init_carry + "\"";

        addLCell(0, target, lut_mask_init_carry, S0(0), S1(0), "'0'", "'0'", "'0'", "'0'", "'0'", "open", CARRY(0));

        for(int i = 0; i < wIn; i++){
            addLCell(i+1, target, lut_mask, S0(0), S1(0), A1(i), B1(i),
                     E(i), A2(i), CARRY(i), Y(i), CARRY(i+1));
        }
        vhdl << tab << Y(wIn) << " <= " << CARRY(wIn) << ";" << endl;



    }

    void IntelRCCM::addLCell(int id, Target *target, string lut_mask, string dataa, string datab,
                             string datac, string datad, string datae, string dataf,
                             string cin, string sumout, string cout){

        auto *lcell = new Intel_LCELL(this, target, lut_mask);
        if(dataa == "'0'")
            inPortMapCst("dataa", dataa);
        else
            inPortMap("dataa", dataa);
        if(datab == "'0'")
            inPortMapCst("datab", datab);
        else
            inPortMap("datab", datab);
        if(datac == "'0'")
            inPortMapCst("datac", datac);
        else
            inPortMap("datac", datac);
        if(datad == "'0'")
            inPortMapCst("datad", datad);
        else
            inPortMap("datad", datad);
        if(datae == "'0'")
            inPortMapCst("datae", datae);
        else
            inPortMap("datae", datae);
        if(dataf == "'0'")
            inPortMapCst("dataf", dataf);
        else
            inPortMap("dataf", dataf);
        inPortMapCst("datag", "'0'");
        if(cin == "'0'")
            inPortMapCst("cin", cin);
        else
            inPortMap("cin", cin);
        inPortMapCst("sharein","'0'");
        outPortMap("sumout", sumout);
        outPortMap("cout", cout);
        outPortMap("combout", "open");
        outPortMap("shareout", "open");

        vhdl << lcell->primitiveInstance("lcell" + to_string(id)) << endl;
    }


    OperatorPtr IntelRCCM::parseArguments(OperatorPtr parentOp, Target *target, vector<string> &args)
    {
        int wIn;
        string type;

        UserInterface::parseInt(args, "wIn", &wIn);
        UserInterface::parseString(args, "type", &type);

        return new IntelRCCM(parentOp,target,wIn,type);
    }

    void IntelRCCM::registerFactory(){
        UserInterface::add("IntelRCCM", // name
                           "Test Operator for Intel ALM Primitives",
                           "Primitives", // categories
                           "",
                           "wIn(int): Size of input in bits; \
                           type(string): 4-digit operation type",
                           "",
                           IntelRCCM::parseArguments,
                           IntelRCCM::unitTest
        ) ;
    }

    void IntelRCCM::emulate (TestCase* tc)
    {
        mpz_class s0 = tc->getInputValue("S0");
        mpz_class s1 = tc->getInputValue("S1");
        mpz_class a1 = tc->getInputValue("A1");
        mpz_class a2 = tc->getInputValue("A2");
        mpz_class b1 = tc->getInputValue("B1");
        mpz_class result;

        switch(type[(2 * s1.get_si() + s0.get_si())]){
            case '0':
                result = a1 + b1;
                break;
            case '1':
                result = a1 - b1;
                break;
            case '2':
                result = a1;
                break;
            case '3':
                result = -a1 + b1;
                break;
            case '4':
                result = -a1;
                break;
            case '5':
                result = a2 + b1;
                break;
            case '6':
                result = a2 - b1;
                break;
            case '7':
                result = a2;
                break;
            case '8':
                result = -a2 + b1;
                break;
            case '9':
                result = -a2;
                break;
            case 'a':
                result = b1;
                break;
            case 'b':
                result = -b1;
                break;
            case 'c':
                result = 0;
                break;
            default:
                result = -1;
        }

        tc->addExpectedOutput("Y", result);

    }

    TestList IntelRCCM::unitTest(int)
    {
        TestList testStateList;
        vector<pair<string,string>> paramList;

        for(int a = 0; a <= 12; a++){
            for(int b = 0; b <= 12; b++){
                if(b == a) continue;
                for(int c = 0; c <= 12; c++){
                    if(c == a || c == b) continue;
                    for(int d = 0; d <= 12; d++){
                        if(d == a || d == b || d == c) continue;

                        paramList.push_back(make_pair("type", string(1, (char)(((a >= 10) ? 'W' : '0') + a)) + string(1, (char)(((b >= 10) ? 'W' : '0') + b)) + string(1, (char)(((c >= 10) ? 'W' : '0') + c)) + string(1, (char)(((d >= 10) ? 'W' : '0') + d))));
                        paramList.push_back(make_pair("wIn", to_string((rand() % 2) + 4)));
                        testStateList.push_back(paramList);
                        paramList.clear();
                    }
                }
            }
        }
        return testStateList;
    }

}