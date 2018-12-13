//
// Created by nfiege on 10/12/18.
//

#include "AddressCalculation.hpp"

#include "NeuralNetworks/Utility/BitheapWrapper.hpp"

#include <vector>

using namespace std;

namespace flopoco {
    AddressCalculation::AddressCalculation(Target *target, unsigned int startAddress, unsigned int innerCounterWidth,
                                           unsigned int outerCounterWidth, int bytesPerAccess, int bytesPerFeature)
		: Operator(target){
    	cout << "#q# 1" << endl;

        // definition of the source file name, used for info and error reporting using REPORT
        srcFileName="WeightFetcher";

        // use numeric_std
        useNumericStd();

        // author
        setCopyrightString("Nicolai Fiege, 2018");
		cout << "#q# 2" << endl;

        // definition of the name of the operator
        ostringstream name;
        name << "AddressCalculation_" << startAddress << "_" << innerCounterWidth << "_" << outerCounterWidth << "_"
             << (bytesPerAccess>0?to_string(bytesPerAccess):"m"+to_string(-bytesPerAccess)) << "_" << bytesPerFeature;
        setName(name.str());
		cout << "#q# 3" << endl;

        if(bytesPerAccess==0) THROWERROR("bytes per memory access == 0 => that makes no sense");

        ///////////
        // ports //
        ///////////
        if(innerCounterWidth > 0)
        	addInput("InnerCounter",innerCounterWidth);
		if(outerCounterWidth > 0)
        	addInput("OuterCounter",outerCounterWidth);
		cout << "#q# 4" << endl;
        int bytesPerAccessWidth;
        if(bytesPerAccess<0)
        {
			bytesPerAccessWidth = -bytesPerAccess;
			addInput("BytesPerAccess",bytesPerAccessWidth);
        }
        else
		{
			bytesPerAccessWidth = (int)ceil(log2(bytesPerAccess+1));
			this->vhdl << declare("BytesPerAccess",bytesPerAccessWidth) << " <= std_logic_vector(to_unsigned("
								  << bytesPerAccess << "," << bytesPerAccessWidth << "));" << endl;
		}
		cout << "#q# 5" << endl;
		auto bytesPerFeatureWidth = (int)ceil(log2(bytesPerFeature+1));
		int p1Width = bytesPerAccessWidth + innerCounterWidth;
		int p2Width = bytesPerFeatureWidth + outerCounterWidth;
		cout << "#q# 6" << endl;

		this->vhdl << declare("BytesPerFeature",bytesPerFeatureWidth) << " <= std_logic_vector(to_unsigned("
				   << bytesPerFeature << "," << bytesPerFeatureWidth << "));" << endl;

		if(innerCounterWidth > 0)
		{
			this->vhdl << declare("P1",p1Width) << " <= std_logic_vector(unsigned(InnerCounter) * unsigned(BytesPerAccess));" << endl;
		} else {
			this->vhdl << declare("P1",1) << " <= \"0\";" << endl;
			p1Width = 1;
		}

		if(outerCounterWidth > 0)
		{
			this->vhdl << declare("P2",p2Width) << " <= std_logic_vector(unsigned(OuterCounter) * unsigned(BytesPerFeature));" << endl;
		}
		else {
			this->vhdl << declare("P2",1) << " <= \"0\";" << endl;
			p2Width = 1;
		}

		cout << "#q# 7" << endl;

		nextCycle();

		this->vhdl << declare("StartAddress",32) << " <= \"" << AddressCalculation::uintTo32BitString(startAddress)
				   << "\";" << endl;
		cout << "#q# 8" << endl;

		vector<int> wIns;
		wIns.emplace_back(32);
		wIns.emplace_back(p1Width);
		wIns.emplace_back(p2Width);
		cout << "#q# 9" << endl;
		auto bh = new BitheapWrapper(target,wIns,3,vector<bool>(3,false),vector<int>(3,0));
		cout << "#q# 10" << endl;
		this->addSubComponent(bh);
		this->inPortMap(bh,"X0","StartAddress");
		this->inPortMap(bh,"X1","P1");
		this->inPortMap(bh,"X2","P2");
		this->outPortMap(bh,"R","NextAddressTemp",true);
		this->vhdl << instance(bh,"AddressAdder");
		this->syncCycleFromSignal("NextAddressTemp");
		cout << "#q# 11" << endl;

		addOutput("NextAddress",32);
		this->vhdl << "NextAddress <= NextAddressTemp(31 downto 0);" << endl;
		cout << "#q# 12" << endl;
    }

	std::string AddressCalculation::uintTo32BitString(unsigned int i) {
		cout << "#q# BOOP 1" << endl;
		string s = "";
		for(int j=31; j>=0; j--)
		{
			s += to_string((i >> j) & 1);
		}
		cout << "#q# BOOP 2" << endl;
		return s;
	}

}