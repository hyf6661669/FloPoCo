#include "decompose_period.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
	try{
		srand48(0);
		
		std::vector<std::pair<int,int> > generators;
		for(int i=0;i<16;i++){
			generators.push_back(std::make_pair(i+16,32+lrand48()%32));
		}
		
		for(int target_r=16;target_r<256;target_r++){
			std::vector<std::pair<int,int> > sol=flopoco::random::decompose_period::Execute(target_r, generators);
			std::cout<<target_r<< " : [";
			for(int i=0;i<sol.size();i++){
				if(i!=0)
					std::cout<<",";
				std::cout<<"("<<sol[i].first<<","<<sol[i].second<<")";
			}
			std::cout<<"], ";
			std::pair<mpz_class,mpz_class> periods=flopoco::random::decompose_period::period(sol);
			std::cout<<"log2(max)="<<mpz_sizeinbase(periods.first.get_mpz_t(),2)<<", log2(got)="<<mpz_sizeinbase(periods.second.get_mpz_t(),2)<<"\n";
		}
	}catch(std::string &msg){
		std::cerr<<"Exception : "<<msg<<"\n";
		return 1;
	}
	
	return 0;
}
