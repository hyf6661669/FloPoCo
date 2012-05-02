#include "correct_distribution.hpp"

#include "random/distributions/gaussian_distribution.hpp"
#include "random/distributions/make_table_approximation.hpp"
#include "random/distributions/clt_distribution.hpp"


using namespace flopoco::random;

	

int main(int argc, char *argv[])
{
	try{
		ContinuousDistribution<double>::TypePtr target=boost::make_shared<GaussianDistribution<double> >();
		
		for(unsigned w=3;w<=16;w++){
			for(unsigned k=2;k<=8;k+=2){
				TableDistribution<double>::TypePtr table=MakeCLTDistribution<double>(w, k);
							
				for(unsigned d=1;d<=9;d+=2){
					std::vector<double> poly=FindPolynomialCorrection<double>(table, target, d);
					
					TableDistribution<double>::TypePtr corrected=table->ApplyPolynomial(poly);
					
					poly.resize(10, 0);	// extend to degree 9 with zeros
					fprintf(stdout, "clt, %d, %d, %d, %.12lg, %.12lg, %.12lg, %.12lg, %.12lg\n",
						w, k, d, poly[1], poly[3], poly[5], poly[7], poly[9]
					);
					fprintf(stderr, "  m2=%g, m4=%lg\n", corrected->StandardMoment(2), corrected->StandardMoment(4));
				}
			}
		}
		
		for(unsigned n=16;n<=65536;n*=2){
			TableDistribution<double>::TypePtr table=MakeTableApproximation<double>(target, n);
			
			for(unsigned d=1;d<=9;d+=2){
				std::vector<double> poly=FindPolynomialCorrection<double>(table, target, d);
				
				TableDistribution<double>::TypePtr corrected=table->ApplyPolynomial(poly);
				
				poly.resize(10, 0);	// extend to degree 9 with zeros
				fprintf(stdout, "table_real, %d, %d, %.12lg, %.12lg, %.12lg, %.12lg, %.12lg\n",
					n, d, poly[1], poly[3], poly[5], poly[7], poly[9]
				);
				
				fprintf(stderr, "  m2=%g, m4=%lg\n", table->StandardMoment(2), table->StandardMoment(4));
			}
		}
		
		return 0;
	}catch(std::exception &e){
		std::cerr<<"Caught Exception : "<<e.what()<<"\n";
		return 1;
	}catch(...){
		std::cerr<<"Caught unexpected exception.";
		return 1;
	}
}
