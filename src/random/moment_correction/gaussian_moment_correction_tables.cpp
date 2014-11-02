#include "correct_distribution.hpp"

#include "random/distributions/gaussian_distribution.hpp"
#include "random/distributions/make_table_approximation.hpp"
#include "random/distributions/clt_distribution.hpp"


using namespace flopoco::random;

	

int main(int argc, char *argv[])
{
	try{
		ContinuousDistribution<double>::TypePtr target=boost::make_shared<GaussianDistribution<double> >();
		
		/*
		for(unsigned w=3;w<=16;w++){
			for(unsigned k=2;k<=8;k+=2){
				TableDistribution<double>::TypePtr table=MakeCLTDistribution<double>(w, k);
							
				for(unsigned d=1;d<=7;d+=2){
					std::vector<double> poly=FindPolynomialCorrection<double>(table, target, d);
					
					TableDistribution<double>::TypePtr corrected=table->ApplyPolynomial(poly);
					
					poly.resize(10, 0);	// extend to degree 9 with zeros
					fprintf(stdout, "clt, %d, %d, %d, %.12lg, %.12lg, %.12lg, %.12lg, %.12lg\n",
						w, k, d, poly[1], poly[3], poly[5], poly[7], poly[9]
					);
					fprintf(stderr, "  m2=%lg -> %g, m4=%lg -> %lg\n", table->RawMoment(2), corrected->RawMoment(2), table->RawMoment(4), corrected->RawMoment(4));
				}
			}
		}*/
		
		for(unsigned n=16;n<=(1<<20);n*=2){
			TableDistribution<double>::TypePtr table=MakeTableApproximation<double>(target, n);
			
			for(unsigned d=1;d<=7;d+=2){
				std::vector<double> poly=FindPolynomialCorrection<double>(table, target, d);
				
				TableDistribution<double>::TypePtr corrected=table->ApplyPolynomial(poly);
				
				double error[9];
				double metric=0.0;
				for(int i=2;i<=8;i+=2){
					if(d+1>=i){
						metric += pow(corrected->StandardMoment(i) - target->StandardMoment(i),2);
					}
					error[i]=(corrected->StandardMoment(i) - target->StandardMoment(i)) / target->StandardMoment(i);
				}
				metric=sqrt(metric);
				
				
				poly.resize(10, 0);	// extend to degree 9 with zeros
				fprintf(stdout, "table_real, %d, %d, %.12lg, %.12lg, %.12lg, %.12lg, %.12lg,  %.12lg, %.12lg, %.12lg, %.12lg, %.12lg\n",
					n, d, poly[1], poly[3], poly[5], poly[7], poly[9],
					metric, error[2], error[4], error[6], error[8]
				);
				
				fprintf(stderr, "  m2=%g->%g, m4=%lg->%g, m6=%g->%g, m8=%lg->%g, \n",
					table->StandardMoment(2), corrected->StandardMoment(2),
					table->StandardMoment(4), corrected->StandardMoment(4),
					table->StandardMoment(6), corrected->StandardMoment(6),
					table->StandardMoment(8), corrected->StandardMoment(8)
				);
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
