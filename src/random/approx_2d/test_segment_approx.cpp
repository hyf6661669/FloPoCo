#include "approx_2d.hpp"

#include "gsl/gsl_sf_bessel.h"

using namespace flopoco::random::approx_2d;

double f_sin_mul_exp(double x, double y)
{
	return sin(x)*exp(y);
}

double f_sin_div_x(double x, double y)
{
	return sin(x*4)*(y+0.1);
}

double f_bessel(double x, double y)
{
	return gsl_sf_bessel_Inu(x,y);
}


/*
	x,y
	
	1, x, y, x^2
    1, x, y, y^2
	1, x, y, x*y
 
    1, x, y, x^2, x^3
	1, x, y, x^2, x^2*y
    1, x, y, y^2, y^2*x
    1, x, y, y^2, y^3
	1, x, y, m1=(x|y)*(x|y), m2=m1*(x|y)

	1, x, y, x^2, x^3
	1, x, y, x^2, x^2*y
    1, x, y, y^2, y^2*x
    1, x, y, y^2, y^3
	1, x, y, m1=x*x, m2=x*y, m3=y*y, 


	


*/

int main(int argc, char *argv[])
{
	std::vector<double> xpoints=MakeLinearSpacing(7, 0.0, 1.0);
	std::vector<double> ypoints=MakeLinearSpacing(7, 0.0, 1.0);
	
	//func_t f=f_bessel;
	
	func_t f=f_sin_mul_exp;
	
	//for(int degree=2;degree<=2;degree++){
	//	monomial_basis_t basis=MakeFullBasis(degree);
		
		/*
	for(int keep=basis.size();keep>=basis.size();keep--){
			monomial_basis_set_t ss=MakeOmitNBasis(basis, keep);
		
			std::shared_ptr<quad_node_t> quad=BuildQuadTree(
				f, ss, 10,
				-8, 8, -8, 8,
				error_target_t{error_metric_worst_absolute, ldexp(1.0,-12)}
			);
			std::cerr<<"degree="<<degree<<", keep="<<keep<<", totalLeaves="<<quad->totalLeaves<<", depth="<<quad->minDepth<<".."<<quad->maxDepth<<"\n";
		}		
	}*/
		
	monomial_basis_t basis1{
		std::make_pair(0,0),
		std::make_pair(0,1),std::make_pair(1,0),
		std::make_pair(2,0), std::make_pair(1,1),
		std::make_pair(3,0),std::make_pair(2,1)
	};
    monomial_basis_t basis2{
		std::make_pair(0,0),
		std::make_pair(0,1),std::make_pair(1,0),
		std::make_pair(2,0), std::make_pair(1,1),std::make_pair(0,2),
		std::make_pair(3,0)
	};
    monomial_basis_t basis3{
		std::make_pair(0,0),
		std::make_pair(0,1),std::make_pair(1,0),
		std::make_pair(2,0), std::make_pair(1,1),std::make_pair(0,2),
		std::make_pair(2,1)
	};
	monomial_basis_set_t bset{basis1,basis2,basis3};
	
	const double PI=3.1415926535897932384626433832795;
	
	for(int rots=1;rots<=64;rots*=2){
		std::set<double> rotations;
		for(int i=0;i<rots;i++){
			rotations.insert((PI*i)/rots);
		}
		
		std::shared_ptr<quad_node_t> quad=BuildQuadTree(
			f, bset, 10,
			-4, 4, -4, 4, rotations,
			error_target_t{error_metric_worst_absolute, ldexp(1.0,-12)}
		);
		std::cerr<<"rots="<<rots<<", totalLeaves="<<quad->totalLeaves<<", depth="<<quad->minDepth<<".."<<quad->maxDepth<<"\n";
	}
	
	
	return 0;
}
