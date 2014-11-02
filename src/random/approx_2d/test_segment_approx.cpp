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

int main(int argc, char *argv[])
{
	std::vector<double> xpoints=MakeLinearSpacing(32, 0.1, 1.0);
	std::vector<double> ypoints=MakeLinearSpacing(32, 0.1, 1.0);
	
	func_t f=f_bessel;
	
	for(int degree=2;degree<=3;degree++){
		monomial_basis_t basis=MakeFullBasis(degree);
		
		for(int keep=basis.size();keep>2;keep--){
			monomial_basis_set_t ss=MakeOmitNBasis(basis, keep);
		
			segment_t seg=SolveMesh(f, xpoints, ypoints, basis);
			
			boost::shared_ptr<quad_node_t> quad=BuildQuadTree(
				f, ss, 10, pow(2.0,-10),
				0.1, 1.0, 0.1, 1.0
			);
			std::cerr<<"  keep="<<keep<<", totalLeaves="<<quad->totalLeaves<<"\n";
		}		
	}
	
	
	
	return 0;
}
