#ifndef flopoco_random_approx_2d_hpp
#define flopoco_random_approx_2d_hpp

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_multifit.h"
#include <vector>
#include <cmath>
#include <set>
#include <boost/function.hpp>
#include <boost/smart_ptr.hpp>
#include <iostream>

namespace flopoco
{
namespace random
{
namespace approx_2d
{

	typedef std::pair<int,int> monomial_t;	// Represents (x^first * y^second)
	typedef std::vector<monomial_t> monomial_basis_t;
	typedef std::vector<monomial_basis_t> monomial_basis_set_t;
	
	typedef boost::function<double(double,double)> func_t;
		
	double EvalMonomial(const monomial_t &m, double x, double y)
	{ return pow(x,m.first) * pow(y,m.second); }
	
	/* Create a basis that consists of all combinations of x^0..x^degree and y^0..y^degree */
	monomial_basis_t MakeFullBasis(int degree)
	{
		monomial_basis_t res;
		for(int xi=0;xi<=degree;xi++){
			for(int yi=0;yi<=degree;yi++){
				res.push_back(monomial_t(xi,yi));
			}
		}
		return res;
	}
	
	monomial_basis_set_t MakeOmitNBasis(monomial_basis_t original, int keep)
	{
		if(keep==original.size())
			return monomial_basis_set_t(1, original);
		
		std::set<monomial_basis_t> acc;
		for(int i=0;i<original.size();i++){
			monomial_basis_t tmp(original);
			tmp.erase(tmp.begin()+i);
			
			monomial_basis_set_t got=MakeOmitNBasis(tmp, keep);
			acc.insert(got.begin(), got.end());
		}
		
		return monomial_basis_set_t(acc.begin(), acc.end());
	}
	
	void PrintPolynomial(std::ostream &dst, const monomial_basis_t &basis, const std::vector<double> &coeffs)
	{
		assert(basis.size()==coeffs.size());
		for(int i=0;i<basis.size();i++){
			if(i!=0)
				dst<<" + ";
			dst<<coeffs[i];
			if(basis[i].first || basis[i].second)
				dst<<"*";
			if(basis[i].first!=0){
				if(basis[i].first==1){
					dst<<"x";
				}else{
					dst<<"x^"<<basis[i].first;
				}
				if(basis[i].second!=0){
					dst<<"*";
				}
			}
			if(basis[i].second==1){
				dst<<"y";
			}else if(basis[i].second>1){
				dst<<"y^"<<basis[i].second;
			}
		}
	}
	
	std::vector<double> MakeLinearSpacing(int npoints, double xmin, double xmax)
	{
		std::vector<double> res(npoints);
		for(int i=0;i<npoints;i++){
			res[i]=xmin + ((xmax-xmin) * double(i)) / (npoints-1);
		}
		return res;
	}
	
	struct segment_t
	{
		std::vector<monomial_t> monomials;
		std::vector<double> coeffs;
		double xmin, xmax, ymin, ymax;
		double worstAbsError;
		double meanSquaredError;
		func_t f;
	};
		
	boost::shared_ptr<gsl_matrix> BuildMeshX(
		const std::vector<double> &xpoints,
		const std::vector<double> &ypoints,
		const std::vector<monomial_t> &monomials
	){
		int n=xpoints.size()*ypoints.size();
		boost::shared_ptr<gsl_matrix> res(gsl_matrix_alloc(n, monomials.size()), gsl_matrix_free);
		
		for(int xx=0;xx<xpoints.size();xx++){
			for(int yy=0;yy<ypoints.size();yy++){
				for(int mm=0;mm<monomials.size();mm++){
					gsl_matrix_set(res.get(), xx*ypoints.size()+yy, mm, EvalMonomial(monomials[mm], xpoints[xx], ypoints[yy]));
				}
			}
		}
		
		return res;
	}
	
	boost::shared_ptr<gsl_vector> BuildMeshY(
		func_t f,
		const std::vector<double> &xpoints,
		const std::vector<double> &ypoints
	){
		int n=xpoints.size()*ypoints.size();
		boost::shared_ptr<gsl_vector> res(gsl_vector_alloc(n), gsl_vector_free);
		
		for(int xx=0;xx<xpoints.size();xx++){
			for(int yy=0;yy<ypoints.size();yy++){
				gsl_vector_set(res.get(), xx*ypoints.size()+yy, f(xpoints[xx], ypoints[yy]));
			}
		}
		
		return res;
	}
		
	segment_t SolveMesh(
		func_t f,
		const std::vector<double> &xpoints,
		const std::vector<double> &ypoints,
		const std::vector<monomial_t> &monomials
	){
		boost::shared_ptr<gsl_matrix> X(BuildMeshX(xpoints, ypoints, monomials));
		boost::shared_ptr<gsl_vector> y(BuildMeshY(f, xpoints,ypoints));
		
		int n=xpoints.size()*ypoints.size();
		boost::shared_ptr<gsl_multifit_linear_workspace> workspace(gsl_multifit_linear_alloc(n, monomials.size()), gsl_multifit_linear_free);
		
		boost::shared_ptr<gsl_matrix> cov(gsl_matrix_alloc(monomials.size(), monomials.size()), gsl_matrix_free);
		double chi2;
		
		boost::shared_ptr<gsl_vector> coeffs(gsl_vector_alloc(monomials.size()), gsl_vector_free);
		if(gsl_multifit_linear (X.get(), y.get(), coeffs.get(), cov.get(), &chi2, workspace.get()))
			throw std::string("Error code from gsl_multifit_linear.");
		
		boost::shared_ptr<gsl_vector> residuals(gsl_vector_alloc(n), gsl_vector_free);
		if(gsl_multifit_linear_residuals (X.get(),y.get(), coeffs.get(), residuals.get()))
			throw std::string("Error code from gsl_multifit_linear_residuals.");
		
		segment_t res;
		res.monomials=monomials;
		res.coeffs.resize(monomials.size());
		for(int i=0;i<monomials.size();i++){
			res.coeffs[i]=gsl_vector_get(coeffs.get(), i);
		}
		res.xmin=xpoints.front();
		res.xmax=xpoints.back();
		res.ymin=ypoints.front();
		res.ymax=ypoints.back();
		res.worstAbsError=0;
		res.meanSquaredError=0;
		for(int i=0;i<n;i++){
			double err=gsl_vector_get(residuals.get(), i);
			res.worstAbsError=std::max(res.worstAbsError, std::abs(err));
			res.meanSquaredError += err*err;
		}
		res.meanSquaredError=res.meanSquaredError/n;
		
		return res;
	}
	
	segment_t SolveMesh(
		func_t f,
		const std::vector<double> &xpoints,
		const std::vector<double> &ypoints,
		const monomial_basis_set_t &basisSet,
		double acceptableError=0.0
	){
		segment_t best;
		best.worstAbsError=DBL_MAX;
		
		for(int i=0;i<basisSet.size();i++){
			segment_t curr=SolveMesh(f, xpoints, ypoints, basisSet[i]);
			if(curr.worstAbsError<best.worstAbsError)
				best=curr;
			if(curr.worstAbsError<acceptableError)
				return curr;
		}
		
		return best;
	}
	
	
	struct quad_node_t
	{
		double xmin, xmax, ymin, ymax;
		
		boost::shared_ptr<segment_t> leaf;	// If leaf is non-null, then no more splitting happens
		
		double splitX, splitY;
		boost::shared_ptr<quad_node_t> q00, q01, q10, q11;
		
		int totalLeaves;
	};
	
	boost::shared_ptr<quad_node_t> BuildQuadTree(
		func_t f, const monomial_basis_set_t &basisSet, int pointsPerAxis, double maxAbsError,
		double xmin, double xmax, double ymin, double ymax
	){
		boost::shared_ptr<quad_node_t> res(new quad_node_t());
		res->xmin=xmin;
		res->xmax=xmax;
		res->ymin=ymin;
		res->ymax=ymax;
		
		std::vector<double> xpoints=MakeLinearSpacing(pointsPerAxis, xmin, xmax);
		std::vector<double> ypoints=MakeLinearSpacing(pointsPerAxis, ymin, ymax);
		
		segment_t sol=SolveMesh(f, xpoints, ypoints, basisSet, maxAbsError);
		
		if(sol.worstAbsError<maxAbsError){
			//fprintf(stderr, "[%lf,%lf] - [%lf,%lf], area=%lg\n", xmin, ymin, xmax, ymax, (ymax-ymin)*(xmax-xmin));
			res->leaf.reset(new segment_t(sol));
			res->totalLeaves=1;
		}else{
			res->splitX=(xmin+xmax)/2;
			res->splitY=(ymin+ymax)/2;
			
			res->q00=BuildQuadTree(f, basisSet, pointsPerAxis, maxAbsError, xmin, res->splitX, ymin, res->splitY);
			res->q01=BuildQuadTree(f, basisSet, pointsPerAxis, maxAbsError, xmin, res->splitX, res->splitY, ymax);
			res->q10=BuildQuadTree(f, basisSet, pointsPerAxis, maxAbsError, res->splitX, xmax, ymin, res->splitY);
			res->q11=BuildQuadTree(f, basisSet, pointsPerAxis, maxAbsError, res->splitX, xmax, res->splitY, ymax);
			res->totalLeaves=res->q00->totalLeaves+res->q01->totalLeaves+res->q10->totalLeaves+res->q11->totalLeaves;
		}
		return res;
	}
		

}; // approx_2d
}; // random
}; // flopoco

#endif
