#ifndef flopoco_random_approx_2d_hpp
#define flopoco_random_approx_2d_hpp

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_multifit.h"
#include "gsl/gsl_blas.h"
#include <vector>
#include <cmath>
#include <set>
#include <functional>
#include <iostream>
#include <cassert>
#include <memory>

namespace flopoco
{
namespace random
{
namespace approx_2d
{

	typedef std::pair<int,int> monomial_t;	// Represents (x^first * y^second)
	typedef std::vector<monomial_t> monomial_basis_t;
	typedef std::vector<monomial_basis_t> monomial_basis_set_t;
	
	typedef std::function<double(double,double)> func_t;
		
	double EvalMonomial(const monomial_t &m, double x, double y, double rotation)
	{
		if(rotation){
			double tx=cos(rotation)*x-sin(rotation)*y;
			double ty=sin(rotation)*x+cos(rotation)*y;
			x=tx;
			y=ty;
		}
		return pow(x,m.first) * pow(y,m.second);
	}
	
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
	
	std::vector<double> MakeLinearSpacingExcl(int npoints, double xmin, double xmax)
	{
		std::vector<double> res(npoints);
		for(int i=0;i<npoints;i++){
			res[i]=xmin + ((xmax-xmin) * double(i+1)) / (npoints+1);
		}
		return res;
	}
	
	enum error_metric_t
	{
		error_metric_worst_absolute=0,
		error_metric_worst_relative=1
	};
	
	struct error_target_t
	{
		error_metric_t metric;
		double target;
	};
	
	struct error_t
	{
		error_t(double weight)
			: worstRelError(0)
			, worstAbsError(0)
			, rawSumSquaredAbsError(0)
			, rawSumSquaredRelError(0)
			, weightedSumSquaredAbsError(0)
			, weightedSumSquaredRelError(0)
			, totalPoints(0)
			, totalWeight(weight)
		{}
			
		void add(double ref, double got)
		{
			double absErr=std::abs(ref-got);
			double relErr=std::abs(absErr/ref);
			
			worstAbsError=std::max(absErr, worstAbsError);
			worstRelError=std::max(relErr, worstRelError);
			rawSumSquaredAbsError+=absErr*absErr;
			rawSumSquaredRelError+=relErr*relErr;
			weightedSumSquaredAbsError+=absErr*absErr;
			weightedSumSquaredRelError+=relErr*relErr;
			totalPoints=totalPoints+1;
		}
		
		double worstRelError;
		double worstAbsError;
		double rawSumSquaredAbsError;
		double rawSumSquaredRelError;
		double weightedSumSquaredAbsError;
		double weightedSumSquaredRelError;
		double totalPoints;		// Number of points considered
		double totalWeight;		// weight (area) of the input area
	};
	
	error_t combine(const error_t &a, const error_t &b)
	{
		error_t res(a.totalWeight+b.totalWeight);
		res.worstRelError=std::max(a.worstRelError, b.worstRelError);
		res.worstAbsError=std::max(a.worstAbsError, b.worstAbsError);
		res.rawSumSquaredAbsError=a.rawSumSquaredAbsError+b.rawSumSquaredAbsError;
		res.rawSumSquaredRelError=a.rawSumSquaredRelError+b.rawSumSquaredRelError;
		res.weightedSumSquaredAbsError=(a.totalWeight*a.weightedSumSquaredAbsError+b.totalWeight*b.weightedSumSquaredAbsError)/(a.totalWeight+b.totalWeight);
		res.weightedSumSquaredRelError=(a.totalWeight*a.weightedSumSquaredRelError+b.totalWeight*b.weightedSumSquaredRelError)/(a.totalWeight+b.totalWeight);
		res.totalPoints=a.totalPoints+b.totalPoints;
		res.totalWeight=a.totalWeight+b.totalWeight;
		return res;
	}
	
	bool better_than(error_target_t target, const error_t &a, const error_t &b)
	{
		if(target.metric==error_metric_worst_absolute){
			return a.worstAbsError < b.worstAbsError;
		}else if(target.metric==error_metric_worst_relative){
			return a.worstRelError < b.worstRelError;
		}else{
			throw std::invalid_argument("Unknown error target.");
		}
	}
	
	bool good_enough(error_target_t target, const error_t &e)
	{
		if(target.metric==error_metric_worst_absolute){
			return e.worstAbsError < target.target;
		}else if(target.metric==error_metric_worst_relative){
			return e.worstRelError < target.target;
		}else{
			throw std::invalid_argument("Unknown error target.");
		}
	}
	
	struct segment_t
	{
		std::vector<monomial_t> monomials;
		std::vector<double> coeffs;
		double xmin, xmax, ymin, ymax;
		error_t error;
		func_t f;
	};
		
	std::shared_ptr<gsl_matrix> BuildMeshX(
		const std::vector<double> &xpoints,
		const std::vector<double> &ypoints,
		const std::vector<monomial_t> &monomials,
		double rotation
	){
		int n=xpoints.size()*ypoints.size();
		std::shared_ptr<gsl_matrix> res(gsl_matrix_alloc(n, monomials.size()), gsl_matrix_free);
		
		for(int xx=0;xx<xpoints.size();xx++){
			for(int yy=0;yy<ypoints.size();yy++){
				for(int mm=0;mm<monomials.size();mm++){
					gsl_matrix_set(res.get(), xx*ypoints.size()+yy, mm, EvalMonomial(monomials[mm], xpoints[xx], ypoints[yy], rotation));
				}
			}
		}
		
		return res;
	}
	
	std::shared_ptr<gsl_vector> BuildMeshY(
		func_t f,
		const std::vector<double> &xpoints,
		const std::vector<double> &ypoints
	){
		int n=xpoints.size()*ypoints.size();
		std::shared_ptr<gsl_vector> res(gsl_vector_alloc(n), gsl_vector_free);
		
		for(int xx=0;xx<xpoints.size();xx++){
			for(int yy=0;yy<ypoints.size();yy++){
				gsl_vector_set(res.get(), xx*ypoints.size()+yy, f(xpoints[xx], ypoints[yy]));
			}
		}
		
		return res;
	}
		
	std::shared_ptr<segment_t> SolveMesh(
		func_t f,
		const std::vector<double> &xpoints,
		const std::vector<double> &ypoints,
		double rotation,
		const std::vector<monomial_t> &monomials,
		error_metric_t metric
	){
		//fprintf(stderr, "Solve\n");
		
		std::shared_ptr<gsl_matrix> X(BuildMeshX(xpoints, ypoints, monomials, rotation));
		std::shared_ptr<gsl_vector> y(BuildMeshY(f, xpoints,ypoints));
		
		int n=xpoints.size()*ypoints.size();
		std::shared_ptr<gsl_multifit_linear_workspace> workspace(gsl_multifit_linear_alloc(n, monomials.size()), gsl_multifit_linear_free);
		
		std::shared_ptr<gsl_matrix> cov(gsl_matrix_alloc(monomials.size(), monomials.size()), gsl_matrix_free);
		double chi2;
		
		std::shared_ptr<gsl_vector> coeffs(gsl_vector_alloc(monomials.size()), gsl_vector_free);
		
		if(metric==error_metric_worst_relative){
			std::shared_ptr<gsl_vector> w(gsl_vector_alloc(n), gsl_vector_free);
			for(unsigned i=0;i<n;i++){
				double r=gsl_vector_get(y.get(),i);
				if(r==0)
					throw std::string("Cannot fit relative error to zero.");
				gsl_vector_set(w.get(), i, 1.0/gsl_vector_get(y.get(), i));
			}
			
			if(gsl_multifit_wlinear (X.get(), w.get(), y.get(), coeffs.get(), cov.get(), &chi2, workspace.get()))
				throw std::string("Error code from gsl_multifit_wlinear.");
		}else if(metric==error_metric_worst_absolute){
			if(gsl_multifit_linear (X.get(), y.get(), coeffs.get(), cov.get(), &chi2, workspace.get()))
				throw std::string("Error code from gsl_multifit_linear.");
		}else{
			throw std::string("Unknown error metric.");
		}
		
		for(unsigned i=0;i<monomials.size();i++){
			for(unsigned j=0;j<monomials.size();j++){
				//fprintf(stderr, " % 8lg", gsl_matrix_get(cov.get(), i, j));
			}
			//fprintf(stderr, "\n");
		}
		
		std::shared_ptr<gsl_vector> est(gsl_vector_alloc(n), gsl_vector_free);
		gsl_blas_dgemv(CblasNoTrans, 1.0, X.get(), coeffs.get(), 0.0, est.get());
		
		
		std::vector<double> coeffsD(monomials.size());
		for(int i=0;i<monomials.size();i++){
			coeffsD[i]=gsl_vector_get(coeffs.get(), i);
			//fprintf(stderr, "c%d = %lg\n", i, coeffsD[i]);
		}
		
		error_t err((xpoints.back()-xpoints.front())*(ypoints.back()-ypoints.front()));
		for(int i=0;i<n;i++){
			double ref=gsl_vector_get(y.get(), i);
			double got=gsl_vector_get(est.get(), i);
			err.add(ref,got);
		}
		//fprintf(stderr, "  points=%d, n=%d\n", err.totalPoints, n);
		
		segment_t res={
			monomials,
			coeffsD,
			xpoints.front(),
			xpoints.back(),
			ypoints.front(),
			ypoints.back(),
			err,
			f
		};
		
		return std::make_shared<segment_t>(res);
	}
	
	std::shared_ptr<segment_t> SolveMesh(
		func_t f,
		const std::vector<double> &xpoints,
		const std::vector<double> &ypoints,
		std::set<double> rotations,
		const monomial_basis_set_t &basisSet,
		error_target_t target
	){
		std::shared_ptr<segment_t> best;
		
		for(double rotation : rotations){
			for(int i=0;i<basisSet.size();i++){
				auto curr=SolveMesh(f, xpoints, ypoints, rotation, basisSet[i], target.metric);
				if(!best){
					best=curr;
				}else if(better_than(target, curr->error, best->error)){
					best=curr;
				}
				
				if(good_enough(target, best->error))
					return best;
			}
		}
		
		return best;
	}
	
	
	struct quad_node_t
	{
		quad_node_t(double _xmin, double _xmax, double _ymin, double _ymax)
			: xmin(_xmin), xmax(_xmax), ymin(_ymin), ymax(_ymax)
			, error((_xmax-_xmin)*(_ymax-_ymin))
		{}
		
		double xmin, xmax, ymin, ymax;
		
		std::shared_ptr<segment_t> leaf;	// If leaf is non-null, then no more splitting happens
		
		double splitX, splitY;
		std::shared_ptr<quad_node_t> q00, q01, q10, q11;
		
		error_t error;
		
		int totalLeaves;
		int maxDepth, minDepth;
	};
	
	std::shared_ptr<quad_node_t> BuildQuadTree(
		func_t f, const monomial_basis_set_t &basisSet, int pointsPerAxis,
		double xmin, double xmax, double ymin, double ymax, std::set<double> rotation,
		error_target_t target
	){
		std::shared_ptr<quad_node_t> res(new quad_node_t(xmin,xmax,ymin,ymax));
		
		std::vector<double> xpoints=MakeLinearSpacingExcl(pointsPerAxis, xmin, xmax);
		std::vector<double> ypoints=MakeLinearSpacingExcl(pointsPerAxis, ymin, ymax);
		
		auto sol=SolveMesh(f, xpoints, ypoints, rotation, basisSet, target);
			
		if(good_enough(target, sol->error)){
			//fprintf(stderr, "[%lf,%lf] - [%lf,%lf], area=%lg, wRel=%lg, wAbs=%lg (points=%lg)\n", xmin, ymin, xmax, ymax, (ymax-ymin)*(xmax-xmin), sol->error.worstRelError, sol->error.worstAbsError, sol->error.totalPoints);
			res->leaf=sol;
			res->totalLeaves=1;
			res->maxDepth=0;
			res->minDepth=0;
		}else{
			res->splitX=(xmin+xmax)/2;
			res->splitY=(ymin+ymax)/2;
			
			res->q00=BuildQuadTree(f, basisSet, pointsPerAxis, xmin, res->splitX, ymin, res->splitY, rotation, target);
			res->q01=BuildQuadTree(f, basisSet, pointsPerAxis, xmin, res->splitX, res->splitY, ymax, rotation, target);
			res->q10=BuildQuadTree(f, basisSet, pointsPerAxis, res->splitX, xmax, ymin, res->splitY, rotation, target);
			res->q11=BuildQuadTree(f, basisSet, pointsPerAxis, res->splitX, xmax, res->splitY, ymax, rotation, target);
			res->totalLeaves=res->q00->totalLeaves+res->q01->totalLeaves+res->q10->totalLeaves+res->q11->totalLeaves;
			res->maxDepth=1+std::max(std::max(res->q00->maxDepth,res->q01->maxDepth),std::max(res->q10->maxDepth,res->q11->maxDepth));
			res->minDepth=1+std::min(std::min(res->q00->maxDepth,res->q01->minDepth),std::max(res->q10->minDepth,res->q11->minDepth));
			res->error=combine(combine(res->q00->error,res->q01->error),combine(res->q10->error,res->q11->error));
		}
		return res;
	}
		

}; // approx_2d
}; // random
}; // flopoco

#endif
