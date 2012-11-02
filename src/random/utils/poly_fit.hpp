#ifndef flopoco_random_poly_fit_hpp
#define flopoco_random_poly_fit_hpp

#include <Eigen/Dense>

namespace flopoco{
namespace random{
	
/*! Note that this is not supposed to replace sollya etc., it's just to get something
	that is somewhat quicker (and dirtier), particularly when deciding whether to
	split segments in adaptive quantisation. Later on the coefficients can be used
	as the seed to sollya's minimax, or thrown into supnorm.
*/

template<class T>
std::pair<std::vector<T>,T> LeastSquaresFit(
	const std::vector<std::pair<T,T> > &points,
	const std::vector<int> &monomials,
	const std::vector<T> &weights
){
	using namespace Eigen;
	
	typedef Matrix<double, Dynamic, Dynamic> Matrix;
	typedef Matrix<double, Dynamic, 1> Vector;
	
	Matrix A(points.size(), monomials.size());
	Vector B(points.size());
	
	for(unsigned i=0;i<points.size();i++){
		const T &w=weights[i];
		for(unsigned j=0;j<monomials.size();i++){
			system(i,j) = w*pow(points[i].first, monomials[j]);
		}
		B(i)=w*points[i].second;
	}
	
	Vector coeffs=A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);
	
	T error=(A*res-b).array().abs.max();
	
	std::vector<T> res(monomials.size());
	for(unsigned i=0;i<monomials.size();i++){
		res[i]=coeffs(i);
	}
	return std::make_pair(res,error);
}

template<class T,class F>
std::pair<std::vector<T>,T> LeastSquaresFit(
	const F &f, T a, T B, const std::vector<T> &nodes,
	const std::vector<int> &monomials
){
	T scale=b-a;
	std::vector<std::pair<T,T> > points(nodes.size());
	std::vector<T> weights(nodes.size());
	for(unsigned i=0;i<nodes.size();i++){
		points[i].second=a+scale*nodes[i];
		points[i].first=f(points[i].first);
		weights[i]=1.0;
	}
	return LeastSquaresFit(points, monomials, weights);
}

}; // random
}; // flopoco

#endif
