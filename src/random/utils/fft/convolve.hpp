#ifndef flopoco_random_utils_convolve_hpp
#define flopoco_random_utils_convolve_hpp

namespace flopoco
{
namespace random
{
	template<class T>
	std::vector<T> convolve(const std::vector<T> &a, const std::vector<T> &b);
	
	template<class T>
	std::vector<T> self_convolve(const std::vector<T> &a, unsigned k);
};
};

#endif
