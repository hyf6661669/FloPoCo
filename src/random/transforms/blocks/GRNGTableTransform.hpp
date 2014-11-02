#ifndef random_grng_table_transform_hpp
#define random_grng_table_transform_hpp

namespace flopoco{
namespace random{

	void GRNGTableTransform_registerFactory();
	
	TableTransform *MakeGRNGTable(Target *target, int k, int wF, mpfr::mpreal stddev, std::string correction, std::string quantisation);
	
}; // random
}; // flopoco
#endif
