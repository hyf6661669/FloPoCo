#ifndef flopoco_random_utils_mpfr_vec_hpp
#define flopoco_random_utils_mpfr_vec_hpp

#include <stdio.h>
#include <mpfr.h>
#include <vector>
#include <limits.h>
#include <stdlib.h>
#include <mpfr.h>
#include <assert.h>

namespace flopoco
{
namespace random
{

class MPFRVec
{
private:
	int len, allocLen;
	mpfr_t *values;
public:
	MPFRVec()
		: len(0)
		, allocLen(0)
		, values(NULL)
	{}

	MPFRVec(int size, int prec)
		: len(size)
		, allocLen(size)
		, values((mpfr_t*)malloc(sizeof(mpfr_t)*len))
	{
		for(int i=0;i<size;i++){
			mpfr_init2(values[i], prec);
		}
	}
	
	MPFRVec(const MPFRVec &o)
		: len(o.len)
		, allocLen(o.len)
		, values((mpfr_t*)malloc(sizeof(mpfr_t)*len))
	{
		for(int i=0;i<len;i++){
			mpfr_init2(values[i], mpfr_get_prec(o.values[i]));
			mpfr_set(values[i], o.values[i], MPFR_RNDN);
		}
	}
	
	MPFRVec &operator=(const MPFRVec &o)
	{
		if(this!=&o){
			clear();
			len=o.len;
			allocLen=o.len;
			values=(mpfr_t*)realloc(values, sizeof(mpfr_t)*len);
			for(int i=0;i<len;i++){
				mpfr_init2(values[i], mpfr_get_prec(o.values[i]));
				mpfr_set(values[i], o.values[i], MPFR_RNDN);
			}
		}
		return *this;
	}
	
	void resize(int newLen, int prec)
	{
		if(newLen<len){
			while(newLen<len){
				mpfr_clear(values[len-1]);
				--len;
			}
		}
		if(newLen>len){
			if(newLen > allocLen){
				values=(mpfr_t*)realloc(values, sizeof(mpfr_t)*newLen);
				allocLen=newLen;
			}
			while(newLen>len){
				mpfr_init2(values[len], prec);
				++len;
			}
		}
	}
	
	void clear()
	{
		resize(0, 0);
	}
	
	void push_back(mpfr_t x)
	{
		if(allocLen <= len){
			allocLen=std::max(8, allocLen*2);
			values=(mpfr_t*)realloc(values, sizeof(mpfr_t)*allocLen);
		}
		assert(len < allocLen);
		mpfr_init2(values[len], mpfr_get_prec(x));
		mpfr_set(values[len], x, MPFR_RNDN);
		++len;
	}
	
	const mpfr_t *begin() const
	{ return values; }
	
	const mpfr_t *end() const
	{ return values+len; }
	
	~MPFRVec()
	{
		clear();
		if(values){
			free(values);
			values=0;
			allocLen=0;
		}
	}
	
	int size() const
	{ return len; }
	
	mpfr_t &operator[](int i)
	{
		if((i<0) || (i>=len))
			throw std::string("MPFRVec - Out of range.");
		return values[i];
	}
	
	const mpfr_t &operator[](int i) const
	{
		if((i<0) || (i>=len))
			throw std::string("MPFRVec - Out of range.");
		return values[i];
	}
};

}; // random
}; // flopoco

#endif
