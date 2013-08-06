#ifndef flopoco_random_fixed_point_exp_residual_type_hpp
#define flopoco_random_fixed_point_exp_residual_type_hpp

#include <stdint.h>
#include <cmath>
#include <climits>
#include <cfloat>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <cassert>

#include "result_type.hpp"

namespace flopoco
{
namespace random
{

namespace detail
{	
	template<class T>
	struct iterate_range_t
	{
		T current, delta;
		
		iterate_range_t(T _current, T _delta)
			: current(_current)
			, delta(_delta)
		{}
		
		bool operator==(const iterate_range_t &o) const
		{ return current==o.current; }
		
		bool operator!=(const iterate_range_t &o) const
		{ return current!=o.current; }
		
		iterate_range_t &operator++()
		{ current+=delta; return *this; }
		
		iterate_range_t &operator--()
		{ current-=delta; return *this; }
		
		const T &operator*() const
		{ return current; }
	};
}; //detail
	
	template<class T>
	struct residual_type
	{
	private:
		bool hasSign;
		int width;	
		int shift;
		T rangeMin;
		T rangeMax;
	
		T power2(int x) const
		{ return ldexp((T)1, x); }
		
		T random() const
		{
			T acc=0;
			for(int i=0;i<width;i++){
				acc=acc*2+(rand()>(RAND_MAX/2) ? 1 : 0);
			}
			return ldexp(acc, -width);
		}
	public:	
		residual_type(bool _hasSign, int _width, int _shift=0)
			: hasSign(_hasSign)
			, width(_width)
			, shift(_shift)
			, rangeMin(DBL_MAX)
			, rangeMax(-DBL_MAX)
		{}
			
		bool operator==(const residual_type &rt) const
		{
			return (hasSign==rt.hasSign) && (width==rt.width) && (shift==rt.shift) &&
				(rangeMin==rt.rangeMin) && (rangeMax==rt.rangeMax);
		}
		
		std::string DescriptionId() const
		{
			std::stringstream dst;
			if(hasSign){
				dst<<"SFix_";
			}else{
				dst<<"UFix_";
			}
			dst<<"w"<<width;
			if(shift<0){
				dst<<"_sN"<<-shift;
			}else if(shift>0){
				dst<<"_sP"<<shift;
			}
			if(DBL_MAX<DBL_MIN){
				dst<<"_min"<<ToBits(RangeMin())<<"_max="<<ToBits(RangeMax());
			}
			return dst.str();
		}
	
		// msb==lsb -> single bit wide, msb==lsb-1 -> zero bits
		int Width() const
		{ return NonSignWidth() + (hasSign?1:0); }
		
		int NonSignWidth() const
		{ return width; }
		
		int Shift() const
		{ return shift; }
		
		bool HasSign() const
		{ return hasSign; }
		
		residual_type take_msbs(int x) const
		{
			if(x<0)
				throw std::invalid_argument("Can't take negative number of bits.");
			if(x>Width())
				throw std::invalid_argument("Can't take more than Width() bits.");
			if(x==0)
				return residual_type(false, 0, shift);
			if(x==Width())
				return *this;
			residual_type res=residual_type(hasSign, hasSign?x-1:x, shift);
			if(rangeMin<=rangeMax){
				res.rangeMin=Floor(rangeMin);
				res.rangeMax=Floor(rangeMax);
			}
			return res;
		}
		
		residual_type drop_msbs(int x) const
		{
			if(x<0)
				throw std::invalid_argument("Can't take negative number of bits.");
			if(x>Width())
				throw std::invalid_argument("Can't take more than Width() bits.");
			if(x==0)
				return *this;
			
			residual_type res=residual_type(false, hasSign?width-x+1:width-x, hasSign?shift-x+1:shift-x);
			
			if(rangeMin<=rangeMax){
				res.rangeMin=std::max(rangeMin, res.TypeMin());
				res.rangeMax=std::min(rangeMax, res.TypeMax());
			}
			return res;
		}
		
		int MsbPower() const
		{ return -1+shift; }
		
		int LsbPower() const
		{ return -1+shift-(width-1); }
		
		T TypeMin() const
		{ return hasSign ? -power2(MsbPower()+1) : (T)0; }
		
		T TypeMax() const
		{ return power2(MsbPower()+1)-power2(LsbPower()); }
		
		T TypeDelta() const
		{ return power2(LsbPower()); }
		
		T RangeMin() const
		{ return (rangeMin<=rangeMax) ? rangeMin : TypeMin(); }
		
		T RangeMax() const
		{ return (rangeMin<=rangeMax) ? rangeMax : TypeMax(); }
		
		unsigned RangeSize() const
		{
			T tmp = ((rangeMin<=rangeMax) ? (RangeMax()-RangeMin()) : (TypeMax()-TypeMin())) / TypeDelta() + 1;
			if(tmp>(1u<<31))
				throw std::logic_error("residual_type::RangeSize - Size is more than 2^31.");
			return (unsigned)tmp;
		}
			
		bool IsInRange(const T &x) const
		{ return (RangeMin()<=x) && (x<=RangeMax()); }
		
		T Round(const T &x) const
		{ return ldexp(round(ldexp(x,width-shift)),-width+shift); }
		
		T Floor(const T &x) const
		{ return ldexp(floor(ldexp(x,width-shift)),-width+shift); }
		
		void ResetRange()
		{
			rangeMin=DBL_MAX;
			rangeMax=-DBL_MAX;
		}
		
		void Add(const T &x)
		{
			assert(TypeMin()<=x && x<=TypeMax());
			rangeMin=std::min(rangeMin, x);
			rangeMax=std::max(rangeMax, x);
		}
		
		T RandomElement() const
		{
			T tmp=random();
			tmp=RangeMin()+(RangeMax()+TypeDelta()-RangeMin())*tmp;
			return Floor(tmp);
		}
		
		uint64_t ToBits(const T &x) const
		{
			if(Width()>63)
				throw std::runtime_error("Type is wider than 63 bits.");
			if(!IsInRange(x))
				throw std::runtime_error("Value is out of range.");
			T ob=ldexp(x,width-shift);
			T b=round(ob);
			if(b!=ob)
				throw std::runtime_error("Value is not rounded.");
			int64_t ix=b;
			assert(ix==(T)b);
			uint64_t mask=(1ULL<<Width())-1;
			return mask&(uint64_t)ix;
		}
		
		mpz_class ToMpz(const T &x) const
		{ return to_mpz_class(ToBits(x)); }
		
		T FromBits(uint64_t x) const
		{
			if(Width()>63)
				throw std::runtime_error("Type is wider than 63 bits.");
			int64_t sx=(int64_t)x;
			if(hasSign){
				sx=(sx<<(63-width))>>(63-width);
			}
			T b=(T)(sx);
			T xx=ldexp(b, -width+shift);
			if(!IsInRange(xx))
				throw std::runtime_error("Value is out of range.");
			return xx;
		}
		
		T FromMpz(const mpz_class &x) const
		{ return FromBits(to_uint64_t(x)); }
		
		typedef detail::iterate_range_t<T> iterator;
		
		iterator begin_msbs(int msbs) const
		{
			assert(msbs <= Width());
			return iterator(TypeMin(), TypeDelta());
		}
			
		iterator end_msbs(int msbs) const
		{
			assert(msbs <= Width());
			return iterator(TypeMax()+TypeDelta(), TypeDelta());
		}
		
		iterator begin() const
		{ return begin_msbs(Width()); }
		
		iterator end() const
		{ return end_msbs(Width()); }
	};
	
	template<class T>
	std::ostream &operator<<(std::ostream &dst, const residual_type<T> &t)
	{
		if(t.HasSign()){
			dst<<"SFix[";
		}else{
			dst<<"UFix[";
		}
		dst<<"w="<<t.NonSignWidth()<<";s="<<t.Shift();
		if(DBL_MAX<DBL_MIN){
			dst<<";min="<<t.RangeMin()<<";max="<<t.RangeMax();
		}
		dst<<"]";
		return dst;
	}
	
}; // random
}; // flopoco

#endif
