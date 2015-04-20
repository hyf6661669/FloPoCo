#include "hls/HLSTypes.hpp"

namespace flopoco
{
  HLSTypePtr makeHLSType(const Signal &sig)
  {
        if(sig.isFP()){
            return HLSTypeFloat::create(sig.wE(), sig.wF());
        }else if( sig.isBus()){
            return HLSTypeInt::create(false, sig.width());
        }else if(sig.width()==1 && !sig.isBus()){
            return HLSTypeBool::create();
        }else{
            throw std::runtime_error("Unknown signal type for signal "+sig.getName());
        }
  }


};
