#ifndef HLSOperator_hpp
#define HLSOperator_hpp

#include "Operator.hpp"

#include "hls/HLSTypes.hpp"
#include "hls/HLSExpr.hpp"

#include <set>

namespace flopoco
{
    class HLSContext;
    class HLSScope;

    class HLSOperator
    {
    private:
        HLSOperator &operator=(const HLSOperator &);
    
    
        // Used to remember if this operator is bound to a context on the
        // current thread.

    protected:
    
        HLSOperator()
        {}  
    public:       
        virtual ~HLSOperator()
        {}

        //! Gets a pointer to the associated operator
        /*! Often this will simply re-cast to the same object as an Operator, but
            there is no requirement that they be the same object. However, the
            object returned should always be the same object.
    
            Any HLSOperator which returns the same object from getOperator must
            be functionally equivalent (i.e. they must be the object, or be immutable)
        */
        virtual const Operator &getOperator() const =0;

      //! Stop the automatic walking of Operators that are used by this one
      /*! Sometimes things like IntMultiplier are used in the VHDL route,
	but it doesn't make sense to output them for HLS */
      virtual std::set<Operator*> suppressedHLSDefinitions() const
      { return std::set<Operator*>(); }
    
        //! Responsible for the operator-specific part of the HLS output. Must be overriden
        virtual void emitHLSBody
        (
            HLSContext &ctxt,
            HLSScope &scope
        )const =0;
        
        //! This should be called once for each call of getHLSOperator
        /*! This due to the general object lifetime problem of FloPoCo,
            where ownership is difficult to deal with. For a proxy object
            this will delete the proxy, but for inherited support the
            lifetime needs to be managed on the main Operator
        */
        virtual void releaseHLS() =0;
    };
    
    //! Try to get a HLSOperator associated with the given operator
    /*! This may simply re-cast the object if it natively supports HLSOperator,
        or may create some sort of proxy. Each call to this function may return
        a different object for the same op, but all should be functionally equivalent.
    
        Each object returned from here should be freed using releaseHLS.

	Any object returned from here should not outlive the Operator passed as op.
    */
    HLSOperator *getHLSOperator(const Operator *op);
    

    /* The following functions are bound up with HLSScope */

	//! Declare a new variable with the given type
	HLSExpr hls_declare(const std::string &name, const HLSType &type);

	//! Declares a new variable with unsigned type of given width
	HLSExpr hls_declare(const std::string &name, int w);

	//! Get some kind of existing declaration (variable, input, or output)
	HLSExpr hls_get(const std::string &name);

}; // flopoco
    
#endif

#include "HLSContext.hpp"
