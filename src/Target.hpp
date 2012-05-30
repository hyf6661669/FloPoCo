/*
  The abstract class that models different chips for delay and area. 
  Classes for real chips inherit from this one. They should be in subdirectory Targets
 
  Authors:   Bogdan Pasca, Florent de Dinechin

  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  

  All Rights Reserved
*/


//TODO a logicTableDelay() that would replace the delay computation in Table.cpp (start from there)


#ifndef TARGET_HPP
#define TARGET_HPP

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <vector>
#include "Targets/DSP.hpp"
#include <string>

using namespace std;

namespace flopoco{


	/** Abstract target Class. All classes which model real chips inherit from this class */
	class Target
	{
	public:
		/** The default constructor. Creates a pipelined target, with 4 inputs/LUT, 
		 * with a desired frequency of 400MHz and which is allowed to use hardware multipliers
		 */ 
		Target()   {
			pipeline_          = true;
			lutInputs_         = 4;
			frequency_         = 400000000.;
			useHardMultipliers_= true;
			id_                = "generic";
		}
	
		/** The destructor */
		virtual ~Target() {}

		/** Returns ID of instantiated target. This ID is represented by the name
		 * @return the ID
		 */
		string getID();

		/** Returns ID of the vendor, currently "Altera" or "Xilinx". 
		 * @return the ID
		 */
		string getVendor(){
			return vendor_;
		}

		// Architecture-related methods
		/** Returns the number of inputs that the LUTs have on the specific device
		 * @return the number of inputs for the look-up tables (LUTs) of the device
		 */
		int lutInputs();
	
		/** Function for determining the submultiplication sizes so that the design is able 
		 * to function at a desired frequency ( multiplicatio is considerd X * Y )
		 * @param[in,out] x the size of the submultiplication for x
		 * @param[in,out] y the size of the submultiplication for y
		 * @param[in] wInX the width of X
		 * @param[in] wInY the width of Y
		 */
		virtual bool suggestSubmultSize(int &x, int &y, int wInX, int wInY)=0;

		/** Function for determining the subadition sizes so that the design is able 
		 * to function at a desired frequency ( addition is considerd X + Y )
		 * @param[in,out] x the size of the subaddition for x and y
		 * @param[in] wIn the widths of X and Y
		 */
		virtual bool suggestSubaddSize(int &x, int wIn)=0; 	

		/** Function for determining the subadition sizes so that the design is able 
		 * to function at a desired frequency ( addition is considerd X + Y )
		 * @param[in,out] x the size of the subaddition for x and y
		 * @param[in] wIn the widths of X and Y
		 * @param[in] slack the time delay consumed out of the input period 
		 */
		virtual bool   suggestSlackSubaddSize(int &x, int wIn, double slack)=0;

		/** Function for determining the subcomparator sizes so that the design is able 
		* to function at a desired frequency 
		* @param[in,out] x the size of the subaddition for x and y
		* @param[in] wIn the widths of x and y
		* @param[in] slack the time delay consumed out of the input period 
		* @param[in] constant if the comparisson is done against a constant or Y
		*/
		virtual bool suggestSlackSubcomparatorSize(int &x, int wIn, double slack, bool constant)=0;
		
		
		/* -------------------- Delay-related methods ------------------------*/

		/* -------------------  logic related  -------------------------------*/
		
		/** Function which returns the lutDelay for this target
		 * @return the LUT delay
		 */
		virtual double lutDelay() =0;
	
		/** Function which returns the flip-flop Delay for this target
			 (not including any net delay)
			 * @return the flip-flop delay
			 */
		virtual double ffDelay() =0;
	
		/** Function which returns the carry propagate delay
		 * @return the carry propagate dealy
		 */
		virtual double carryPropagateDelay() =0;

		/** Function which returns addition delay for an n bit addition
		 * @param n the number of bits of the addition (n-bit + n-bit )
		 * @return the addition delay for an n bit addition
		 */
		virtual double adderDelay(int n) =0;

		/** Function which returns the delay for an n bit equality comparison 
		* @param n the number of bits of the comparisson
		* @return the delay of the comparisson between two vectors
		*/
		virtual double eqComparatorDelay(int n) =0;

		/** Function which returns the delay for an n bit equality comparison with a 
		* constant
		* @param n the number of bits of the comparisson
		* @return the delay of the comparisson between a vector and a constant
		*/
		virtual double eqConstComparatorDelay(int n) =0;
		
		
		/** Function which returns the distant wire delay.
		 * @return distant wire delay 
		 */
		virtual double distantWireDelay(int n) =0;

		/** Function which returns the local wire delay (local routing)
		 * @return local wire delay 
		 */
		virtual double localWireDelay(int fanout = 1) =0;

		/* --------------------  DSP related  --------------------------------*/

		/** Function which returns the delay of the multiplier of the DSP block
		 * @return delay of the DSP multiplier
		 */
		virtual double DSPMultiplierDelay() =0;

		/** Function which returns the delay of the adder of the DSP block
		 * @return delay of the DSP adder
		 */
		virtual double DSPAdderDelay()=0;
		
		/** Function which returns the delay of the wiring between neighboring 
		 * DSPs
		 * @return delay of the interconnect between DSPs
		 */
		virtual double DSPCascadingWireDelay()=0;

		/** Function which returns the delay of the wiring between  
		 * DSPs output and logic (slices)
		 * @return delay DSP -> slice
		 */
		virtual double DSPToLogicWireDelay()=0;

		/** Function which returns the delay of the wiring between 
		 * slices and DSPs input
		 * @return delay slice -> DSP
		 */
		virtual double LogicToDSPWireDelay()=0;

		/* -------------------  BRAM related  --------------------------------*/

		/** Function which returns the delay between
		 * RAMs and slices
		 * @return delay RAMout to slice
		 */
		virtual double RAMToLogicWireDelay()=0;

		/** Function which returns the delay between
		 * slices and RAM input
		 * @return delay slice to RAMin
		 */
		virtual double LogicToRAMWireDelay()=0;

		/** Function which returns the delay of the RAM
		 * @return delay RAM
		 */
		virtual double RAMDelay()=0;


		/** Function which returns the size of a primitive memory block,which could be recognized by the synthesizer as a dual block.
		 * @return the size of the primitive memory block
		 */	
		virtual long sizeOfMemoryBlock() = 0 ;

	
		/** Function which returns the Equivalence between slices and a DSP.
		 * @return X ,  where X * Slices = 1 DSP
		 */
		virtual int getEquivalenceSliceDSP() = 0;
	
		/** Function which returns the number of DSPs that exist in FPGA
		 * @return number of DSPs
		 */
	
		virtual int getNumberOfDSPs() = 0;
	
		/** Function which returns the maximum widths of the operands of a DSP
		 * @return widths
		 */
		virtual void getDSPWidths(int &x, int &y, bool sign = false) = 0;
	
		virtual void getAdderParameters(double &k1, double &k2, int size) = 0;

		// Methods related to target behaviour and performance
		/** Sets the target to pipelined */
		void setPipelined();                
	
		/**< Sets the target to combinatorial */    
		void setNotPipelined();                 
	
		/** Returns true if the target is pipelined, otherwise false
		 * @return if the target is pipelined
		 */
		bool isPipelined();
	
		/** Returns the desired frequency for this target in Hz
		 * @return the frequency
		 */
		double frequency();

		/** Returns the desired frequency for this target in MHz
		 * @return the frequency
		 */
		double frequencyMHz();


		/** Returns the target frequency, normalized between 0 and 1 in a target-independent way.
			 1 means maximum practical frequency (400MHz on Virtex4, 500MHz on Virtex-5, etc)
			 This method is intended to make it easier to write target-independent frequency-directed operators
		 */
		double normalizedFrequency();

	
		/** Sets the desired frequency for this target
		 * @param f the desired frequency
		 */	
		void setFrequency(double f);

		/** Sets the use of hardware multipliers 
		 * @param v use or not harware multipliers
		 */
		void setUseHardMultipliers(bool v);

		/** Returns true if the operator for this target is allowed to use hardware multipliers
		 * @return the status of hardware multipliers usage
		 */
		bool getUseHardMultipliers(); 

		/** Function which returns the number of LUTs needed to implement
		 *	 a multiplier of the given width
		 * @param	wInX the width (in bits) of the first operand
		 * @param	wInY the width (in bits) of the second operand
		 */
		virtual int getIntMultiplierCost(int wInX, int wInY) =0;
	
	
		/** Function which return the number of Slices needed to implement
		 * an adder that has 2 or more operands.
		 * @param wIn the width (in bits) of the adder operands
		 * @param n the number of operands
		 * @return the cost of the adder in Slices/ALMs
		 */
		virtual int getIntNAdderCost(int wIn, int n) =0;
	
		/** Constructs a specific DSP to each target */
		virtual DSP* createDSP() = 0;
	
		
		
		/*-------------- Resource Estimation related items	----------*/
		
		/**
		 * Return the number of LUT elements that should be added to the 
		 * total resource estimation by adding @count LUTs of type @type.
		 * According to the LUT type and number of the LUTs, the final 
		 * number of LUTs that are inferred after synthesis might be 
		 * different from the initial estimations.
		 * @param count (by default 1) the number of LUTs that the user 
		 * wishes to add to the design
		 * @param nrInputs (by default 0, meaning the value returned by 
		 * the @suggestLUTType() function) the number of inputs of the LUT
		 * @return the type of the LUT to be used for resource estimation
		 */ 
		virtual int lutCount(int nrInputs = 0, int count = 1);
		
		/**
		 * Return the number of Multipliers that should be added to the 
		 * total resource count, when the user wishes to add @count 
		 * multipliers that each have inputs of widths @widthX and
		 * @widthY bits respectively.
		 * The function should take into account the DSP that are only 
		 * partially filled, and that might be completed after another 
		 * call to the function.
		 * @param count (by default 1) the number of such multipliers 
		 * that need to be added
		 * @param widthX the width of the first term to the operation
		 * @param widthY the width of the second term to the operation
		 * @return the number of multipliers to add
		 */
		virtual int multiplierCount(int widthX, int widthY, int count = 1);
		
		/**
		 * Return the number of DSP that should be added to the total 
		 * resource count, after adding count multipliers that each have 
		 * inputs of widths @widthX and @widthY bits respectively.
		 * The function should take into account the DSPs that are only 
		 * partially filled, and that might be completed after another 
		 * call to the function.
		 * @param count (by default 1) the number of such multipliers 
		 * that need to be added
		 * @param widthX the width of the first term to the operation
		 * @param widthY the width of the second term to the operation
		 * @return the number of DSP blocks to add
		 */
		virtual int dspForMultiplier(int widthX, int widthY, int count = 1);
		
		/**
		 * Return the number of LUTs that should be added to the total 
		 * resource count, after adding count multipliers that each have 
		 * inputs of widths @widthX and @widthY bits respectively.
		 * @param count (by default 1) the number of such multipliers 
		 * that need to be added
		 * @param widthX the width of the first term to the operation
		 * @param widthY the width of the second term to the operation
		 * @return the number of LUTs to add
		 */
		virtual int lutForMultiplier(int widthX, int widthY, int count = 1);
		
		/**
		 * Return the number of FFs that should be added to the total 
		 * resource count, after adding count multipliers that each have 
		 * inputs of widths @widthX and @widthY bits respectively.
		 * @param count (by default 1) the number of such multipliers 
		 * that need to be added
		 * @param widthX the width of the first term to the operation
		 * @param widthY the width of the second term to the operation
		 * @return the number of FFs to add
		 */
		virtual int ffForMultiplier(int widthX, int widthY, int count = 1);
		
		/**
		 * Return the number of LUTs that should be added to the total 
		 * resource count, after adding @count adders that each have 
		 * inputs of widths @widthX and @widthY bits respectively.
		 * @param count (by default 1) the number of such adders 
		 * that need to be added
		 * @param widthX the width of the first term to the operation
		 * @param widthY the width of the second term to the operation
		 * @return the number of LUTs to add
		 */
		virtual int lutForAdderSubtracter(int widthX, int widthY, int count = 1);
		
		/**
		 * Return the number of FFs that should be added to the total 
		 * resource count, after adding @count adders that each have 
		 * inputs of widths @widthX and @widthY bits respectively.
		 * @param count (by default 1) the number of such adders 
		 * that need to be added
		 * @param widthX the width of the first term to the operation
		 * @param widthY the width of the second term to the operation
		 * @return the number of FFs to add
		 */
		virtual int ffForAdderSubtracter(int widthX, int widthY, int count = 1);
		
		/**
		 * Return the number of Memory elements that should be added to 
		 * the total resource count, when the user wishes to add @count 
		 * memories, each having @size words of width @width. The 
		 * memories can be RAMs or ROMs, depending on the @type parameter
		 * @param count (by default 1) the number of memories with the 
		 * specified characteristics to add
		 * @param size the size of the memory, measured in words
		 * @param width the width of the memory words
		 * @param type (by default 0, RAM) the memory type (RAM or ROM)
		 */
		virtual int memoryCount(int size, int width, int type = 0, int count = 1);
		
		/**
		 * Return the number of Shift Registers that should be added to the 
		 * total resource count, when the user wishes to add @count 
		 * shift registers that each have inputs of width @width bits.
		 * @param count (by default 1) the number of such SRLs that need 
		 * to be added
		 * @param width the width of the shift registers
		 * @param the depth of the shift register
		 * @return the number of SRLs to add
		 */
		virtual int srlCount(int width, int depth, int count = 1);
		
		/**
		 * Return the number of LUTs that should be added to the 
		 * total resource count when the user wishes to add @count 
		 * shift registers each having the input of width @width bits.
		 * @param count (by default 1) the number of such SRLs that need 
		 * to be added
		 * @param width the width of the shift register(s)
		 * @param depth the depth of the shift register(s)
		 * @return the number of LUTs to add to the resource estimation
		 */
		virtual int lutForSRL(int width, int depth, int count = 1);
		
		/**
		 * Return the number of FFs that should be added to the 
		 * total resource count, when the user wishes to add @count 
		 * shift registers that each have inputs of width @width bits.
		 * @param count (by default 1) the number of such SRLs that need 
		 * to be added
		 * @param width the width of the shift register(s)
		 * @param depth the depth of the shift register(s)
		 * @return the number of FFs to add to the resource estimation
		 */
		virtual int ffForSRL(int width, int depth, int count = 1);
		
		/**
		 * Return the number of RAM blocks that should be added to the 
		 * total resource count, when the user wishes to add @count 
		 * shift registers that each have inputs of width @width bits.
		 * @param count (by default 1) the number of such SRLs that need 
		 * to be added
		 * @param width the width of the shift register(s)
		 * @param depth the depth of the shift register(s)
		 * @return the number of RAM blocks to add to the resource 
		 * estimation
		 */
		virtual int ramForSRL(int width, int depth, int count = 1);
		
		/**
		 * Return the number of Multiplexers that should be added to the 
		 * total resource count, when the user wishes to add @count 
		 * multiplexers that each have inputs of width @width bits.
		 * @param count (by default 1) the number of such Muxs that need 
		 * to be added
		 * @param nrInputs the number of inputs of the multiplexer
		 * @param width the width of the multiplexers
		 * @return the number of Multiplexers to add
		 */
		virtual int muxCount(int nrInputs, int width, int count = 1);
		
		/**
		 * Return the number of LUTs that should be added to the 
		 * total resource count, when the user wishes to add @count 
		 * multiplexers that each have inputs of width @width bits.
		 * @param count (by default 1) the number of such Multiplexers 
		 * that need to be added
		 * @param nrInputs the number of inputs of the multiplexer
		 * @param width the width of the multiplexers
		 * @return the number of LUTs to add to the resource estimation
		 */
		virtual int lutForMux(int nrInputs, int width, int count = 1);
		
		/**
		 * Return the number of LUTs that should be added to the 
		 * total resource count, when the user wishes to add @count 
		 * counters that each have inputs of width @width bits.
		 * @param count (by default 1) the number of such counters
		 * that need to be added
		 * @param width the width of the counters
		 * @return the number of LUTs to add to the resource estimation
		 */
		virtual int lutForCounter(int width, int count = 1);
		
		/**
		 * Return the number of FFs that should be added to the 
		 * total resource count, when the user wishes to add @count 
		 * counters that each have inputs of width @width bits.
		 * @param count (by default 1) the number of such counters
		 * that need to be added
		 * @param width the width of the counters
		 * @return the number of FFs to add to the resource estimation
		 */
		virtual int ffForCounter(int width, int count = 1);
		
		/**
		 * Return the number of LUTs that should be added to the 
		 * total resource count, when the user wishes to add @count 
		 * accumulators that each have inputs of width @width bits.
		 * @param count (by default 1) the number of such accumulators
		 * that need to be added
		 * @param width the width of the accumulators
		 * @param useDSP (by default false) whether DSPs are used, or not
		 * @return the number of LUTs to add to the resource estimation
		 */
		virtual int lutForAccumulator(int width, bool useDSP = false, int count = 1);
		
		/**
		 * Return the number of FFs that should be added to the 
		 * total resource count, when the user wishes to add @count 
		 * accumulators that each have inputs of width @width bits.
		 * @param count (by default 1) the number of such accumulators
		 * that need to be added
		 * @param width the width of the accumulators
		 * @param useDSP whether DSPs are used, or not
		 * @return the number of FFs to add to the resource estimation
		 */
		virtual int ffForAccumulator(int width, bool useDSP = false, int count = 1);
		
		/**
		 * Return the number of DSPs that should be added to the 
		 * total resource count, when the user wishes to add @count 
		 * accumulators that each have inputs of width @width bits.
		 * @param count (by default 1) the number of such accumulators
		 * that need to be added
		 * @param width the width of the accumulators
		 * @param useDSP whether the use of DSPs is allowed
		 * @return the number of DSPs to add to the resource estimation
		 */
		virtual int dspForAccumulator(int width, bool useDSP, int count = 1);
		
		/**
		 * Return the number of LUTs that should be added to the 
		 * total resource count, when the user wishes to add @count 
		 * encoders/decoders that each have inputs of width @width bits.
		 * @param count (by default 1) the number of such encoders/decoders
		 * that need to be added
		 * @param width the width of the encoders/decoders
		 * @return the number of LUTs to add to the resource estimation
		 */
		virtual int lutForDecoder(int width, int count = 1);
		
		/**
		 * Return the number of FFs that should be added to the 
		 * total resource count, when the user wishes to add @count 
		 * encoders/decoders that each have inputs of width @width bits.
		 * @param count (by default 1) the number of such encoders/decoders
		 * that need to be added
		 * @param width the width of the encoders/decoders
		 * @return the number of FFs to add to the resource estimation
		 */
		virtual int ffForDecoder(int width, int count = 1);
		
		/**
		 * Return the number of LUTs that should be added to the 
		 * total resource count, when the user wishes to add @count 
		 * arithmetic operators that each have inputs of width @width 
		 * bits.
		 * @param count (by default 1) the number of such arithmetic 
		 * operators that need to be added
		 * @param nrInputs the number of inputs of the arithmetic 
		 * operator
		 * @param width the width of the arithmetic operators' inputs
		 * @return the number of LUTs to add to the resource estimation
		 */
		virtual int lutForArithmeticOperator(int nrInputs, int width, int count = 1);
		
		/**
		 * Return the number of LUTs that should be added to the 
		 * total resource count, when the user wishes to add @count 
		 * FSMs.
		 * @param count (by default 1) the number of such FSM that need 
		 * to be added
		 * @param nrStates the number of states of the finite state machine
		 * @param nrTransitions the number of transitions of the state 
		 * machine
		 * @return the number of LUTs to add to the resource estimation
		 */
		virtual int lutForFSM(int nrStates, int nrTransitions, int count = 1);
		
		/**
		 * Return the number of FFs that should be added to the 
		 * total resource count, when the user wishes to add @count 
		 * FSMs.
		 * @param count (by default 1) the number of such FSM that need 
		 * to be added
		 * @param nrStates the number of states of the finite state machine
		 * @param nrTransitions the number of transitions of the state 
		 * machine
		 * @return the number of LUTs to add to the resource estimation
		 */
		virtual int ffForFSM(int nrStates, int nrTransitions, int count = 1);
		
		/**
		 * Return the number of RAMs that should be added to the 
		 * total resource count, when the user wishes to add @count 
		 * FSM that each have inputs of width @width bits.
		 * @param count (by default 1) the number of such FSM
		 * that need to be added
		 * @param states the number of states of the finite state machine
		 * @param transitions the number of transitions of the state machine
		 * @param width the width of the FSM
		 * @return the number of RAMs to add to the resource estimation
		 */
		virtual int ramForFSM(int nrStates, int nrTransitions, int count = 1);
		
		/*------------ Resource Estimation - target specific ---------*/
		
		
		/**
		 * IMPORTANT: These functions return values that are specific 
		 * to the target FPGA. This is an implementation that tries to 
		 * give as good as an approximation as possible, but if the 
		 * results are not satisfactory, the class that implements Target
		 * can just as well override the functions to provide better 
		 * ones.
		 */
		
		
		/**
		 * Determine whether a LUT can be split to implement multiple 
		 * independent functions. The decision is taken based on the type 
		 * of LUT being added (currently the type is the number of inputs 
		 * of the LUT).
		 * @param nrInputs the number of inputs of the LUT being added
		 * @return whether the LUTs in the specific architecture can/cannot 
		 * be split to perform multiple independent functions
		 */
		virtual bool lutSplitInputs(int nrInputs);
		
		/**
		 * Determine whether the DSP can be used to implement multiple
		 * independent multiplications. The decision is taken based on 
		 * the target FPGA's characteristics.
		 * @return whether or not the specific architecture can be use a
		 * DSP block to implement several independent multiplications
		 */
		virtual bool dspSplitMult();
		
		/**
		 * The DSP can house several independent multiplications, so this
		 * function determines the required number of DSPs for creating
		 * @count multipliers having the inputs of widths @widthX and
		 * @widthY, respectively.
		 * @param widthX the width of the first input to the multiplier
		 * @param widthY the width of the second input to the multiplier
		 * @return the number of required DSP blocks
		 */
		virtual double getMultPerDSP(int widthX, int widthY);
		
		/**
		 * Determine the required number of LUTs for creating
		 * @count multipliers having the inputs of widths @widthX and
		 * @widthY, respectively.
		 * @param widthX the width of the first input to the multiplier
		 * @param widthY the width of the second input to the multiplier
		 * @return the number of required LUTs
		 */
		virtual double getLUTPerMultiplier(int widthX, int widthY);
		
		/**
		 * Determine the required number of FFs for creating
		 * @count multipliers having the inputs of widths @widthX and
		 * @widthY, respectively.
		 * @param widthX the width of the first input to the multiplier
		 * @param widthY the width of the second input to the multiplier
		 * @return the number of required FFs
		 */
		virtual double getFFPerMultiplier(int widthX, int widthY);
		
		/**
		 * Determine the required number of LUTs for creating @count
		 * adders/subtracters having the inputs of widths @widthX and
		 * @widthY, respectively.
		 * @param widthX the width of the first input to the adder
		 * @param widthY the width of the second input to the adder
		 * @return the number of required LUTs
		 */
		virtual double getLUTPerAdderSubtracter(int widthX, int widthY);
		
		/**
		 * Determine the required number of FFs for creating @count
		 * adders/subtracters having the inputs of widths @widthX and
		 * @widthY, respectively.
		 * @param widthX the width of the first input to the adder
		 * @param widthY the width of the second input to the adder
		 * @return the number of required FFs
		 */
		virtual double getFFPerAdderSubtracter(int widthX, int widthY);
		
		/**
		 * Determine the number of words in a memory block that has the
		 * width given by the @width parameter. Memory blocks can have 
		 * multiple configurations, according to the required width and 
		 * depth. This function is useful for further determining the 
		 * necessary memory blocks.
		 * @param width the width of the memory words in bits
		 * @return the standard number of words per block
		 */
		virtual int wordsPerBlock(int width);
		
		/**
		 * Determine the depth of the shift register, according to the 
		 * required depth, given by the @depth parameter. Shift registers 
		 * can be configured to different depths on the target FPGA, thus 
		 * allowing the use of several SRLs with less resources.
		 * @param depth the required depth of the shift registers
		 * @return the default depth for a shift register, best suiting 
		 * a shift register of @depth depth
		 */
		virtual int getSRLDepth(int depth);
		
		/**
		 * Determine the required number of LUTs for a shift register 
		 * having the depth given by the @depth parameter. The number of
		 * LUTs depends on the target FPGA. A return value of 0 symbolizes
		 * that no LUTs are used to build a SRL of the given depth on the
		 * target FPGA.
		 * @param depth the required depth of the shift register
		 * @return the number of LUTs required for creating a shift 
		 * register of depth @depth; return value of 0 means that the 
		 * resource is not required
		 */
		virtual double getLUTPerSRL(int depth);
		
		/**
		 * Determine the required number of FFs for a shift register 
		 * having the depth given by the @depth parameter. The number of
		 * FFs depends on the target FPGA. A return value of 0 symbolizes
		 * that no FFs are used to build a SRL of the given depth on the
		 * target FPGA.
		 * @param depth the required depth of the shift register
		 * @return the number of FFs required for creating a shift 
		 * register of depth @depth; return value of 0 means that the 
		 * resource is not required
		 */
		virtual double getFFPerSRL(int depth);
		
		/**
		 * Determine the required number of RAM elements for a shift 
		 * register having the depth given by the @depth parameter. The 
		 * number of RAM elements depends on the target FPGA. A return 
		 * value of 0 symbolizes that no RAM elements are used to build 
		 * a SRL of the given depth on the target FPGA.
		 * @param depth the required depth of the shift register
		 * @return the number of RAM elements required for creating a 
		 * shift register of depth @depth; return value of 0 means that 
		 * the resource is not required
		 */
		virtual double getRAMPerSRL(int depth);
		
		/**
		 * Determine the required number of LUTs for a multiplexer having 
		 * @nrInputs inputs. The number of LUTs depends on the target
		 * FPGA (the function generator architecture and other available 
		 * resources on the chip might influence the resource count).
		 */
		virtual double getLUTFromMux(int nrInputs);
		
		/**
		 * Determine the required number of LUTs for a counter having a
		 * bitwidth of @width bits. The number of LUTs can vary according 
		 * to the properties of the FPGA and the width of the counter.
		 * This function computes a factor that symbolizes the ratio of
		 * LUTs per counter bits.
		 * @param width the width of the counter
		 * @return the ratio of LUTs per counter bits
		 */
		virtual double getLUTPerCounter(int width);
		
		/**
		 * Determine the required number of FFs for a counter having a
		 * bitwidth of @width bits. The number of FFs can vary according 
		 * to the properties of the FPGA and the width of the counter.
		 * This function computes a factor that symbolizes the ratio of
		 * FFs per counter bits.
		 * @param width the width of the counter
		 * @return the ratio of LUTs per counter bits
		 */
		virtual double getFFPerCounter(int width);
		
		/**
		 * Determine the required number of LUTs for an accumulator having a
		 * bitwidth of @width bits. The number of LUTs can vary according 
		 * to the properties of the FPGA, the width of the accumulator and 
		 * the use of DSP blocks.
		 * This function computes a factor that symbolizes the ratio of
		 * LUTs per accumulator bits.
		 * @param width the width of the counter
		 * @return the ratio of LUTs per counter bits
		 */
		virtual double getLUTPerAccumulator(int width, bool useDSP);
		
		/**
		 * Determine the required number of FFs for an accumulator having a
		 * bitwidth of @width bits. The number of FFs can vary according 
		 * to the properties of the FPGA, the width of the accumulator and 
		 * the use of DSP blocks.
		 * This function computes a factor that symbolizes the ratio of
		 * FFs per accumulator bits.
		 * @param width the width of the counter
		 * @return the ratio of FFs per counter bits
		 */
		virtual double getFFPerAccumulator(int width, bool useDSP);
		
		/**
		 * Determine the required number of DSPs for an accumulator having a
		 * bitwidth of @width bits. The number of DSPs can vary according 
		 * to the properties of the FPGA and the width of the accumulator.
		 * This function computes a factor that symbolizes the ratio of
		 * DSPs per accumulator bits.
		 * @param width the width of the counter
		 * @return the ratio of DSPs per counter bits
		 */
		virtual double getDSPPerAccumulator(int width);
		
		/**
		 * Determine the required number of LUTs for a decoder having a
		 * bitwidth of @width bits. The number of LUTs can vary according 
		 * to the properties of the FPGA, the width of the decoder.
		 * This function computes a factor that symbolizes the ratio of
		 * LUTs per decoder bits.
		 * @param width the width of the decoder
		 * @return the ratio of LUTs per decoder bits
		 */
		virtual double getLUTPerDecoder(int width);
		
		/**
		 * Determine the required number of FFs for a decoder having a
		 * bitwidth of @width bits. The number of FFs can vary according 
		 * to the properties of the FPGA, the width of the decoder.
		 * This function computes a factor that symbolizes the ratio of
		 * FFs per decoder bits.
		 * @param width the width of the decoder
		 * @return the ratio of FFs per decoder bits
		 */
		virtual double getFFPerDecoder(int width);
		/*------------------------------------------------------------*/
		
		
		/*------------ Floorplanning Related Items -------------------*/
		/**
		 * NOTE: These variables should be set for each different FPGA 
		 * architecture, in their corresponding constructor.
		 */
		vector<int> multiplierPosition;			/**< The position of the columns of multipliers. The Position represents the neighboring LUT column, on the left. */
		vector<int> memoryPosition;				/**< The position of the columns of memories. The Position represents the neighboring LUT column, on the left. */
		int topSliceX;							/**< The x coordinate of the top right slice. */
		int topSliceY;							/**< The y coordinate of the top right slice. */
		int lutPerSlice;						/**< The number of function generators per slice. */
		int ffPerSlice;							/**< The number of registers per slice. */
		int dspHeightInLUT;						/**< The height of a DSP cell, expressed using the height of one LUT as unit of measure */
		int ramHeightInLUT;						/**< The height of a RAM block, expressed using the height of one LUT as unit of measure */
		int dspPerColumn;						/**< The number of DSP blocks in a column of DSPs */
		int ramPerColumn;						/**< The number of RAM blocks in a column of RAMs */
		/*------------------------------------------------------------*/
		
		
		//todo

	protected:
		string id_;
		string vendor_;
		int    lutInputs_;          /**< The number of inputs for the LUTs */
		bool   pipeline_;           /**< True if the target is pipelined/ false otherwise */
		double frequency_;          /**< The desired frequency for the operator in Hz */
		int    multXInputs_;        /**< The size for the X dimmension of the hardware multipliers */
		int    multYInputs_;        /**< The size for the Y dimmension of the hardware multipliers */
		bool   useHardMultipliers_; /**< If true the operator design can use the hardware multipliers */
		long   sizeOfBlock_;		/**< The size of a primitive memory block */
		double maxFrequencyMHz_ ;/**< The maximum practical frequency attainable on this target. An indicator of relative performance of FPGAs. 400 is for Virtex4 */
		
	
	};

}
#endif
