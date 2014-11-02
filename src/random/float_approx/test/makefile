
FLOPOCO = ../../../../cygwin/flopoco

fa_sincos_w16_deg3.vhdl :
	$(FLOPOCO) -outputfile=$@ -verbose=2 FloatApprox 8 16  1 8  8 16  "cos(x)-sin(x*1.3)+1.1" 3

fa_sincos_w24_deg3.vhdl :
	$(FLOPOCO) -outputfile=$@ -verbose=2 FloatApprox 8 24  1 8  8 24  "cos(x)-sin(x*1.3)+1.1" 3
	
# Loop in PolyEval
#fa_sincos_w16_deg2.vhdl :
#	$(FLOPOCO) -outputfile=$@ -verbose=2 FloatApprox 8 16  1 8  8 16  "cos(x)-sin(x*1.3)+1.1" 2

fa_sincos_w24_deg2.vhdl :
	$(FLOPOCO) -outputfile=$@ -verbose=2 FloatApprox 8 23  1 8  8 23  "cos(x)-sin(x*1.3)+1.1" 2
	

# Loop in PolyEval
#fa_sqrtlog_w16_deg3.vhdl :
#	$(FLOPOCO) -outputfile=$@ -verbose=2 FloatApprox 8 16  "2^-32" "1-2^-16"  8 16  "sqrt(-2*log(x))" 3

fa_sqrtlog_w23_deg3.vhdl :
	$(FLOPOCO) -outputfile=$@ -verbose=2 FloatApprox 8 23  "2^-32" "1-2^-23"  8 23  "sqrt(-2*log(x))" 3
	
fa_erf_w23_deg2.vhdl :
	$(FLOPOCO) -outputfile=$@ -verbose=2 FloatApprox 8 23  "2^-32" "4"  8 23  "erf(x)" 2

all : fa_sincos_w16_deg3.vhdl fa_sincos_w24_deg3.vhdl fa_sincos_w24_deg2.vhdl fa_sincos_w53_deg4.vhdl \
	fa_sqrtlog_w23_deg3.vhdl
