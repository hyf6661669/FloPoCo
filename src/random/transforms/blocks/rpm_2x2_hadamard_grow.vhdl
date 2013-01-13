library ieee; 
use ieee.std_logic_1164.all; 
use ieee.std_logic_arith.all;
use ieee.numeric_std.all; 

library unisim;
use unisim.vcomponents.all; 

entity rpm_2x2_hadamard_grow is
	generic (W:integer:=32; X:integer:=0; Y:integer:=0; SLICE_HEIGHT:integer:=4);
	port(
		clk : in std_logic;
		left_in : in std_logic_vector (W-1 downto 0);
		right_in : in std_logic_vector (W-1 downto 0);
		left_out : out std_logic_vector(W downto 0);
		right_out : out std_logic_vector(W downto 0)
	);
end entity;

architecture RTL_xilinx of rpm_2x2_hadamard_grow is
--	component MUXCY 
--		port(O : out std_logic; DI : in std_logic; CI : in std_logic; S : in std_logic); 
--	end component; 
	
--	component XORCY 
--		port(O : out std_logic; LI : in std_logic; CI : in std_logic); 
--	end component; 
	
--	component LUT2 
--		generic (INIT : bit_vector := "0110"); 
--		port(O : out std_logic; I0 : in std_logic; I1 : in std_logic); 
--	end component; 
	
--	component FD
--		port(C : in std_logic; D : in std_logic; Q : out std_logic); 
--	end component; 

	attribute RLOC : string;
	
	signal left_in_d1, right_in_d1 : std_logic_vector(W downto 0);
	signal carryin_left, carryin_right : std_logic_vector(W+1 downto 0);
	signal carryout_left, carryout_right : std_logic_vector(W downto 0);
	signal sumsig_left, sumsig_right : std_logic_vector(W downto 0);
	signal slct_left, slct_right : std_logic_vector(W downto 0);
begin

	carryin_left(0) <=  '0';
	carryin_right(0) <=  '1';	-- Extra carry in for twos complement of a

	gen_fadds : for i in 0 to W generate 	
		attribute RLOC of fadd_left_in_FD: label is "X" & integer'image(integer(X)) & "Y" & integer'image(integer(i/SLICE_HEIGHT+Y)); 
		attribute RLOC of fadd_left_LUT2: label is "X" & integer'image(integer(X)) & "Y" & integer'image(integer(i/SLICE_HEIGHT+Y)); 
		attribute RLOC of fadd_left_MUXCY: label is "X" & integer'image(integer(X)) & "Y" & integer'image(integer(i/SLICE_HEIGHT+Y)); 
		attribute RLOC of fadd_left_XORCY: label is "X" & integer'image(integer(X)) & "Y" & integer'image(integer(i/SLICE_HEIGHT+Y)); 
		attribute RLOC of fadd_left_out_FD: label is "X" & integer'image(integer(X)) & "Y" & integer'image(integer(i/SLICE_HEIGHT+Y)); 
		
		attribute RLOC of fadd_right_in_FD: label is "X" & integer'image(integer(X+1)) & "Y" & integer'image(integer(i/SLICE_HEIGHT+Y)); 
		attribute RLOC of fadd_right_LUT2: label is "X" & integer'image(integer(X+1)) & "Y" & integer'image(integer(i/SLICE_HEIGHT+Y)); 
		attribute RLOC of fadd_right_MUXCY: label is "X" & integer'image(integer(X+1)) & "Y" & integer'image(integer(i/SLICE_HEIGHT+Y)); 
		attribute RLOC of fadd_right_XORCY: label is "X" & integer'image(integer(X+1)) & "Y" & integer'image(integer(i/SLICE_HEIGHT+Y)); 
		attribute RLOC of fadd_right_out_FD: label is "X" & integer'image(integer(X+1)) & "Y" & integer'image(integer(i/SLICE_HEIGHT+Y)); 
		
		signal src_left, src_right : std_logic;
	begin
		bb : if i=W generate
			src_left <= left_in(W-1);
			src_right <= right_in(W-1);
		end generate;
		aa : if i<W generate
			src_left <= left_in(i);
			src_right <= right_in(i);
		end generate;
	
		fadd_left_in_FD : FD port map (C => clk, Q => left_in_d1(i), D => src_left  );
		fadd_right_in_FD : FD port map (C => clk, Q => right_in_d1(i), D => src_right );
	
		fadd_left_LUT2 : LUT2  generic map (INIT => "0110") port map (O => slct_left(i), I0 => left_in_d1(i), I1 => right_in_d1(i)); 
		-- Original = (a^b), now we want (~a ^ b)
		fadd_right_LUT2 : LUT2  generic map (INIT => "1001") port map (O => slct_right(i), I0 => left_in_d1(i), I1 => right_in_d1(i)); 
		
		fadd_left_MUXCY : MUXCY port map (O => carryout_left(i), DI => left_in_d1(i), CI => carryin_left(i), S => slct_left(i)); 
		fadd_right_MUXCY : MUXCY port map (O => carryout_right(i), DI => left_in_d1(i), CI => carryin_right(i), S => slct_right(i)); 
		
		fadd_left_XORCY : XORCY port map (O => sumsig_left(i), LI => slct_left(i), CI => carryin_left(i)); 
		fadd_right_XORCY : XORCY port map (O => sumsig_right(i), LI => slct_right(i), CI => carryin_right(i)); 
		
		carryin_left(i+1) <= carryout_left(i); 
		carryin_right(i+1) <= carryout_right(i); 
		
		fadd_left_out_FD : FD port map (C => clk, Q => left_out(i), D => sumsig_left(i) );
		fadd_right_out_FD : FD port map (C => clk, Q => right_out(i), D => sumsig_right(i) );
	end generate;
end RTL_xilinx;

--architecture RTL_generic of rpm_adder is
--	signal tmpSum : std_logic_vector(W downto 0);
--begin
--	process (clk) 
--	begin
		--if( rising_edge(clk)) then
			--tmpSum := std_logic_vector(signed('0' & a) + signed('0' & b));
			--tmpSum <= signed('0' & a) + signed('0' & b);
		--end if;
	--end process;
	--sum<=tmpSum(W-1 downto 0);
	--cout<=tmpSum(W);
--end RTL_generic;

