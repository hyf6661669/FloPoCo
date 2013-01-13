library ieee; 
use ieee.std_logic_1164.all; 
use ieee.std_logic_arith.all;
use ieee.numeric_std.all; 

library unisim;
use unisim.vcomponents.all; 

entity rpm_nxn_hadamard_grow is
	generic (LOG2N:integer:=3; W:integer:=32; SLICE_HEIGHT:integer:=4);
	port(
		clk : in std_logic;
		src : in std_logic_vector (W*(2**LOG2N)-1 downto 0);
		dst : out std_logic_vector((W+LOG2N)*(2**LOG2N)-1 downto 0)
	);
end entity;

architecture RTL_xilinx of rpm_nxn_hadamard_grow is
begin

	base_case : if LOG2N=1 generate
	begin
		rpm_had_2x2 : entity work.rpm_2x2_hadamard_grow
			generic map(W=>W, SLICE_HEIGHT=>SLICE_HEIGHT)
			port map(clk=>clk,
				left_in=>src(2*W-1 downto W), right_in=>src(W-1 downto 0),
				left_out=>dst(2*(W+1)-1 downto (W+1)), right_out=>dst((W+1)-1 downto 0)
			);
	end generate;

	inductive_case : if LOG2N>1 generate
		constant K : natural := 2**(LOG2N-1);
		signal local : std_logic_vector((W+LOG2N-1)*(2**LOG2N)-1 downto 0);
		signal local_buff : std_logic_vector((W+LOG2N-1)*(2**LOG2N)-1 downto 0);
	begin
		had_left : entity work.rpm_nxn_hadamard_grow
			generic map(W=>W, LOG2N=>LOG2N-1, SLICE_HEIGHT=>SLICE_HEIGHT)
			port map(clk=>clk,
				src=>src(W*2*K-1 downto W*K),
				dst=>local((W+LOG2N-1)*2*k-1 downto (W+LOG2N-1)*K)
			);
			
		had_right : entity work.rpm_nxn_hadamard_grow
			generic map(W=>W, LOG2N=>LOG2N-1, SLICE_HEIGHT=>SLICE_HEIGHT)
			port map(clk=>clk,
				src=>src(W*K-1 downto 0),
				dst=>local((W+LOG2N-1)*K-1 downto 0)
			);
			
		no_buff : if LOG2N <= 3 generate
			local_buff <= local;
		end generate;
		
		with_buff : if LOG2N > 3 generate
			process(clk) begin if(rising_edge(clk)) then
				local_buff <= local;
			end if; end process;
		end generate;
			
		zip : for i in 0 to 2**(LOG2N-1)-1 generate
			had : entity work.rpm_2x2_hadamard_grow
				generic map(W=>W+LOG2N-1, SLICE_HEIGHT=>SLICE_HEIGHT)
				port map(clk=>clk,
					left_in=>local_buff((i+1)*(W+LOG2N-1)-1 downto i*(W+LOG2N-1)),right_in=>local_buff((i+K+1)*(W+LOG2N-1)-1 downto (i+K)*(W+LOG2N-1)),
					left_out=>dst((2*i+1)*(W+LOG2N)-1 downto (2*i)*(W+LOG2N)), right_out=>dst((2*i+1+1)*(W+LOG2N)-1 downto (2*i+1)*(W+LOG2N))
				);
		end generate;
	end generate;
end RTL_xilinx;

