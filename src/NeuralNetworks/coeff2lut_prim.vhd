-- start of coeff2lut_prim converter 
library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity coeff2lut_prim is
 generic (Bc : integer := 3);
  port(
    clk_i: in std_logic;
    rst_i: in std_logic;
    ce_i: in std_logic;
    coeff_i: in std_logic_vector(Bc-1 downto 0);
    ctrl_int : in std_logic_vector(2 downto 0);
    cdo_lsb_o: out std_logic_vector(Bc+4-1 downto 0);
    cdo_msb_o: out std_logic_vector(Bc+4-1 downto 0)
  );
end coeff2lut_prim;
 
architecture conversion of coeff2lut_prim is
signal coeff_internal, cdo_temp_lsb, cdo_temp_msb,cdo_int_lsb, cdo_int_msb, add_int_lsb, add_int_msb : std_logic_vector(Bc+5-1 downto 0);

component mux_add1 is
  generic(
	input_a1_word_size : integer := Bc+5;
	input_b1_word_size : integer := Bc+5;
	input_b2_word_size : integer := Bc+5;
	output_word_size : integer := Bc+5;
	input_a1_left_shift : integer := 0;
	input_b1_left_shift : integer := 0;
	input_b2_left_shift : integer := 0;
	use_output_ff : boolean := false
  );

  port (
	rst_i: in std_logic;
	clk_i: in std_logic;
	c_in_addsub: in std_logic; --add = '0', sub = '1'
	a1_i: in std_logic_vector(input_a1_word_size-1 downto 0);
	b1_i: in std_logic_vector(input_b1_word_size-1 downto 0);
	b2_i: in std_logic_vector(input_b2_word_size-1 downto 0);
	mux_sel_0: in std_logic;
	mux_sel_1: in std_logic;
	sum_o  : out std_logic_vector(output_word_size-1 downto 0)
  );

end component;

begin

coeff_internal <=  std_logic_vector(shift_left(resize(signed(coeff_i),Bc+5),1)) when (ctrl_int(2) = '1') else std_logic_vector(resize(signed(coeff_i),Bc+5));

mux_add_lsb: mux_add1
generic map(
	input_b1_left_shift => 0,
	input_b2_left_shift => 4
)
port map(
	rst_i => rst_i,
	clk_i=> clk_i,
	c_in_addsub => '1',
	a1_i => coeff_internal,
	b1_i => cdo_int_lsb,
	b2_i => coeff_internal,
	mux_sel_0 => ctrl_int(0),
	mux_sel_1 => '0',
	sum_o => cdo_temp_lsb
);

mux_add_msb: mux_add1
generic map(
	input_b1_left_shift => 0,
	input_b2_left_shift => 3
)
port map(
	rst_i => rst_i,
	clk_i=> clk_i,
	c_in_addsub => '1',
	a1_i => coeff_internal,
	b1_i => cdo_int_msb,
	b2_i => coeff_internal,
	mux_sel_0 => ctrl_int(1),
	mux_sel_1 => ctrl_int(0),
	sum_o => cdo_temp_msb
);


save_int: process(clk_i)
begin
  if clk_i'event and clk_i = '1' then
    if rst_i = '1' then
      cdo_int_lsb<= (others => '0');
      cdo_int_msb<= (others => '0');
    else
      if ce_i = '1' then
        cdo_int_lsb <= cdo_temp_lsb;
        cdo_int_msb <= cdo_temp_msb;
    end if;  
    end if;
  end if;
end process;

cdo_lsb_o <= cdo_temp_lsb(cdo_temp_lsb'length-1 downto 1);
cdo_msb_o <= cdo_temp_msb(cdo_temp_msb'length-1 downto 1);


-- define control signals

end conversion;

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

library unisim;
use unisim.vcomponents.all;

entity mux_add1 is
  generic(
	input_a1_word_size : integer := 8;
	input_b1_word_size : integer := 8;
	input_b2_word_size : integer := 8;
	output_word_size : integer := 9;
	input_a1_left_shift : integer := 0;
	input_b1_left_shift : integer := 0;
	input_b2_left_shift : integer := 0;
	use_output_ff : boolean := false
  );

  port (
	rst_i: in std_logic;
	clk_i: in std_logic;
	c_in_addsub: in std_logic; --add = '0', sub = '1'
	a1_i: in std_logic_vector(input_a1_word_size-1 downto 0);
	b1_i: in std_logic_vector(input_b1_word_size-1 downto 0);
	b2_i: in std_logic_vector(input_b2_word_size-1 downto 0);
	mux_sel_0: in std_logic;
	mux_sel_1: in std_logic;
	sum_o  : out std_logic_vector(output_word_size-1 downto 0)
  );

end mux_add1;

-- ########################################################################################

architecture mux_add1_arch of mux_add1 is

  constant no_of_slices : integer := ((output_word_size - 1) / 4) + 1;
  constant bits_used_in_msb_slice : integer := 4-4*no_of_slices+output_word_size;
  signal a1 : std_logic_vector(output_word_size-1 downto 0);
  signal b1 : std_logic_vector(output_word_size-1 downto 0);
  signal b2 : std_logic_vector(output_word_size-1 downto 0);
  signal z : std_logic_vector(output_word_size-1 downto 0); --result before FFs
  signal sum_sh : std_logic_vector(output_word_size-1 downto 0); --result after FFs, if used

  type carry_vec is array (0 to no_of_slices-2) of std_logic_vector (3 downto 0);
  signal carry : carry_vec; --carry signals between carry4 elements
  signal D : std_logic_vector(4*no_of_slices-1 downto 0); --carry-MUX data input (vector for all slices)
  signal S : std_logic_vector(4*no_of_slices-1 downto 0); --carry-MUX select line (vector for all slices)
  signal sum_msb : std_logic_vector(3 downto 0);

begin

  a1 <= std_logic_vector(shift_left(resize(signed(a1_i),output_word_size),input_a1_left_shift));
  b1 <= std_logic_vector(shift_left(resize(signed(b1_i),output_word_size),input_b1_left_shift));
  b2 <= std_logic_vector(shift_left(resize(signed(b2_i),output_word_size),input_b2_left_shift));

  lut_gen: for i in 0 to output_word_size-1 generate
      lut_bit_i: LUT6_2
        generic map(
          INIT => X"A555A59955555555" --init hex value for LUT6_2 prim
        )
        port map(
          I0  => a1(i),
          I1  => b1(i),
          I2  => b2(i),
          I3  => mux_sel_0,
          I4  => mux_sel_1,
          I5  => '1',
          O5 => D(i),
          O6 => S(i)
        );
  end generate;

  zero_gen: for i in output_word_size to 4*no_of_slices-1 generate
    D(i) <= '0';
    S(i) <= '0';
  end generate;

  --CARRY CHAIN ##############################################################

  if_more_than_one_slice: if no_of_slices > 1 generate
	--first carry4:
	carry4_inst_0 : carry4
	  port map (
		co => carry(0),             -- 4-bit carry out
		o => z(3 downto 0),         -- 4-bit carry chain xor data out
		ci => '0',                -- 1-bit carry cascade input
		cyinit => c_in_addsub,      -- 1-bit carry initialization
		di => D(3 downto 0),        -- 4-bit carry-mux data in
		s => S(3 downto 0)          -- 4-bit carry-mux select input
	);

    --middle carry4:
    carry4_gen: for i in 1 to no_of_slices-2 generate
      carry4_inst_i : carry4
        port map (
        co => carry(i),                       -- 4-bit carry out
        o => z(4*(i+1)-1 downto 4*i),         -- 4-bit carry chain xor data out
        ci => carry(i-1)(3),                  -- 1-bit carry cascade input
        cyinit => '0',                      -- 1-bit carry initialization
        di => D(4*(i+1)-1 downto 4*i),        -- 4-bit carry-mux data in
        s => S(4*(i+1)-1 downto 4*i)          -- 4-bit carry-mux select input
      );
    end generate;

    --last carry4:
    carry4_inst_N : carry4
      port map (
      co => open,                                          -- 4-bit carry out
      o => sum_msb,                                        -- 4-bit carry chain xor data out
      ci => carry(no_of_slices-2)(3),                      -- 1-bit carry cascade input
      cyinit => '0',                                     -- 1-bit carry initialization
      di => D(4*no_of_slices-1 downto 4*(no_of_slices-1)), -- 4-bit carry-mux data in
      s => S(4*no_of_slices-1 downto 4*(no_of_slices-1))   -- 4-bit carry-mux select input
    );
  end generate;

--#######################################################################################

  z(output_word_size-1 downto output_word_size-bits_used_in_msb_slice) <= sum_msb(bits_used_in_msb_slice-1 downto 0);
  sum_sh <= z(output_word_size-1 downto 0);

  --implement flipflops when use_output_ff
  use_ff: if use_output_ff = true generate
    ff_gen: for i in 0 to output_word_size-1 generate
      ff_i: fdce
        generic map(
          -- initialize all flip flops with '0'
          init => '0'
        )
        port map(
          clr => rst_i,
          ce  => '1',
          d   => sum_sh(i),
          c   => clk_i,
          q   => sum_o(i)
        );
    end generate;
  end generate;
  -- otherwise bypass the flipflops
  bypass_ff: if use_output_ff = false generate
    sum_o <= sum_sh;
  end generate;
end architecture;
