-- start of coeff2lut_no_prim converter 
library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity coeff2lut_no_prim is
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
end coeff2lut_no_prim;

architecture conversion of coeff2lut_no_prim is
signal coeff_internal, cdo_temp_lsb, cdo_temp_msb,cdo_int_lsb, cdo_int_msb, add_int_lsb, add_int_msb : std_logic_vector(Bc+5-1 downto 0);

begin

coeff_internal <=  std_logic_vector(shift_left(resize(signed(coeff_i),Bc+5),1)) when (ctrl_int(2) = '1') else std_logic_vector(resize(signed(coeff_i),Bc+5));

add_int_lsb <= std_logic_vector(shift_left(signed(coeff_internal),4)) when (ctrl_int(0) = '1') else cdo_int_lsb;
add_int_msb <= std_logic_vector(shift_left(signed(coeff_internal),3)) when (ctrl_int(1 downto 0) = "10") else
            cdo_int_msb when (ctrl_int(1 downto 0) = "00") else
            (others => '0');

cdo_temp_lsb <= std_logic_vector(signed(add_int_lsb) - signed(coeff_internal));
cdo_temp_msb <= std_logic_vector(signed(add_int_msb) - signed(coeff_internal));

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
