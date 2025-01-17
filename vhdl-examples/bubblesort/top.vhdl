library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use work.array_types.all;

entity top_e is
   generic(
       depth_g : positive := 4;
       width_g : positive := 4
   );
   port(
       clk  : in std_logic := '0';
       led : out std_logic_vector(15 downto 0) := (others => '0');
       sw : in std_logic_vector(15 downto 0) := (others => '0')
   );
end entity top_e;

architecture test of top_e is
   signal start_i_s, done_o_s : std_logic := '0';
   signal array_i_s, array_o_s : slv_1d_t(0 to depth_g - 1)(width_g - 1 downto 0) := (others => (others => '0'));
begin
   DUT : entity work.bubblesort_e
       generic map(
           depth_g => depth_g,
           width_g => width_g
       )
       port map(
           clk_i => clk,
           start_i => start_i_s,
           array_i => array_i_s,
           array_o => array_o_s,
           done_o => done_o_s
       );

   -- input from switches, output sorted results to LEDs
   led <= array_o_s(0) & array_o_s(1) & array_o_s(2) & array_o_s(3);
   array_i_s <= (sw(15 downto 12), sw(11 downto 8), sw(7 downto 4), sw(3 downto 0));

   CONTROL : process(clk) is
       variable sw_prev : std_logic_vector(15 downto 0) := (others => '0');
   begin
       if rising_edge(clk) then
           if sw_prev /= sw then
               start_i_s <= '1';
           else
               start_i_s <= '0';
           end if;
           sw_prev := sw;
       end if;
   end process CONTROL;
end architecture test;
