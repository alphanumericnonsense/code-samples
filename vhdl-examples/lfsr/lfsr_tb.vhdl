library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity lfsr_tb is
    generic(
        width_g : natural := 16;
        poly_g : std_logic_vector(width_g -1 downto 0) :=  "0000010010001001";
        steps_g : natural := 1;
        trunc_g : natural := 16;
        seed_in_g : std_logic_vector(width_g - 1 downto 0) := (0 => '1', others => '0')
    );
end entity lfsr_tb;

architecture test of lfsr_tb is
    signal clk, rstn, reseed : std_logic := '0';
    signal seed_in : std_logic_vector(width_g - 1 downto 0) := seed_in_g;
    signal trunc_out : std_logic_vector(trunc_g - 1 downto 0) := (others => '0');
begin

    TEST : process is
        variable counter : natural := 0;
    begin
        rstn <= '0';
        wait until rising_edge(clk);
        rstn <= '1';
        reseed <= '1';
        wait until rising_edge(clk);
        reseed <= '0';
        report "Time to get schwifty, I mean shifty";
        while counter < 100 loop
            report to_string(trunc_out);
            counter := counter + 1;
            wait until rising_edge(clk);
        end loop;
        std.env.stop;
    end process TEST;

    DUT : entity work.lfsr
        generic map(
            width_g => width_g,
            poly_g => poly_g,
            steps_g => steps_g,
            trunc_g => trunc_g
        )
        port map(
            clk_in => clk,
            rstn_in => rstn,
            reseed_in => reseed,
            seed_in => seed_in_g,
            trunc_out => trunc_out
        );

    CLOCK : process is

    begin
        clk <= NOT clk;
        wait for 1 ns;
    end process CLOCK;

end architecture test;
