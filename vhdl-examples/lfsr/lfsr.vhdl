library ieee;
use ieee.std_logic_1164.all;

-- http://poincare.matf.bg.ac.rs/~ezivkovm/publications/primpol1.pdf

entity lfsr is
    generic(
        width_g : natural := 16;
        poly_g : std_logic_vector(width_g -1 downto 0) :=  "0000010010001001";
        steps_g : natural := 8;
        trunc_g : natural := 8
    );
    port(
        clk_in: in std_logic;
        rstn_in : in std_logic;
        reseed_in : in std_logic;
        seed_in : in std_logic_vector(width_g-1 downto 0);
        trunc_out : out std_logic_vector(trunc_g - 1 downto 0)
    );
end entity lfsr;

architecture fibonacci of lfsr is

begin
    fib_main : process(clk_in) is
        variable state : std_logic_vector(width_g - 1 downto 0) := (others => '0');
        variable new_bit : std_logic := '0';
    begin
        if rising_edge(clk_in) then
            if NOT rstn_in then
                state := (others => '0');
                new_bit := '0';
            else
                if reseed_in then
                    state := seed_in; -- seed_in XOR state;
                    --report "new seed: " & to_string(seed_in);
                else
                    for i in 1 to steps_g loop
                        new_bit := xor (poly_g and state);
                        state := state(width_g - 2 downto 0) & new_bit;
                    end loop;
                end if;
            end if;
            trunc_out <= state(trunc_g-1 downto 0);
        end if;
    end process fib_main;
end architecture fibonacci;
