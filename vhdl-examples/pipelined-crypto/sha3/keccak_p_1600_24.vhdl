library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use work.sha3_pkg.all;

entity keccak_p_1600_24 is
    port(
        clk_i : in std_logic;
        arstn_i : in std_logic;
        state_i : in std_logic_vector(0 to 1600-1);
        state_o : out std_logic_vector(0 to 1600-1)
    );
end entity keccak_p_1600_24;

architecture func_pipe of keccak_p_1600_24 is
    type pipe_state_t is array(0 to 24) of std_logic_vector(0 to 1600-1);
    signal state_s : pipe_state_t := (others => (others => '0'));
begin
    process(clk_i, arstn_i) is

    begin
        if NOT arstn_i then
            state_s <= (others => (others => '0'));
            state_o <= (others => '0');
        elsif rising_edge(clk_i) then
            state_s(0) <= state_i;
            for i in 1 to 24 loop
                state_s(i) <= keccak_round_f(state_s(i-1), i-1);
            end loop;
            state_o <= state_s(24);
        else
            NULL;
        end if;
    end process;
end architecture func_pipe;
