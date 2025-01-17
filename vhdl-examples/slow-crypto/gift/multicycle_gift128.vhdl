library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use work.gift128_pkg.all;

entity multicycle_gift128 is
    port(
        clk_in : in std_logic;
        rstn_in : in std_logic;
        key_in : in std_logic_vector(0 to 128-1);
        plain_in : in std_logic_vector(0 to 128-1);
        cipher_out : out std_logic_vector(0 to 128-1);
        cipher_valid_out : out std_logic
    );
end entity multicycle_gift128;

architecture tworoundspercycle of multicycle_gift128 is

begin
    MAIN : process(clk_in, rstn_in) is
        variable state_v : std_logic_vector(128-1 downto 0) := (others => '0');
        variable key_v : std_logic_vector(256-1 downto 0) := (others => '0');
        variable cipher_out_v : std_logic_vector(128-1 downto 0) := (others => '0');
        variable round_v : natural := 0;
    begin
        if rising_edge(clk_in) then
            if NOT rstn_in then
                round_v := 0;
                state_v := (others => '0');
                key_v := (others => '0');
                -- cipher_out <= (others => '0');
                -- cipher_valid_out <= '0';
            else
                --report to_hstring(state_v) & " " & to_hstring(key_v);
                case round_v is
                    when 0 =>
                        state_v := plain_in;
                        key_v := key_in & key_round_f(key_in);
                        round_v := round_v + 1;
                        cipher_valid_out <= '0';
                    when 20 =>
                        state_v := double_gift128_round_f(state_v, key_v(128+6*16-1 downto 128+4*16) & key_v(6*16-1 downto 4*16), key_v(128+32-1 downto 128+0) & key_v(32-1 downto 0), 2*round_v-2);
                        key_v := double_key_round_f(key_v(128-1 downto 0));
                        cipher_out_v := state_v;
                        round_v := 0; -- auto-reset
                        cipher_valid_out <= '1';
                        --report "gift done";
                    when others =>
                        state_v := double_gift128_round_f(state_v, key_v(128+6*16-1 downto 128+4*16) & key_v(6*16-1 downto 4*16), key_v(128+32-1 downto 128+0) & key_v(32-1 downto 0), 2*round_v-2);
                        key_v := double_key_round_f(key_v(128-1 downto 0));
                        round_v := round_v + 1;
                        cipher_valid_out <= '0';
                end case;
                cipher_out <= cipher_out_v;
            end if;
        end if;
    end process MAIN;
end architecture tworoundspercycle;
