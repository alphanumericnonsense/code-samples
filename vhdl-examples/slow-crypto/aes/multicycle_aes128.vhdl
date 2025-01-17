library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use work.aes128_enc_pkg.all;

entity multicycle_aes128 is
    port(
        clk_in : in std_logic;
        rstn_in : in std_logic;
        key_in : in std_logic_vector(0 to 128-1);
        plain_in : in std_logic_vector(0 to 128-1);
        cipher_out : out std_logic_vector(0 to 128-1);
        cipher_valid_out : out std_logic
    );
end entity multicycle_aes128;

architecture oneroundpercycle of multicycle_aes128 is

begin
    MAIN : process(clk_in, rstn_in) is
        variable state_v : aes_state_t := (others => (others => (others => '0')));
        variable key_v : aes_state_t := (others => (others => (others => '0')));
        variable cipher_out_v : std_logic_vector(0 to 128-1) := (others => '0');
        variable round_v : natural := 0;
    begin
        if rising_edge(clk_in) then
            if NOT rstn_in then
                round_v := 0;
                state_v := (others => (others => (others => '0')));
                key_v := (others => (others => (others => '0')));
                -- cipher_out <= (others => '0');
                -- cipher_valid_out <= '0';
            else
                --report "state, key: " & to_hstring(state_to_slv_f(state_v)) & ", " & to_hstring(state_to_slv_f(key_v));
                case round_v is
                    when 0 =>
                        state_v := slv_to_state_f(plain_in);
                        key_v := slv_to_state_f(key_in);
                        state_v := aes128_round_f(state_v, key_v, 0);
                        round_v := round_v + 1;
                        cipher_valid_out <= '0';
                    when 10 =>
                        key_v := aes128_key_exp_round_f(key_v, 10);
                        state_v := aes128_round_f(state_v, key_v, 10);
                        cipher_out_v := state_to_slv_f(state_v);
                        round_v := 0; -- auto-reset
                        cipher_valid_out <= '1';
                        --report "aes done";
                    when others =>
                        key_v := aes128_key_exp_round_f(key_v, round_v);
                        state_v := aes128_round_f(state_v, key_v, round_v);
                        round_v := round_v + 1;
                        cipher_valid_out <= '0';
                end case;
                cipher_out <= cipher_out_v;
            end if;
        end if;
    end process MAIN;
end architecture oneroundpercycle;
