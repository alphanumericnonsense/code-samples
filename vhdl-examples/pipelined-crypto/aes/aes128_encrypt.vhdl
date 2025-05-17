----------------------------------------------
-- aes encryption with parallel key expansion
-- latency 13
----------------------------------------------
library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use work.aes_pkg.all;
use std.textio.all;

entity aes128_encrypt is
    port(
        clk_i : in std_logic := '0';
        key_i : in std_logic_vector(0 to 128-1) := (others => '0');
        plain_i : in std_logic_vector(0 to 128-1) := (others => '0');
        cipher_o : out std_logic_vector(0 to 128-1) := (others => '0');
        arstn_i : std_logic := '0'
    );
end entity aes128_encrypt;

architecture pipe of aes128_encrypt is
    signal state0, state1, state2,  state3, state4, state5, state6, state7, state8, state9, state10 : aes_state_t := (others => (others => (others => '0')));
    signal rk0, rk1, rk2, rk3, rk4, rk5, rk6, rk7, rk8, rk9, rk10 : aes_state_t := (others => (others => (others => '0')));
begin
    PIPE : process(clk_i, arstn_i) is

    begin
        if NOT arstn_i then
            cipher_o <= (others => '0');

            state0 <= (others => (others => (others => '0')));
            state1 <= (others => (others => (others => '0')));
            state2 <= (others => (others => (others => '0')));
            state3 <= (others => (others => (others => '0')));
            state4 <= (others => (others => (others => '0')));
            state5 <= (others => (others => (others => '0')));
            state6 <= (others => (others => (others => '0')));
            state7 <= (others => (others => (others => '0')));
            state8 <= (others => (others => (others => '0')));
            state9 <= (others => (others => (others => '0')));
            state10 <= (others => (others => (others => '0')));

            rk0 <= (others => (others => (others => '0')));
            rk1 <= (others => (others => (others => '0')));
            rk2 <= (others => (others => (others => '0')));
            rk3 <= (others => (others => (others => '0')));
            rk4 <= (others => (others => (others => '0')));
            rk5 <= (others => (others => (others => '0')));
            rk6 <= (others => (others => (others => '0')));
            rk7 <= (others => (others => (others => '0')));
            rk8 <= (others => (others => (others => '0')));
            rk9 <= (others => (others => (others => '0')));
            rk10 <= (others => (others => (others => '0')));
        elsif rising_edge(clk_i) then
            -- encryption
            state0 <= slv_to_state_f(plain_i);
            state1 <= aes128_round_f(state0, rk0, 0);
            state2 <= aes128_round_f(state1, rk1, 1);
            state3 <= aes128_round_f(state2, rk2, 2);
            state4 <= aes128_round_f(state3, rk3, 3);
            state5 <= aes128_round_f(state4, rk4, 4);
            state6 <= aes128_round_f(state5, rk5, 5);
            state7 <= aes128_round_f(state6, rk6, 6);
            state8 <= aes128_round_f(state7, rk7, 7);
            state9 <= aes128_round_f(state8, rk8, 8);
            state10 <= aes128_round_f(state9, rk9, 9);

            cipher_o <= state_to_slv_f(aes128_round_f(state10, rk10, 10));

            -- key expansion
            rk0 <= slv_to_state_f(key_i);
            rk1 <= aes128_key_exp_round_f(rk0, 1);
            rk2 <= aes128_key_exp_round_f(rk1, 2);
            rk3 <= aes128_key_exp_round_f(rk2, 3);
            rk4 <= aes128_key_exp_round_f(rk3, 4);
            rk5 <= aes128_key_exp_round_f(rk4, 5);
            rk6 <= aes128_key_exp_round_f(rk5, 6);
            rk7 <= aes128_key_exp_round_f(rk6, 7);
            rk8 <= aes128_key_exp_round_f(rk7, 8);
            rk9 <= aes128_key_exp_round_f(rk8, 9);
            rk10 <= aes128_key_exp_round_f(rk9, 10);

            -- report "round0 " & to_hstring(state_to_slv_f(state0)) & " " & to_hstring(state_to_slv_f(rk0));
            -- report "round1 " & to_hstring(state_to_slv_f(state1)) & " " & to_hstring(state_to_slv_f(rk1));
            -- report "round2 " & to_hstring(state_to_slv_f(state2)) & " " & to_hstring(state_to_slv_f(rk2));
            -- report "round3 " & to_hstring(state_to_slv_f(state3)) & " " & to_hstring(state_to_slv_f(rk3));
            -- report "round4 " & to_hstring(state_to_slv_f(state4)) & " " & to_hstring(state_to_slv_f(rk4));
            -- report "round5 " & to_hstring(state_to_slv_f(state5)) & " " & to_hstring(state_to_slv_f(rk5));
            -- report "round6 " & to_hstring(state_to_slv_f(state6)) & " " & to_hstring(state_to_slv_f(rk6));
            -- report "round7 " & to_hstring(state_to_slv_f(state7)) & " " & to_hstring(state_to_slv_f(rk7));
            -- report "round8 " & to_hstring(state_to_slv_f(state8)) & " " & to_hstring(state_to_slv_f(rk8));
            -- report "round9 " & to_hstring(state_to_slv_f(state9)) & " " & to_hstring(state_to_slv_f(rk9));
            -- report "round10 " & to_hstring(state_to_slv_f(state10)) & " " & to_hstring(state_to_slv_f(rk10));
            -- report to_hstring(key_i) & " " & to_hstring(plain_i) & " " & to_hstring(cipher_o) & LF;
            --
            --
            -- report "arg to cipher_o " & to_hstring(state_to_slv_f(aes128_round_f(state10, rk10, 10)));
            -- report "cipher_o " & to_hstring(cipher_o);
            -- report to_hstring(aes128_round_f(state10, rk10, 10)(0,0));
        else
            NULL;
        end if;
    end process PIPE;
end architecture pipe;
