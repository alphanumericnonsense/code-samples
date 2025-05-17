library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use work.chacha20_pkg.all;

entity chacha20_block is
    port(
        clk_i : in std_logic;
        arstn_i : in std_logic;
        key_i : in std_logic_vector(0 to 256-1);
        ctr_i : in std_logic_vector(0 to 64-1);
        nonce_i : in std_logic_vector(0 to 64-1);
        state_o : out std_logic_vector(0 to 512-1)
    );
end entity chacha20_block;

architecture slow of chacha20_block is
    signal state0, state1, state2, state3, state4, state5, state6, state7, state8, state9, state10, state11, state12, state13, state14, state15, state16, state17, state18, state19, state20 : chacha20_state_t := (others => (others => (others => '0')));
    constant buflen_c : integer := 20;
    type chacha20_buffer_t is array(1 to buflen_c) of chacha20_state_t;
    signal chacha20_buf : chacha20_buffer_t := (others => (others => (others => (others => '0'))));
begin
    process(clk_i, arstn_i) is
        variable init_state_v : chacha20_state_t := (others => (others => (others => '0')));
    begin
        if NOT arstn_i then
            init_state_v := (others => (others => (others => '0')));
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
            state11 <= (others => (others => (others => '0')));
            state12 <= (others => (others => (others => '0')));
            state13 <= (others => (others => (others => '0')));
            state14 <= (others => (others => (others => '0')));
            state15 <= (others => (others => (others => '0')));
            state16 <= (others => (others => (others => '0')));
            state17 <= (others => (others => (others => '0')));
            state18 <= (others => (others => (others => '0')));
            state19 <= (others => (others => (others => '0')));
            state20 <= (others => (others => (others => '0')));
            chacha20_buf <= (others => (others => (others => (others => '0'))));
            state_o <= (others => '0');
        elsif rising_edge(clk_i) then
            init_state_v := ((r1_c(0), r1_c(1), r1_c(2), r1_c(3)),
                       (key_i(0 to 32-1), key_i(32 to 64-1), key_i(64 to 96-1), key_i(96 to 128-1)),
                       (key_i(128 to 160-1), key_i(160 to 192-1), key_i(192 to 224-1), key_i(224 to 256-1)),
                       (ctr_i(0 to 32-1), ctr_i(32 to 64-1), nonce_i(0 to 32-1), nonce_i(32 to 64-1))
                       );
            for i in 1 to 3 loop -- constants aren't little endian
                for j in 0 to 3 loop
                    init_state_v(i,j) := biglittleswap_f(init_state_v(i,j));
                end loop;
            end loop;
            state0 <= init_state_v;
            chacha20_buf(1) <= state0;
            for i in 1 to buflen_c-1 loop
                chacha20_buf(i+1) <= chacha20_buf(i);
            end loop;
            state1 <= col_round_f(state0);
            state2 <= diag_round_f(state1);
            state3 <= col_round_f(state2);
            state4 <= diag_round_f(state3);
            state5 <= col_round_f(state4);
            state6 <= diag_round_f(state5);
            state7 <= col_round_f(state6);
            state8 <= diag_round_f(state7);
            state9 <= col_round_f(state8);
            state10 <= diag_round_f(state9);
            state11 <= col_round_f(state10);
            state12 <= diag_round_f(state11);
            state13 <= col_round_f(state12);
            state14 <= diag_round_f(state13);
            state15 <= col_round_f(state14);
            state16 <= diag_round_f(state15);
            state17 <= col_round_f(state16);
            state18 <= diag_round_f(state17);
            state19 <= col_round_f(state18);
            state20 <= diag_round_f(state19);
            state_o <= serialize_state_f(add_states_f(state20, chacha20_buf(20)));-- mod 2^32 add initial state chacha20_buf(20)?
        else
            NULL;
        end if;
    end process;
end architecture slow;

architecture fast of chacha20_block is
    signal state0, state2, state4, state6, state8, state10, state12, state14, state16, state18, state20 : chacha20_state_t := (others => (others => (others => '0')));
    constant buflen_c : integer := 10;
    type chacha20_buffer_t is array(1 to buflen_c) of chacha20_state_t;
    signal chacha20_buf : chacha20_buffer_t:= (others => (others => (others => (others => '0'))));
begin
    process(clk_i, arstn_i) is
        variable init_state_v : chacha20_state_t := (others => (others => (others => '0')));
    begin
        if NOT arstn_i then
            init_state_v := (others => (others => (others => '0')));
            state0 <= (others => (others => (others => '0')));
            state2 <= (others => (others => (others => '0')));
            state4 <= (others => (others => (others => '0')));
            state6 <= (others => (others => (others => '0')));
            state8 <= (others => (others => (others => '0')));
            state10 <= (others => (others => (others => '0')));
            state12 <= (others => (others => (others => '0')));
            state14 <= (others => (others => (others => '0')));
            state16 <= (others => (others => (others => '0')));
            state18 <= (others => (others => (others => '0')));
            state20 <= (others => (others => (others => '0')));
            chacha20_buf <= (others => (others => (others => (others => '0'))));
            state_o <= (others => '0');
        elsif rising_edge(clk_i) then
            init_state_v := ((r1_c(0), r1_c(1), r1_c(2), r1_c(3)),
                       (key_i(0 to 32-1), key_i(32 to 64-1), key_i(64 to 96-1), key_i(96 to 128-1)),
                       (key_i(128 to 160-1), key_i(160 to 192-1), key_i(192 to 224-1), key_i(224 to 256-1)),
                       (ctr_i(0 to 32-1), ctr_i(32 to 64-1), nonce_i(0 to 32-1), nonce_i(32 to 64-1))
                       );
            for i in 1 to 3 loop -- constants aren't little endian
                for j in 0 to 3 loop
                    init_state_v(i,j) := biglittleswap_f(init_state_v(i,j));
                end loop;
            end loop;
            state0 <= init_state_v;
            chacha20_buf(1) <= state0;
            for i in 1 to buflen_c-1 loop
                chacha20_buf(i+1) <= chacha20_buf(i);
            end loop;
            state2 <= diag_round_f(col_round_f(state0));
            state4 <= diag_round_f(col_round_f(state2));
            state6 <= diag_round_f(col_round_f(state4));
            state8 <= diag_round_f(col_round_f(state6));
            state10 <= diag_round_f(col_round_f(state8));
            state12 <= diag_round_f(col_round_f(state10));
            state14 <= diag_round_f(col_round_f(state12));
            state16 <= diag_round_f(col_round_f(state14));
            state18 <= diag_round_f(col_round_f(state16));
            state20 <= diag_round_f(col_round_f(state18));
            state_o <= serialize_state_f(add_states_f(state20, chacha20_buf(10)));-- mod 2^32 add initial state chacha20_buf(10)?
            --report "state_o " & to_hstring(state_o);
            --report "arg to state_o " & to_hstring(serialize_state_f(add_states_f(state20, chacha20_buf(10))));
        else
            NULL;
        end if;
    end process;
end architecture fast;
