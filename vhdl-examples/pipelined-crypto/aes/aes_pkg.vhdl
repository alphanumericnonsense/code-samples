library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

package aes_pkg is
    -- sbox LUT
    type sbox_lut_t is array(0 to 256-1) of std_logic_vector(8-1 downto 0);
    constant sbox_c : sbox_lut_t := (
        X"63", X"7c", X"77", X"7b", X"f2", X"6b", X"6f", X"c5", X"30", X"01", X"67", X"2b", X"fe", X"d7", X"ab", X"76",
        X"ca", X"82", X"c9", X"7d", X"fa", X"59", X"47", X"f0", X"ad", X"d4", X"a2", X"af", X"9c", X"a4", X"72", X"c0",
        X"b7", X"fd", X"93", X"26", X"36", X"3f", X"f7", X"cc", X"34", X"a5", X"e5", X"f1", X"71", X"d8", X"31", X"15",
        X"04", X"c7", X"23", X"c3", X"18", X"96", X"05", X"9a", X"07", X"12", X"80", X"e2", X"eb", X"27", X"b2", X"75",
        X"09", X"83", X"2c", X"1a", X"1b", X"6e", X"5a", X"a0", X"52", X"3b", X"d6", X"b3", X"29", X"e3", X"2f", X"84",
        X"53", X"d1", X"00", X"ed", X"20", X"fc", X"b1", X"5b", X"6a", X"cb", X"be", X"39", X"4a", X"4c", X"58", X"cf",
        X"d0", X"ef", X"aa", X"fb", X"43", X"4d", X"33", X"85", X"45", X"f9", X"02", X"7f", X"50", X"3c", X"9f", X"a8",
        X"51", X"a3", X"40", X"8f", X"92", X"9d", X"38", X"f5", X"bc", X"b6", X"da", X"21", X"10", X"ff", X"f3", X"d2",
        X"cd", X"0c", X"13", X"ec", X"5f", X"97", X"44", X"17", X"c4", X"a7", X"7e", X"3d", X"64", X"5d", X"19", X"73",
        X"60", X"81", X"4f", X"dc", X"22", X"2a", X"90", X"88", X"46", X"ee", X"b8", X"14", X"de", X"5e", X"0b", X"db",
        X"e0", X"32", X"3a", X"0a", X"49", X"06", X"24", X"5c", X"c2", X"d3", X"ac", X"62", X"91", X"95", X"e4", X"79",
        X"e7", X"c8", X"37", X"6d", X"8d", X"d5", X"4e", X"a9", X"6c", X"56", X"f4", X"ea", X"65", X"7a", X"ae", X"08",
        X"ba", X"78", X"25", X"2e", X"1c", X"a6", X"b4", X"c6", X"e8", X"dd", X"74", X"1f", X"4b", X"bd", X"8b", X"8a",
        X"70", X"3e", X"b5", X"66", X"48", X"03", X"f6", X"0e", X"61", X"35", X"57", X"b9", X"86", X"c1", X"1d", X"9e",
        X"e1", X"f8", X"98", X"11", X"69", X"d9", X"8e", X"94", X"9b", X"1e", X"87", X"e9", X"ce", X"55", X"28", X"df",
        X"8c", X"a1", X"89", X"0d", X"bf", X"e6", X"42", X"68", X"41", X"99", X"2d", X"0f", X"b0", X"54", X"bb", X"16");

    -- 4x4 matrix of bytes
    -- used for state and round keys
    type aes_state_t is array (0 to 4-1, 0 to 4-1) of std_logic_vector(8-1 downto 0);

    -- functions
    function rijmult02_f(a : std_logic_vector(8-1 downto 0)) return std_logic_vector;
    function rijmult03_f(a : std_logic_vector(8-1 downto 0)) return std_logic_vector;
    function add_round_key_f(state_i, round_key : aes_state_t) return aes_state_t;
    function sub_bytes_f(state_i : aes_state_t) return aes_state_t;
    function shift_rows_f(state_i : aes_state_t) return aes_state_t;
    function mix_cols_f(state_i : aes_state_t) return aes_state_t;
    function aes128_round_f(state_i, round_key : aes_state_t; round : integer) return aes_state_t;
    function aes128_key_exp_round_f(round_key_prev : aes_state_t; round : integer) return aes_state_t;
    function slv_to_state_f(slv_i : std_logic_vector(0 to 128-1)) return aes_state_t;
    function state_to_slv_f(state_i : aes_state_t) return std_logic_vector;
    -- aes128 encrypt test function
    function aes128_encrypt_f(plain_i, key_i : std_logic_vector(0 to 128-1)) return std_logic_vector;
end package aes_pkg;

package body aes_pkg is
    function rijmult02_f(a : std_logic_vector(8-1 downto 0)) return std_logic_vector is
        variable b : std_logic_vector(8-1 downto 0) := (others => '0');
    begin
        ------------------
        -- multiply by x02
        ------------------
        b(0) := a(7);
        b(1) := a(0) XOR a(7);
        b(2) := a(1);
        b(3) := a(2) XOR a(7);
        b(4) := a(3) XOR a(7);
        b(5) := a(4);
        b(6) := a(5);
        b(7) := a(6);
        return b;
    end function rijmult02_f;

    function rijmult03_f(a : std_logic_vector(8-1 downto 0)) return std_logic_vector is
        variable b : std_logic_vector(8-1 downto 0) := (others => '0');
    begin
        ------------------
        -- multiply by x03
        ------------------
        b(0) := a(0) XOR a(7);
        b(1) := a(0) XOR a(1) XOR a(7);
        b(2) := a(1) XOR a(2);
        b(3) := a(2) XOR a(3) XOR a(7);
        b(4) := a(3) XOR a(4) XOR a(7);
        b(5) := a(4) XOR a(5);
        b(6) := a(5) XOR a(6);
        b(7) := a(6) XOR a(7);
        return b;
    end function rijmult03_f;

    function add_round_key_f(state_i, round_key : aes_state_t) return aes_state_t is
        variable state_o : aes_state_t := (others => (others => (others => '0')));
    begin
        for i in 0 to 3 loop
            for j in 0 to 3 loop
                state_o(i,j) := state_i(i,j) XOR round_key(i,j);
                --report to_string(round_key(i,j)) & " " & to_string(state_i(i,j)) & " " & to_string(state_o(i,j));
            end loop;
        end loop;
        return state_o;
    end function add_round_key_f;

    function sub_bytes_f(state_i : aes_state_t) return aes_state_t is
        variable state_o : aes_state_t := (others => (others => (others => '0')));
    begin
        for i in 0 to 3 loop
            for j in 0 to 3 loop
                state_o(i,j) := sbox_c(to_integer(unsigned(state_i(i,j))));
            end loop;
        end loop;
        return state_o;
    end function sub_bytes_f;

    function shift_rows_f(state_i : aes_state_t) return aes_state_t is
        variable state_o : aes_state_t := (others => (others => (others => '0')));
    begin
        for i in 0 to 3 loop
            for j in 0 to 3 loop
                state_o(i,j) := state_i(i, (j+i) mod 4);
            end loop;
        end loop;
        return state_o;
    end function shift_rows_f;

    function mix_cols_f(state_i : aes_state_t) return aes_state_t is
        variable state_o : aes_state_t := (others => (others => (others => '0')));
    begin
        for j in 0 to 3 loop
            state_o(0,j) := rijmult02_f(state_i(0,j)) XOR rijmult03_f(state_i(1,j)) XOR state_i(2,j) XOR state_i(3,j);
            state_o(1,j) := state_i(0,j) XOR rijmult02_f(state_i(1,j)) XOR rijmult03_f(state_i(2,j)) XOR state_i(3,j);
            state_o(2,j) := state_i(0,j) XOR state_i(1,j) XOR rijmult02_f(state_i(2,j)) XOR rijmult03_f(state_i(3,j));
            state_o(3,j) := rijmult03_f(state_i(0,j)) XOR state_i(1,j) XOR state_i(2,j) XOR rijmult02_f(state_i(3,j));
        end loop;
        return state_o;
    end function mix_cols_f;

    function aes128_round_f(state_i, round_key : aes_state_t; round : integer) return aes_state_t is
        variable new_state : aes_state_t := (others => (others => (others => '0')));
    begin
        if round = 9 then
            new_state := shift_rows_f(sub_bytes_f(add_round_key_f(state_i, round_key)));
        elsif round = 10 then
            new_state := add_round_key_f(state_i, round_key);
        else
            new_state := mix_cols_f(shift_rows_f(sub_bytes_f(add_round_key_f(state_i, round_key))));
        end if;

        return new_state;
    end function aes128_round_f;

    function aes128_key_exp_round_f(round_key_prev : aes_state_t; round : integer) return aes_state_t is
        variable round_key_next : aes_state_t := (others => (others => (others => '0')));
        -- could trim down as mostly zero and easy-to-generate powers of x
        type rc_t is array(1 to 10, 0 to 3) of std_logic_vector(8-1 downto 0);
        constant rc_c : rc_t := ((X"01", X"00", X"00", X"00"),
                                (X"02", X"00", X"00", X"00"),
                                (X"04", X"00", X"00", X"00"),
                                (X"08", X"00", X"00", X"00"),
                                (X"10", X"00", X"00", X"00"),
                                (X"20", X"00", X"00", X"00"),
                                (X"40", X"00", X"00", X"00"),
                                (X"80", X"00", X"00", X"00"),
                                (X"1B", X"00", X"00", X"00"),
                                (X"36", X"00", X"00", X"00"));
    begin
        for j in 0 to 3 loop
            for i in 0 to 3 loop
                if j = 0 then
                    round_key_next(i,0) := round_key_prev(i,0) XOR rc_c(round, i) XOR sbox_c(to_integer(unsigned(round_key_prev((i+1) mod 4, 3))));
                else
                    round_key_next(i,j) := round_key_prev(i,j) XOR round_key_next(i, j - 1);
                end if;
            end loop;
        end loop;
        return round_key_next;
    end function aes128_key_exp_round_f;

    function slv_to_state_f(slv_i : std_logic_vector(0 to 128-1)) return aes_state_t is
        variable state_o : aes_state_t := (others => (others => (others => '0')));
    begin
        -- column major
        for i in 0 to 3 loop
            for j in 0 to 3 loop
                state_o(i,j) := slv_i(8*(4*j+i) to 8*(4*j+i+1)-1);
            end loop;
        end loop;
        return state_o;
    end function slv_to_state_f;

    function state_to_slv_f(state_i : aes_state_t) return std_logic_vector is
        variable slv_o : std_logic_vector(0 to 128-1) := (others => '0');
    begin
        -- column major
        for i in 0 to 3 loop
            for j in 0 to 3 loop
                slv_o(8*(4*j+i) to 8*(4*j+i+1)-1) := state_i(i,j);
            end loop;
        end loop;
        return slv_o;
    end function state_to_slv_f;

    function aes128_encrypt_f(plain_i, key_i : std_logic_vector(0 to 128-1)) return std_logic_vector is
        variable state_v : aes_state_t := (others => (others => (others => '0')));
        variable key_v : aes_state_t := (others => (others => (others => '0')));
    begin
        state_v := slv_to_state_f(plain_i);
        key_v := slv_to_state_f(key_i);
        state_v := aes128_round_f(state_v, key_v, 0);
        for i in 1 to 10 loop
            key_v := aes128_key_exp_round_f(key_v, i);
            state_v := aes128_round_f(state_v, key_v, i);
        end loop;
        return state_to_slv_f(state_v);
    end function aes128_encrypt_f;

end package body aes_pkg;
