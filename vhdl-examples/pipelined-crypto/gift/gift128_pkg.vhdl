-- https://eprint.iacr.org/2017/622.pdf
-- https://giftcipher.github.io/gift/

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

package gift128_pkg is
    type sbox_t is array(0 to 16-1) of std_logic_vector(4-1 downto 0);
    constant sbox_c : sbox_t := (x"1", x"a", x"4", x"c", x"6", x"f", x"3", x"9",
                                 x"2", x"d", x"b", x"7", x"5", x"0", x"8", x"e"); -- correct
    type perm_t is array(0 to 128-1) of natural;
    constant perm_c : perm_t := (0, 33, 66, 99, 96, 1, 34, 67, 64, 97, 2, 35, 32, 65, 98, 3,
                                4, 37, 70, 103, 100, 5, 38, 71, 68, 101, 6, 39, 36, 69, 102, 7,
                                8, 41, 74, 107, 104, 9, 42, 75, 72, 105, 10, 43, 40, 73, 106, 11,
                                12, 45, 78, 111, 108, 13, 46, 79, 76, 109, 14, 47, 44, 77, 110, 15,
                                16, 49, 82, 115, 112, 17, 50, 83, 80, 113, 18, 51, 48, 81, 114, 19,
                                20, 53, 86, 119, 116, 21, 54, 87, 84, 117, 22, 55, 52, 85, 118, 23,
                                24, 57, 90, 123, 120, 25, 58, 91, 88, 121, 26, 59, 56, 89, 122, 27,
                                28, 61, 94, 127, 124, 29, 62, 95, 92, 125, 30, 63, 60, 93, 126, 31); -- correct
    type rc_t is array(0 to 39) of std_logic_vector(5 downto 0);
    constant rc_c : rc_t :=     ("000001",
                                "000011",
                                "000111",
                                "001111",
                                "011111",
                                "111110",
                                "111101",
                                "111011",
                                "110111",
                                "101111",
                                "011110",
                                "111100",
                                "111001",
                                "110011",
                                "100111",
                                "001110",
                                "011101",
                                "111010",
                                "110101",
                                "101011",
                                "010110",
                                "101100",
                                "011000",
                                "110000",
                                "100001",
                                "000010",
                                "000101",
                                "001011",
                                "010111",
                                "101110",
                                "011100",
                                "111000",
                                "110001",
                                "100011",
                                "000110",
                                "001101",
                                "011011",
                                "110110",
                                "101101",
                                "011010"); -- correct
    function key_round_f(key_prev : std_logic_vector(128-1 downto 0)) return std_logic_vector;
    function add_round_key_f(U, V : std_logic_vector(32-1 downto 0); state_i : std_logic_vector(128-1 downto 0); round : integer)
        return std_logic_vector;
    function perm_f(state_i : std_logic_vector(128-1 downto 0)) return std_logic_vector;
    function sbox_f(state_i : std_logic_vector(128-1 downto 0)) return std_logic_vector;
    function gift128_round_f(state_i : std_logic_vector(128-1 downto 0); U, V : std_logic_vector(32-1 downto 0) ; round : integer)
        return std_logic_vector;
    -- experimental, trying to do two rounds per cycle
    function double_key_round_f(key_prev : std_logic_vector(128-1 downto 0)) return std_logic_vector;
    function double_gift128_round_f(state_i : std_logic_vector(128-1 downto 0); U, V : std_logic_vector(64-1 downto 0) ; round : integer)
        return std_logic_vector;
end package gift128_pkg;

package body gift128_pkg is
    function key_round_f(key_prev : std_logic_vector(128-1 downto 0)) return std_logic_vector is
        variable key_next : std_logic_vector(128-1 downto 0) := (others => '0');
    begin
        for i in 0 to 6-1 loop
            key_next(16*(i+1)-1 downto 16*i) := key_prev(16*(i+3)-1 downto 16*(i+2));
        end loop;
        key_next(128-1 downto 128-16) := key_prev(32-1 downto 16) ror 2;
        key_next(128-16-1 downto 128-32) := key_prev(16-1 downto 0) ror 12;
        return key_next;
    end function key_round_f;

    -- experimental, two round keys per cycle
    function double_key_round_f(key_prev : std_logic_vector(128-1 downto 0)) return std_logic_vector is
        variable x, y : std_logic_vector(128-1 downto 0) := (others => '0');
    begin
         x := key_round_f(key_prev);
         y := key_round_f(x);
         return x & y;
    end function double_key_round_f;

    function add_round_key_f(U, V : std_logic_vector(32-1 downto 0); state_i : std_logic_vector(128-1 downto 0); round : integer)
    return std_logic_vector is
        variable state_o : std_logic_vector(128-1 downto 0) := (others => '0');
    begin
        state_o := state_i;
        for i in 0 to 32-1 loop
            state_o(4*i+1) := state_o(4*i+1) XOR V(i);
            state_o(4*i+2) := state_o(4*i+2) XOR U(i);
        end loop;
        state_o(23) := state_o(23) XOR rc_c(round)(5);
        state_o(19) := state_o(19) XOR rc_c(round)(4);
        state_o(15) := state_o(15) XOR rc_c(round)(3);
        state_o(11) := state_o(11) XOR rc_c(round)(2);
        state_o(7) := state_o(7) XOR rc_c(round)(1);
        state_o(3) := state_o(3) XOR rc_c(round)(0);
        state_o(127) := state_o(127) XOR '1';
        return state_o;
    end function add_round_key_f;

    function perm_f(state_i : std_logic_vector(128-1 downto 0)) return std_logic_vector is
        variable state_o : std_logic_vector(128-1 downto 0) := (others => '0');
    begin
        for i in 0 to 128-1 loop
            state_o(perm_c(i)) := state_i(i);
        end loop;
        return state_o;
    end function perm_f;

    function sbox_f(state_i : std_logic_vector(128-1 downto 0)) return std_logic_vector is
        variable state_o : std_logic_vector(128-1 downto 0) := (others => '0');
    begin
        for i in 0 to 32-1 loop
            state_o(4*(i+1)-1 downto 4*i) := sbox_c(to_integer(unsigned(state_i(4*(i+1)-1 downto 4*i))));
        end loop;
        return state_o;
    end function sbox_f;

    function gift128_round_f(state_i : std_logic_vector(128-1 downto 0); U, V : std_logic_vector(32-1 downto 0) ; round : integer)
    return std_logic_vector is

    begin
        return add_round_key_f(U, V, perm_f(sbox_f(state_i)), round);
    end function gift128_round_f;

    -- experimental, two rounds per cycle
    function double_gift128_round_f(state_i : std_logic_vector(128-1 downto 0); U, V : std_logic_vector(64-1 downto 0) ; round : integer)
    return std_logic_vector is
        variable x, y : std_logic_vector(128-1 downto 0) := (others => '0');
    begin
        x := gift128_round_f(state_i, U(64-1 downto 32), V(64-1 downto 32), round);
        y := gift128_round_f(x, U(32-1 downto 0), V(32-1 downto 0), round + 1);
        return y;
    end function double_gift128_round_f;
end package body gift128_pkg;
