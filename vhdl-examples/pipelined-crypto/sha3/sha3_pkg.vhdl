library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
-- some test vectors with intermediate values
-- https://csrc.nist.gov/projects/cryptographic-standards-and-guidelines/example-values
package sha3_pkg is
    function theta_f(state_i : std_logic_vector(0 to 1600 - 1)) return std_logic_vector;
    function rho_f(state_i : std_logic_vector(0 to 1600 - 1)) return std_logic_vector;
    function pi_f(state_i : std_logic_vector(0 to 1600 - 1)) return std_logic_vector;
    function chi_f(state_i : std_logic_vector(0 to 1600 - 1)) return std_logic_vector;
    function iota_f(state_i : std_logic_vector(0 to 1600 - 1); rnd : integer) return std_logic_vector;
    function keccak_round_f(state_i : std_logic_vector(0 to 1600 - 1); rnd : integer) return std_logic_vector;
    function byte_reverse_f(a : std_logic_vector) return std_logic_vector;
    function pad_f(m : std_logic_vector) return std_logic_vector;
end package sha3_pkg;

package body sha3_pkg is
    -- pads a small message (less than 1088-4 bits) with 011(0*)1
    function pad_f(m : std_logic_vector) return std_logic_vector is
        variable p : std_logic_vector(0 to 1600-1) := (others => '0');
    begin
        p(0 to m'length - 1) := m;
        p(m'length to m'length + 2 -1) := "01";
        p(m'length + 2) := '1';
        p(1088-1) := '1';
        return p;
    end function pad_f;

    function byte_reverse_f(a : std_logic_vector) return std_logic_vector is
        variable ar8 : std_logic_vector(0 to a'length-1) := (others => '0');
    begin
        for i in 0 to a'length/8 - 1 loop
            for j in 0 to 7 loop
                ar8(8*i + j) := a(8*i + 7 - j);
            end loop;
        end loop;
        return ar8;
    end function byte_reverse_f;

    function theta_f(state_i : std_logic_vector(0 to 1600 - 1)) return std_logic_vector is
        variable state_o : std_logic_vector(0 to 1600 - 1) := (others => '0');
        type plane_t is array(0 to 4, 0 to 63) of std_logic;
        variable C, D : plane_t;
    begin
        for x in 0 to 4 loop
            for z in 0 to 63 loop
                C(x,z) := state_i(64*(5*0+x)+z) xor
                          state_i(64*(5*1+x)+z) xor
                          state_i(64*(5*2+x)+z) xor
                          state_i(64*(5*3+x)+z) xor
                          state_i(64*(5*4+x)+z);
            end loop;
        end loop;

        for x in 0 to 4 loop
            for z in 0 to 63 loop
                D(x,z) := C((x-1) mod 5, z) xor C((x+1) mod 5, (z-1) mod 64);
            end loop;
        end loop;

        for x in 0 to 4 loop
            for z in 0 to 63 loop
                for y in 0 to 4 loop
                    state_o(64*(5*y+x)+z) := state_i(64*(5*y+x)+z) xor D(x,z);
                end loop;
            end loop;
        end loop;
        --report "theta " & to_hstring(byte_reverse_f(state_o));
        return state_o;
    end function theta_f;

    function rho_f(state_i : std_logic_vector(0 to 1600 - 1)) return std_logic_vector is
        variable state_o : std_logic_vector(0 to 1600 - 1) := (others => '0');
        type rho_table_t is array(0 to 4, 0 to 4) of natural;
        constant rotate_c : rho_table_t := ((0,36,3,105,210),
                                            (1,300,10,45,66),
                                            (190,6,171,15,253),
                                            (28,55,153,21,120),
                                            (91,276,231,136,78));
    begin
        for x in 0 to 4 loop
            for z in 0 to 63 loop
                for y in 0 to 4 loop
                    state_o(64*(5*y+x)+z) := state_i(64*(5*y+x)+((z-rotate_c(x,y)) mod 64));
                end loop;
            end loop;
        end loop;
        --report "rho " & to_hstring(byte_reverse_f(state_o));
        return state_o;
    end function rho_f;

    function pi_f(state_i : std_logic_vector(0 to 1600 - 1)) return std_logic_vector is
        variable state_o : std_logic_vector(0 to 1600 - 1) := (others => '0');
    begin
        for x in 0 to 4 loop
            for z in 0 to 63 loop
                for y in 0 to 4 loop
                    state_o(64*(5*y+x)+z) := state_i(64*(5*x+((x+3*y) mod 5))+z);
                end loop;
            end loop;
        end loop;
        --report "pi " & to_hstring(byte_reverse_f(state_o));
        return state_o;
    end function pi_f;

    function chi_f(state_i : std_logic_vector(0 to 1600 - 1)) return std_logic_vector is
        variable state_o : std_logic_vector(0 to 1600 - 1) := (others => '0');
    begin
        for x in 0 to 4 loop
            for z in 0 to 63 loop
                for y in 0 to 4 loop
                    state_o(64*(5*y+x)+z) := state_i(64*(5*y+x)+z)
                    xor (   state_i(64*(5*y+((x+2) mod 5))+z) and not state_i(64*(5*y+((x+1) mod 5))+z)   );
                end loop;
            end loop;
        end loop;
        --report "chi " & to_hstring(byte_reverse_f(state_o));
        return state_o;
    end function chi_f;

    function iota_f(state_i : std_logic_vector(0 to 1600 - 1); rnd : integer) return std_logic_vector is
        variable state_o : std_logic_vector(0 to 1600 - 1) := (others => '0');
        type rc_t is array(0 to 23) of std_logic_vector(0 to 63);
        constant round_c : rc_t :=
            (("1000000000000000000000000000000000000000000000000000000000000000"),
            ("0100000100000001000000000000000000000000000000000000000000000000"),
            ("0101000100000001000000000000000000000000000000000000000000000001"),
            ("0000000000000001000000000000000100000000000000000000000000000001"),
            ("1101000100000001000000000000000000000000000000000000000000000000"),
            ("1000000000000000000000000000000100000000000000000000000000000000"),
            ("1000000100000001000000000000000100000000000000000000000000000001"),
            ("1001000000000001000000000000000000000000000000000000000000000001"),
            ("0101000100000000000000000000000000000000000000000000000000000000"),
            ("0001000100000000000000000000000000000000000000000000000000000000"),
            ("1001000000000001000000000000000100000000000000000000000000000000"),
            ("0101000000000000000000000000000100000000000000000000000000000000"),
            ("1101000100000001000000000000000100000000000000000000000000000000"),
            ("1101000100000000000000000000000000000000000000000000000000000001"),
            ("1001000100000001000000000000000000000000000000000000000000000001"),
            ("1100000000000001000000000000000000000000000000000000000000000001"),
            ("0100000000000001000000000000000000000000000000000000000000000001"),
            ("0000000100000000000000000000000000000000000000000000000000000001"),
            ("0101000000000001000000000000000000000000000000000000000000000000"),
            ("0101000000000000000000000000000100000000000000000000000000000001"),
            ("1000000100000001000000000000000100000000000000000000000000000001"),
            ("0000000100000001000000000000000000000000000000000000000000000001"),
            ("1000000000000000000000000000000100000000000000000000000000000000"),
            ("0001000000000001000000000000000100000000000000000000000000000001"));
    begin
        state_o(0 to 63) := state_i(0 to 63) xor round_c(rnd);
        state_o(64 to 1600-1) := state_i(64 to 1600-1);
        --report "iota " & to_hstring(byte_reverse_f(state_o));
        return state_o;
    end function iota_f;

    function keccak_round_f(state_i : std_logic_vector(0 to 1600 - 1); rnd : integer) return std_logic_vector is
        --variable state_o : std_logic_vector(0 to 1600 - 1) := (others => '0');
    begin
        return iota_f(chi_f(pi_f(rho_f(theta_f(state_i)))), rnd);
    end function keccak_round_f;
end package body sha3_pkg;
